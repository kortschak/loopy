// Copyright ©2015 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// reefer performs blasr alignment and analysis of internal mismatches to
// identify candidate structural variation features.
package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/biogo/biogo/align"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"

	"github.com/kortschak/loopy/blasr"
)

type mat [3]int

var alnmat = mat{1, -1, -1}

func (v *mat) Set(s string) error {
	fields := strings.Split(s, ",")
	if len(fields) != 3 {
		return fmt.Errorf("invalid number of fields: %q", s)
	}
	var err error
	for i, f := range fields {
		v[i], err = strconv.Atoi(f)
		if err != nil {
			return fmt.Errorf("invalid fields: %v", err)
		}
	}
	return nil
}

func (v *mat) String() string { return fmt.Sprintf("%d,%d,%d", v[0], v[1], v[2]) }

var (
	reads       = flag.String("reads", "", "input fasta sequence read file name (required)")
	ref         = flag.String("reference", "", "input reference sequence file name (required)")
	suff        = flag.String("suff", "", "input reference suffix array path")
	useBam      = flag.Bool("bam", false, "use bam file inputs if not running blasr")
	refine      = flag.Bool("refine", true, "use paired SW alignment to refine breakpoints")
	refWindow   = flag.Int("ref-window", 300, "window for refinement around middle of reference indel")
	queryWindow = flag.Int("read-window", 500, "window for refinement beyond ends of of read indel")
	minQueryGap = flag.Int("min-read-gap", 50, "minimum distance between read breakpoints")
	minRefFlank = flag.Int("min-ref-flank", 10, "minimum distance from end of reference window")
	verbose     = flag.Bool("v", false, "verbose logging of breakpoint adjustment")
	blasrPath   = flag.String("blasr", "", "path to blasr if not in $PATH")
	procs       = flag.Int("procs", 1, "number of blasr threads")
	window      = flag.Int("window", 50, "smoothing window")
	minSize     = flag.Int("min", 300, "minimum feature size")
	run         = flag.Bool("run-blasr", true, `actually run blasr
    	false is useful to reconstruct output from fasta input
    	and reefer .blasr outputs`,
	)

	errFile   = flag.String("err", "", "output file name (default to stderr)")
	errStream = os.Stderr
)

func main() {
	flag.Var(&alnmat, "align", "specify the match, mismatch and gap parameters for breakpoint refinement")
	flag.Parse()
	if *reads == "" || (*ref == "" && *run) {
		fmt.Fprintln(os.Stderr, "invalid argument: must have reads, reference and block size set")
		flag.Usage()
		os.Exit(1)
	}

	var err error
	if *errFile != "" {
		errStream, err = os.Create(*errFile)
		if err != nil {
			// Oh, the irony.
			log.Fatalf("failed to create log file: %v", err)
		}
		defer errStream.Close()
		log.SetOutput(errStream)
	}

	// Set up breakpoint refiner.
	var br *refiner
	if *refine {
		refSeq, err := readContigs(*ref)
		if err != nil {
			log.Fatalf("failed to read reference sequences: %v", err)
		}
		br = &refiner{
			refWindow:   *refWindow,
			queryWindow: *queryWindow,
			minQueryGap: *minQueryGap,
			minRefFlank: *minRefFlank,
			ref:         refSeq,
			sw:          makeTable(alnmat),
		}
	}

	out := filepath.Base(*reads)
	f, err := os.Create(out + ".gff")
	if err != nil {
		log.Fatalf("failed to create GFF outfile: %q", out+".gff")
	}
	w := gff.NewWriter(f, 60, true)
	defer f.Close()
	log.Printf("finding alignments for reads in %q", *reads)
	ext := "sam"
	if *useBam && !*run {
		ext = "bam"
	}
	err = deletions(*reads, *ref, *suff, ext, *procs, *run, *window, *minSize, br, w)
	if err != nil {
		log.Fatalf("failed mapping: %v", err)
	}
}

// deletions analyses *sam.Records from mapping reads to the given reference
// using the suffix array file if provided. If run is false, blasr is not
// run and the existing blasr output is used to provide the *sam.Records.
// procs specifies the number of blasr threads to use.
func deletions(reads, ref, suff, ext string, procs int, run bool, window, min int, br *refiner, w *gff.Writer) error {
	base := filepath.Base(reads)
	b := blasr.BLASR{
		Cmd: *blasrPath,

		Reads: reads, Genome: ref, SuffixArray: suff,
		BestN: 1,

		SAM:           true,
		Clipping:      "soft",
		SAMQV:         true,
		CIGARSeqMatch: true,

		Aligned:   base + ".blasr." + ext,
		Unaligned: base + ".blasr.unmapped.fasta",

		Procs: procs,
	}
	if run {
		cmd, err := b.BuildCommand()
		if err != nil {
			return err
		}
		cmd.Stdout = errStream
		cmd.Stderr = errStream
		err = cmd.Run()
		if err != nil {
			return err
		}
	}

	f, err := os.Open(b.Aligned)
	if err != nil {
		return err
	}
	defer f.Close()

	cost := [...]float64{
		sam.CigarInsertion: -2,
		sam.CigarDeletion:  -2,
		sam.CigarEqual:     1,
		sam.CigarMismatch:  -1,

		// Included for explicitness
		sam.CigarSoftClipped: 0,

		// Included to ensure no bounds panic.
		// All CIGAR operations not listed above
		// are given a zero cost.
		sam.CigarBack: 0,
	}

	_, err = w.WriteComment(fmt.Sprintf("smoothing window=%d", window))
	if err != nil {
		return nil
	}
	_, err = w.WriteComment(fmt.Sprintf("minimum feature length=%d", min))
	if err != nil {
		return nil
	}
	gf := &gff.Feature{
		Source:         "reefer",
		Feature:        "discordance",
		FeatFrame:      gff.NoFrame,
		FeatAttributes: gff.Attributes{{Tag: "Read"}, {Tag: "Dup"}},
	}
	var sr interface {
		Read() (*sam.Record, error)
	}
	switch ext {
	case "sam":
		sr, err = sam.NewReader(f)
		if err != nil {
			return err
		}
	case "bam":
		var br *bam.Reader
		br, err = bam.NewReader(f, 0)
		if err != nil {
			return err
		}
		defer br.Close()
		sr = br
	default:
		panic("reefer: invalid extension")
	}
	for {
		r, err := sr.Read()
		if err != nil {
			if err != io.EOF {
				return err
			}
			break
		}

		var (
			scores []costPos
			ref    = r.Start()
			query  int
		)
		for _, co := range r.Cigar {
			for i := 0; i < co.Len(); i++ {
				scores = append(scores, costPos{
					ref:   ref,
					query: query,
					cost:  cost[co.Type()],
				})
				consume := co.Type().Consumes()
				ref += consume.Reference
				query += consume.Query
			}
		}
		if len(scores) <= window {
			continue
		}
		smoothed := make([]costPos, len(scores)-window)
		for i := range scores[:len(scores)-window] {
			smoothed[i] = mean(scores[i : i+window])
		}

		var d deletion
		for i, v := range smoothed[1:] {
			switch {
			case d.record == nil && v.cost < 0 && smoothed[i].cost >= 0:
				d = deletion{record: r, rstart: v.ref + 1, qstart: v.query + 1}
			case d.record != nil && v.cost >= 0 && smoothed[i].cost < 0:
				d.rend = v.ref
				d.qend = v.query
				if d.rend-d.rstart >= min || d.qend-d.qstart >= min {
					gf.SeqName = d.record.Ref.Name()
					gf.FeatStrand = strandFor(d.record)
					if gf.FeatStrand == seq.Minus {
						len := d.record.Seq.Length
						d.qstart, d.qend = len-d.qend, len-d.qstart
					}

					// Adjust ends based on paired SW alignments.
					var refined bool
					d, refined, err = br.adjust(d)
					if err != nil && *verbose {
						log.Printf("failed alignment %s: %v", d.record.Name, err)
					}

					gf.FeatStart = d.rstart
					gf.FeatEnd = d.rend
					if gf.FeatStart == gf.FeatEnd {
						// This is disgusting garbage resulting from
						// GFF not allowing zero length features.
						gf.FeatEnd++
					}

					if refined {
						gf.FeatAttributes = gf.FeatAttributes[:2]
						gf.FeatAttributes[1].Value = strconv.Itoa(d.dup)
					} else {
						gf.FeatAttributes = gf.FeatAttributes[:1]
					}
					gf.FeatAttributes[0].Value = fmt.Sprintf("%s %d %d", d.record.Name, feat.ZeroToOne(d.qstart), d.qend)
					_, err = w.Write(gf)
					if err != nil {
						return err
					}
				}
				d.record = nil
			}
		}
	}
	return nil
}

type deletion struct {
	record *sam.Record

	rstart, rend, dup int
	qstart, qend      int
}

type costPos struct {
	ref, query int
	cost       float64
}

func mean(c []costPos) costPos {
	var mean costPos
	for _, v := range c {
		mean.cost += v.cost
		mean.ref += v.ref
		mean.query += v.query
	}
	scale := float64(len(c))
	mean.cost /= scale
	mean.ref = int(float64(mean.ref)/scale + 0.5)
	mean.query = int(float64(mean.query)/scale + 0.5)
	return mean
}

func strandFor(r *sam.Record) seq.Strand {
	if r.Flags&sam.Reverse != 0 {
		return seq.Minus
	}
	return seq.Plus
}

type refiner struct {
	refWindow   int
	queryWindow int
	minQueryGap int
	minRefFlank int

	ref map[string]*linear.Seq
	sw  align.SW
}

func readContigs(file string) (map[string]*linear.Seq, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	seqs := make(map[string]*linear.Seq)
	sc := seqio.NewScanner(fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNAgapped)))
	for sc.Next() {
		s := sc.Seq().(*linear.Seq)
		seqs[s.ID] = s
	}
	if err != nil {
		return nil, err
	}
	return seqs, nil
}

func makeTable(alnmat mat) align.SW {
	alpha := alphabet.DNAgapped
	match := alnmat[0]
	mismatch := alnmat[1]
	gap := alnmat[2]
	sw := make(align.SW, alpha.Len())
	for i := range sw {
		row := make([]int, alpha.Len())
		for j := range row {
			row[j] = mismatch
		}
		row[i] = match
		sw[i] = row
	}
	for i := range sw {
		sw[0][i] = gap
		sw[i][0] = gap
	}
	return sw
}

// adjustDeletion performs a deletion ends refinement based on a
// pair of Smith-Waterman alignments.
//
//                    l      s   e      r
//  ref:         -----|------+~~~+------|----------
//
//  query_left:  ----|-----------+~~~~~~|~~~~~~+---------------
//                   l           s      m      e
//  query_right: ----------------+~~~~~~|~~~~~~+-----------|---
//                               s      m      e           r
//
//  where ~~ is the region found by CIGAR score walking above in the
//  deletions function.
//
//  align ref(l..r) with query_left(l..m) -> ref(s)-query_left(s)
//  align ref(l..r) with query_right(m..r) -> ref(e)-query_left(e)
//
// This can give either of two outcomes:
//  1. ref(s) < ref(e)
//  2. ref(e) <= ref(s)
//
// The first case is a standard colinear alignment:
//
//                              s   e
//  ref:             -----------+---+-----------------
//                             /     \
//                            /       \
//                           /         \
//                          /           \
//  query: ----------------+-------------+---------------
//                         s             e
//
//
// The second case is a non-colinear alignment:
//
//                              e   s
//  ref:             -----------+---+-----------------
//                               \ /
//                                /
//                               / \
//                              /   \
//                             /     \
//                            /       \
//                           /         \
//                          /           \
//  query: ----------------+-------------+---------------
//                         s             e
//
//
// which has a potential target site duplication interpretation:
//
//                              e   s
//  ref:             -----------+---+-----------------
//                             / \ / \
//                            /   /   \
//                           /   / \   \
//                          /   /   \   \
//                         /   /     \   \
//                        /   /       \   \
//                       /   /         \   \
//                      /   /           \   \
//  query: ------------+---+-------------+---+-----------
//                         s             e
//
// adjustDeletions handles the second case by making ref(s=e) for the
// reference and adding annotation for the length of the duplication
// (d) in ref:
//
//                             s|e s+d
//  ref:             -----------+---+-----------------
//                             / \ / \
//                            /   /   \
//                           /   / \   \
//                          /   /   \   \
//                         /   /     \   \
//                        /   /       \   \
//                       /   /         \   \
//                      /   /           \   \
//  query: ------------+---+-------------+---+-----------
//                    s-d  s             e  e+d
//
func (r *refiner) adjust(d deletion) (refined deletion, ok bool, err error) {
	if r == nil {
		return d, false, nil
	}
	if d.qend-d.qstart < d.rend-d.rstart {
		// Do not do any work for deletions.
		return d, false, fmt.Errorf("not an insertion: len(q)=%d len(r)=%d", d.qend-d.qstart, d.rend-d.rstart)
	}

	name := d.record.Ref.Name()
	ref, ok := r.ref[name]
	if !ok {
		return d, false, fmt.Errorf("no reference sequence for %q", name)
	}

	rs := *ref
	rOff := max(0, d.rstart-r.refWindow/2)
	rs.Seq = ref.Seq[rOff:min(d.rend+r.refWindow/2, len(ref.Seq))]

	q := alphabet.BytesToLetters(d.record.Seq.Expand())

	// Align the left junction of the qeuery to
	// the reference around the indel site.
	qsl := linear.NewSeq(d.record.Name, nil, alphabet.DNAgapped)
	qOffLeft := max(0, d.qstart-r.queryWindow)
	qsl.Seq = q[qOffLeft : (d.qstart+d.qend)/2]
	alnl, err := r.sw.Align(&rs, qsl)
	if err != nil {
		return d, false, err
	}

	// Align the right junction of the qeuery to
	// the reference around the indel site.
	qsr := linear.NewSeq(d.record.Name, nil, alphabet.DNAgapped)
	qOffRight := (d.qstart + d.qend) / 2
	qsr.Seq = q[qOffRight:min(d.qend+r.queryWindow, len(q))]
	alnr, err := r.sw.Align(&rs, qsr)
	if err != nil {
		return d, false, err
	}

	// Get left and right ends of insertion in read
	// and the aligned segment of the reference.
	left := alnl[len(alnl)-1].Features()
	right := alnr[0].Features()

	// Bail out if the alignment extends too far.
	// We might have continued alignment.
	if flank := right[0].Start(); flank < r.minRefFlank {
		return d, false, fmt.Errorf("skipping: right ref flank less than %d from left: len(flank)=%v",
			r.minRefFlank, flank)
	}
	if flank := left[0].End(); len(rs.Seq)-flank < r.minRefFlank {
		return d, false, fmt.Errorf("skipping: left ref flank less than %d from right: len(flank)=%v",
			r.minRefFlank, len(rs.Seq)-flank)
	}

	centrel := r.queryWindow + (d.qend-d.qstart)/2
	centrer := 0

	// Bail out if the insertion is too short.
	// We might have continued alignment.
	if gap := centrel - left[1].End(); gap < r.minQueryGap {
		return d, false, fmt.Errorf("skipping left: left query gap less than %d from centre: len(gap)=%v",
			r.minQueryGap, gap)
	}
	if gap := right[1].Start() - centrer; gap < r.minQueryGap {
		return d, false, fmt.Errorf("skipping right: right query gap less than %d from centre: len(gap)=%v",
			r.minQueryGap, gap)
	}

	d.rstart = rOff + left[0].End()
	d.rend = rOff + right[0].Start()
	if d.rend <= d.rstart {
		d.dup = d.rstart - d.rend
		d.rstart = d.rend
	}

	d.qstart = qOffLeft + left[1].End()
	d.qend = qOffRight + alnr[0].Features()[1].Start()

	return d, true, nil
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
