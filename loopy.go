// Copyright ©2015 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// loopy performs blasr alignment and unmapped flank remapping to identify candidate
// structural variation features.
//
// The program is based on the original python code by Steve Turner.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
)

var (
	reads     = flag.String("reads", "", "input fasta sequence read file name (required)")
	ref       = flag.String("reference", "", "input reference sequence file name (required)")
	suff      = flag.String("suff", "", "input reference suffix array path")
	blasrPath = flag.String("blasr", "", "path to blasr if not in $PATH")
	procs     = flag.Int("procs", 1, "number of blasr threads")
	flank     = flag.Int("flank", 50, "minimum flank length")
	length    = flag.Int("length", 200, "minimum blasr search alignment length")
	discords  = flag.Bool("discords", false, "output GFF file of discordant features")
	run       = flag.Bool("run-blasr", true, `actually run blasr
    	false is useful to reconstruct output from fasta input
    	and loopy .blasr outputs`,
	)

	outFile = flag.String("out", "", "output file name (default to stdout)")
	errFile = flag.String("err", "", "output file name (default to stderr)")
)

func main() {
	flag.Parse()
	if *reads == "" || *ref == "" {
		fmt.Fprintln(os.Stderr, "invalid argument: must have reads, reference and block size set")
		flag.Usage()
		os.Exit(1)
	}

	if *errFile != "" {
		w, err := os.Create(*errFile)
		if err != nil {
			// Oh, the irony.
			log.Fatalf("failed to create log file: %v", err)
		}
		defer w.Close()
		log.SetOutput(w)
	}
	outStream := os.Stdout
	if *outFile != "" {
		outStream, err := os.Create(*errFile)
		if err != nil {
			log.Fatalf("failed to create out file: %v", err)
		}
		defer outStream.Close()
	}

	log.Printf("finding flanks of reads in %q", *reads)
	core, err := hitSetFrom(*reads, *ref, *suff, *procs, *run)
	if err != nil {
		log.Fatalf("failed initial mapping: %v", err)
	}

	// Prepare flank sequences and remap them.
	out := filepath.Base(*reads)
	leftSeqs := out + ".left.in.fa"
	rightSeqs := out + ".right.in.fa"

	log.Printf("writing flanks to %q and %q", leftSeqs, rightSeqs)
	err = writeFlankSeqs(*reads, core, *flank, leftSeqs, rightSeqs)
	if err != nil {
		log.Fatalf("failed to write flanks: %v", err)
	}

	log.Printf("remapping left flanks of reads from %q", leftSeqs)
	left, err := hitSetFrom(leftSeqs, *ref, *suff, *procs, *run)
	if err != nil {
		log.Fatalf("failed left flank remapping: %v", err)
	}

	log.Printf("remapping right flanks of reads from %q", rightSeqs)
	right, err := hitSetFrom(rightSeqs, *ref, *suff, *procs, *run)
	if err != nil {
		log.Fatalf("failed right flank remapping: %v", err)
	}

	var w *gff.Writer
	if *discords {
		f, err := os.Create(out + ".gff")
		if err != nil {
			log.Fatalf("failed to create GFF outfile: %q", out+".gff")
		}
		w = gff.NewWriter(f, 60, true)
		defer f.Close()
	}
	err = writeResults(core, left, right, outStream, *length, *flank, w)
	if err != nil {
		log.Fatalf("failed to write results: %v", err)
	}
}

// hitSet represents a collection of blasr mapping results.
type hitSet map[string]*blasrHit

// hitSetFrom returns a hitSet from mapping reads to the given reference
// using the suffix array file if provided. If run is false, blasr is not
// run and the existing blasr output is used to reconstruct the hitSet.
// procs specifies the number of blasr threads to use.
func hitSetFrom(reads, ref, suff string, procs int, run bool) (hitSet, error) {
	base := filepath.Base(reads)
	b := BLASR{
		Cmd: *blasrPath,

		Reads: reads, Genome: ref, SuffixArray: suff,
		BestN: 1, Format: 4,

		Aligned:   base + ".blasr",
		Unaligned: base + ".blasr.unmapped",

		Procs: procs,
	}
	if run {
		cmd, err := b.BuildCommand()
		if err != nil {
			return nil, err
		}
		cmd.Stdout = os.Stderr
		cmd.Stderr = os.Stderr
		err = cmd.Run()
		if err != nil {
			return nil, err
		}
	}

	f, err := os.Open(b.Aligned)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	hits := make(hitSet)
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		b, err := newBlasrHit(sc.Text())
		if err != nil {
			return nil, err
		}
		hits[b.qName] = b
	}

	return hits, sc.Err()
}

// writeFlankSeqs writes fasta files containing the sequence of unmapped flanks
// identified in the primary hits provided. cutoff specifies the minimum sequence
// length to consider. left and right specify the filenames for the left and right
// flank fasta sequence files.
func writeFlankSeqs(reads string, hits hitSet, cutoff int, left, right string) error {
	f, err := os.Open(reads)
	if err != nil {
		return err
	}
	defer f.Close()

	lf, err := os.Create(left)
	if err != nil {
		return err
	}
	rf, err := os.Create(right)
	if err != nil {
		return err
	}

	r := fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNA))
	sc := seqio.NewScanner(r)
	for sc.Next() {
		seq := sc.Seq().(*linear.Seq)
		h, ok := hits[seq.Name()]
		if !ok {
			continue
		}

		all := seq.Seq
		if h.qStart >= cutoff {
			seq.Seq = all[:h.qStart]
			_, err := fmt.Fprintf(lf, "%60a\n", seq)
			if err != nil {
				return err
			}
		}
		if h.qLen-h.qEnd >= cutoff {
			seq.Seq = all[h.qEnd:]
			_, err := fmt.Fprintf(rf, "%60a\n", seq)
			if err != nil {
				return err
			}
		}
	}
	err = sc.Error()
	if err != nil {
		return err
	}
	err = lf.Close()
	if err != nil {
		return err
	}
	return rf.Close()
}

// writeResults writes out the results of the analysis in a format similar to the
// Pacific Biosciences bridgemapper program (29 tab separated fields). It also writes
// candidate discordances to the discords gff.Writer if it is not nil. Flanks less than
// flank long are not considered and primay mappings less than length long are omitted.
func writeResults(core, left, right hitSet, out io.Writer, length, flank int, discords *gff.Writer) error {
	for id, c := range core {
		if c.qEnd-c.qStart < length {
			continue
		}
		l, ok := left[id]
		if ok && abs(l.tEnd-l.tStart) < flank {
			l = nil
		}
		r, ok := right[id]
		if ok && abs(r.tEnd-r.tStart) < flank {
			r = nil
		}
		if l == nil && r == nil {
			continue
		}
		_, err := fmt.Fprintf(out, "%s\t%d\t%v\t%v\t%v\n", id, c.qLen, l, c, r)
		if err != nil {
			return err
		}
		if discords != nil {
			for _, f := range [2]*blasrHit{l, r} {
				if f == nil {
					continue
				}
				if f.tName != c.tName {
					_, err = discords.Write(&gff.Feature{
						SeqName:    f.tName,
						Feature:    "flank",
						Source:     "loopy",
						FeatStart:  f.tStart,
						FeatEnd:    f.tEnd,
						FeatScore:  floatPtr(float64(f.score)),
						FeatStrand: f.qStrand,
						FeatFrame:  gff.NoFrame,
					})
					if err != nil {
						return err
					}
				} else if f.tStrand == c.tStrand {
					for _, g := range gapOrOverlap(f, c, flank) {
						_, err = discords.Write(g)
						if err != nil {
							return err
						}
					}
				}
			}
		}
	}
	return nil
}

func abs(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

func floatPtr(f float64) *float64 {
	return &f
}

// gapOrOverlap returns features that describe insertion or deletion events
// in the reads relative to the reference. Only features cutoff or longer are
// returned and pairs of read insertion/reference deletion that are within
// cutoff in length are discarded.
func gapOrOverlap(flank, core *blasrHit, cutoff int) []*gff.Feature {
	if flank.tName != core.tName {
		panic("bad hit pair")
	}

	var (
		qGapStart, qGapEnd int
		tGapStart, tGapEnd int
	)
	if flank.qStart < core.qStart {
		qGapStart = flank.qEnd
		qGapEnd = core.qStart

		tGapStart = flank.tEnd
		tGapEnd = core.tStart
	} else {
		qGapStart = core.qEnd
		qGapEnd = core.qEnd + flank.qStart

		tGapStart = core.tEnd
		tGapEnd = flank.tStart
	}
	if tGapEnd < tGapStart {
		tGapEnd, tGapStart = tGapStart, tGapEnd
	}

	if abs((qGapEnd-qGapStart)-(tGapEnd-tGapStart)) < cutoff {
		return nil
	}

	f := make([]*gff.Feature, 0, 2)
	if qGapEnd-qGapStart >= cutoff {
		f = append(f, &gff.Feature{
			SeqName:    flank.tName,
			Feature:    "insertion",
			Source:     "loopy",
			FeatStart:  tGapStart,
			FeatEnd:    tGapEnd,
			FeatStrand: flank.qStrand,
			FeatFrame:  gff.NoFrame,
			FeatAttributes: gff.Attributes{{
				Tag:   "Query",
				Value: fmt.Sprintf("%s %d %d", flank.qName, qGapStart, qGapEnd),
			}},
		})
	}
	if tGapEnd-tGapStart >= cutoff {
		f = append(f, &gff.Feature{
			SeqName:    flank.tName,
			Feature:    "deletion",
			Source:     "loopy",
			FeatStart:  tGapStart,
			FeatEnd:    tGapEnd,
			FeatStrand: flank.qStrand,
			FeatFrame:  gff.NoFrame,
		})
	}
	return f
}

const (
	qnameField = iota
	tnameField
	scoreField
	pctsimilarityField
	qstrandField
	qstartField
	qendField
	qseqlengthField
	tstrandField
	tstartField
	tendField
	tseqlengthField
	mapqvField
	ncellsField
	clusterScoreField
	probscoreField
	numSigClustersField

	numFields
)

// blasrHits is a blasr mapping event.
type blasrHit struct {
	qName   string
	qStrand seq.Strand
	qStart  int
	qEnd    int
	qLen    int

	tName   string
	tStrand seq.Strand
	tStart  int
	tEnd    int
	tLen    int

	score      int
	similarity float64
	mapQV      int
}

func handlePanic(err *error) {
	r := recover()
	if r != nil {
		switch r := r.(type) {
		case error:
			*err = r
		default:
			panic(r)
		}
	}
}

// newBlasrHit returns a blasrHit parsed from a blasr format 4 line.
func newBlasrHit(line string) (b *blasrHit, err error) {
	defer handlePanic(&err)
	fields := strings.Fields(line)
	return &blasrHit{
		// The original code strips the subread start and end from the qname.
		// This is incorrect since multiple movies may exists in the read file,
		// resulting in clobbered map entries (this is also true in the
		// original python).
		// The consequence of this may be miscalculation of query start, end
		// and length values resulting in index out of range or silent sequence
		// truncation.
		// The alternative is to group by read, but I can't see the benefit of
		// that here.
		qName: fields[qnameField],

		qStrand: mustStrand(mustAtoi(fields[qstrandField])),
		qStart:  mustAtoi(fields[qstartField]),
		qEnd:    mustAtoi(fields[qendField]),
		qLen:    mustAtoi(fields[qseqlengthField]),

		tName:   fields[tnameField],
		tStrand: mustStrand(mustAtoi(fields[tstrandField])),
		tStart:  mustAtoi(fields[tstartField]),
		tEnd:    mustAtoi(fields[tendField]),
		tLen:    mustAtoi(fields[tseqlengthField]),

		score:      mustAtoi(fields[scoreField]),
		similarity: mustAtof(fields[pctsimilarityField]),
		mapQV:      mustAtoi(fields[mapqvField]),
	}, nil
}

func mustAtoi(s string) int {
	i, err := strconv.Atoi(s)
	if err != nil {
		panic(err)
	}
	return i
}

func mustAtof(s string) float64 {
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		panic(err)
	}
	return f
}

func mustStrand(s int) seq.Strand {
	switch s {
	case 0:
		return seq.Minus
	case 1:
		return seq.Plus
	default:
		panic(fmt.Sprintf("bad strand value: %d", s))
	}
}

func (b *blasrHit) String() string {
	const empty = "_\t_\t_\t_\t_\t_\t_\t_\t_"
	if b == nil {
		return empty
	}

	start := b.tStart
	end := b.tEnd
	if b.tStrand == 1 {
		start = b.tLen - start
		end = b.tLen - end
	}
	return fmt.Sprintf("%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d",
		b.qStart,
		b.qEnd,
		b.tName,
		b.tStrand,
		start,
		end,
		b.score,
		b.similarity,
		b.mapQV,
	)
}
