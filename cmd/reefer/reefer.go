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

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/hts/sam"

	"github.com/kortschak/loopy/blasr"
)

var (
	reads     = flag.String("reads", "", "input fasta sequence read file name (required)")
	ref       = flag.String("reference", "", "input reference sequence file name (required)")
	suff      = flag.String("suff", "", "input reference suffix array path")
	blasrPath = flag.String("blasr", "", "path to blasr if not in $PATH")
	procs     = flag.Int("procs", 1, "number of blasr threads")
	window    = flag.Int("window", 50, "smoothing window")
	min       = flag.Int("min", 300, "minimum feature size")
	run       = flag.Bool("run-blasr", true, `actually run blasr
    	false is useful to reconstruct output from fasta input
    	and reefer .blasr outputs`,
	)

	discords = flag.Bool("discords", false, "output GFF file of discordant features")
	outFile  = flag.String("out", "", "output file name (default to stdout)")
	errFile  = flag.String("err", "", "output file name (default to stderr)")
)

var (
	outStream = os.Stdout
	errStream = os.Stderr
)

func main() {
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
	if *outFile != "" {
		outStream, err = os.Create(*outFile)
		if err != nil {
			log.Fatalf("failed to create out file: %v", err)
		}
		defer outStream.Close()
	}

	var w *gff.Writer
	if *discords {
		out := filepath.Base(*reads)
		f, err := os.Create(out + ".gff")
		if err != nil {
			log.Fatalf("failed to create GFF outfile: %q", out+".gff")
		}
		w = gff.NewWriter(f, 60, true)
		defer f.Close()
	}
	log.Printf("finding alignments for reads in %q", *reads)
	err = deletions(*reads, *ref, *suff, *procs, *run, *window, *min, w)
	if err != nil {
		log.Fatalf("failed mapping: %v", err)
	}
}

// deletions returns a []*sam.Record from mapping reads to the given reference
// using the suffix array file if provided. If run is false, blasr is not
// run and the existing blasr output is used to reconstruct the []*sam.Record.
// procs specifies the number of blasr threads to use.
func deletions(reads, ref, suff string, procs int, run bool, window, min int, w *gff.Writer) error {
	base := filepath.Base(reads)
	b := blasr.BLASR{
		Cmd: *blasrPath,

		Reads: reads, Genome: ref, SuffixArray: suff,
		BestN: 1,

		SAM:           true,
		Clipping:      "soft",
		SAMQV:         true,
		CIGARSeqMatch: true,

		Aligned:   base + ".blasr.sam",
		Unaligned: base + ".blasr.unmapped.sam",

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
		FeatAttributes: gff.Attributes{{Tag: "Read"}},
	}

	sr, err := sam.NewReader(f)
	if err != nil {
		return err
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

		type deletion struct {
			record *sam.Record

			rstart, rend int
			qstart, qend int
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
					gf.FeatStart = d.rstart
					gf.FeatEnd = d.rend
					gf.FeatAttributes[0].Value = fmt.Sprintf("%s %d %d", d.record.Name, feat.ZeroToOne(d.qstart), d.qend)
					_, err = w.Write(gf)
					if err != nil {
						return nil
					}
				}
				d.record = nil
			}
		}
	}
	return nil
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
