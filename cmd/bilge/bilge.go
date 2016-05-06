// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// bilge filters a set of sequences for low complexity.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/complexity"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

var (
	in     = flag.String("in", "", "specify input gff file (required)")
	thresh = flag.Float64("thresh", 6, "specify minimum total sequence complexity")
	dist   = flag.Bool("dist", false, "only calculate complexity distribution")
)

func main() {
	flag.Parse()
	if *in == "" {
		flag.Usage()
		os.Exit(1)
	}

	f, err := os.Open(*in)
	if err != nil {
		log.Fatalf("failed to open %q: %v", *in, err)
	}
	defer f.Close()

	sc := seqio.NewScanner(fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNAgapped)))
	for sc.Next() {
		seq := sc.Seq().(*linear.Seq)

		// err is always nil for a linear.Seq Start() and End().
		cz, _ := complexity.WF(seq, seq.Start(), seq.End())

		if *dist {
			fmt.Printf("%s\t%v\n", seq.Name(), cz)
			continue
		}
		if cz >= *thresh {
			fmt.Printf("%60a\n", seq)
		}
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during fasta read: %v", err)
	}
	f.Close()
}
