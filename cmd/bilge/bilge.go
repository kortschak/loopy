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
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
)

var (
	in     = flag.String("in", "", "specify input gff file (required)")
	thresh = flag.Float64("thresh", 6, "specify minimum total sequence complexity")
	dist   = flag.Bool("dist", false, "only calculate complexity distribution")
	typ    = flag.Int("type", 0, "specify complexity calculation function (0 - WF, 1 - entropic, 2 - Z)")
)

func main() {
	flag.Parse()
	if *in == "" || *typ < 0 || 2 < *typ {
		flag.Usage()
		os.Exit(1)
	}

	cfn := []func(s seq.Sequence, start, end int) (float64, error){
		0: complexity.WF,
		1: complexity.Entropic,
		2: complexity.Z,
	}[*typ]

	f, err := os.Open(*in)
	if err != nil {
		log.Fatalf("failed to open %q: %v", *in, err)
	}
	defer f.Close()

	sc := seqio.NewScanner(fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNAgapped)))
	for sc.Next() {
		seq := sc.Seq().(*linear.Seq)

		// err is always nil for a linear.Seq Start() and End().
		c, _ := cfn(seq, seq.Start(), seq.End())

		if *dist {
			fmt.Printf("%s\t%v\t%d\n", seq.Name(), c, seq.Len())
			continue
		}
		if c >= *thresh {
			fmt.Printf("%60a\n", seq)
		}
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during fasta read: %v", err)
	}
	f.Close()
}
