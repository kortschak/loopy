// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// dedup-ccs breaks fasta sequences from a PB sequencing run into
// uniquely identified and non-uniquely identified lists.
//
// uniquely - not CCS reads
// non-uniqu - CCS reads
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

var (
	in = flag.String("in", "", "specify input fasta file (required)")
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

	names := make(map[string][]string)

	sc := seqio.NewScanner(fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNAgapped)))
	for sc.Next() {
		seq := sc.Seq().(*linear.Seq)
		idx := strings.LastIndex(seq.ID, "/")
		names[seq.ID[:idx]] = append(names[seq.ID[:idx]], seq.ID[idx+1:])
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during fasta read: %v", err)
	}
	f.Close()

	base := filepath.Base(*in)
	unique, err := os.Create(base + ".unique.text")
	if err != nil {
		log.Fatalf("failed to create %q: %v", base+".unique.text", err)
	}
	defer unique.Close()
	nonUnique, err := os.Create(base + ".non-unique.text")
	if err != nil {
		log.Fatalf("failed to create %q: %v", base+".non-unique.text", err)
	}
	defer nonUnique.Close()
	for name, coords := range names {
		switch len(coords) {
		case 0:
		case 1:
			fmt.Fprintln(unique, name)
		default:
			fmt.Fprintf(nonUnique, "%s\t%v\n", name, coords)
		}
	}
}
