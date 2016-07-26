// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// keelhaul drops fasta sequences from stdin containing IDs in
// the exclude parameter file.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

var exclude = flag.String("exclude", "", "specify file containing excluded reads")

func main() {
	flag.Parse()
	if *exclude == "" {
		flag.Usage()
		os.Exit(1)
	}

	nameSet := make(map[string]struct{})
	f, err := os.Open(*exclude)
	if err != nil {
		log.Fatalf("failed to open exclude file %q: %v", *exclude, err)
	}
	ls := bufio.NewScanner(f)
	for ls.Scan() {
		nameSet[ls.Text()] = struct{}{}
	}
	err = ls.Err()
	if err != nil {
		log.Fatalf("failed to read exclude file: %v", err)
	}

	sc := seqio.NewScanner(fasta.NewReader(os.Stdin, linear.NewSeq("", nil, alphabet.DNA)))
	for sc.Next() {
		s := sc.Seq().(*linear.Seq)
		if _, ok := nameSet[s.ID]; ok {
			continue
		}
		fmt.Printf("%60a\n", s)
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during gff read: %v", err)
	}
}
