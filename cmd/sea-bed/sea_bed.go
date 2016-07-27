// Copyright ©2016 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// sea-bed outputs a set of fasta sequences based on a reference and
// set of bed files.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/bed"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

var (
	ref   = flag.String("ref", "", "genome fasta file")
	flank = flag.Int("flank", 0, "genome fasta file")
)

func main() {
	flag.Parse()
	if flag.NArg() == 0 {
		fmt.Fprintln(os.Stderr, "need at least one bed3 file input")
		os.Exit(0)
	}
	if *ref == "" {
		flag.Usage()
		os.Exit(0)
	}

	seqs, err := readContigs(*ref)
	if err != nil {
		log.Fatalf("failed to read reference file: %v", err)
	}

	for _, in := range flag.Args() {
		bf, err := os.Open(in)
		if err != nil {
			log.Fatalf("failed to open bed file: %v", err)
		}

		br, err := bed.NewReader(bf, 3)
		if err != nil {
			log.Fatalf("failed to read bed file: %v", err)
		}

		out, err := os.Create(basename(in) + ".mfa")
		if err != nil {
			log.Fatalf("failed to create fasta file: %v", err)
		}

		sc := featio.NewScanner(br)
		for sc.Next() {
			f := sc.Feat().(*bed.Bed3)
			s := *seqs[f.Chrom]
			start := max(0, f.ChromStart-*flank)
			end := min(f.ChromEnd+*flank, len(s.Seq))
			s.Seq = s.Seq[start:end]
			s.ID = fmt.Sprintf("%s[%d,%d)", s.ID, start, end)
			if *flank != 0 {
				s.Desc = fmt.Sprintf("flanking [%d,%d)", f.ChromStart, f.ChromEnd)
			}
			_, err := fmt.Fprintf(out, "%60a\n", &s)
			if err != nil {
				log.Fatalf("failed to write fasta sequence: %v", err)
			}
		}
		err = sc.Error()
		if err != nil {
			log.Fatalf("failed to read bed file: %v", err)
		}
		out.Close()
		bf.Close()
	}
}

func readContigs(file string) (map[string]*linear.Seq, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	seqs := make(map[string]*linear.Seq)
	sc := seqio.NewScanner(fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNA)))
	for sc.Next() {
		s := sc.Seq().(*linear.Seq)
		seqs[s.ID] = s
	}
	if err != nil {
		log.Fatalf("error during fasta read: %v", err)
	}
	return seqs, nil
}

func basename(path string) string {
	path = filepath.Base(path)
	ext := filepath.Ext(path)
	return path[:len(path)-len(ext)]
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
