// Copyright ©2015 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// mangle does name mangling on a multiple fasta sequence file.
// It replaces the fasta ID with the sha1 of the fasta descline,
// failing if there is a sha1 collision. mangle is required for
// censor analysis of sequences with long fasta IDs (~80 columns).
package main

import (
	"bufio"
	"crypto/sha1"
	"flag"
	"fmt"
	"log"
	"os"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

var apply = flag.String("unmangle", "", "apply the inverse name mangling to the specified map file")

func main() {
	flag.Parse()
	if *apply != "" {
		unmangle(*apply)
		return
	}
	mangle()
}

func mangle() {
	seen := make(map[string]bool)
	hash := sha1.New()
	sc := seqio.NewScanner(fasta.NewReader(os.Stdin, linear.NewSeq("", nil, alphabet.DNA)))
	for sc.Next() {
		s := sc.Seq().(*linear.Seq)
		if s.Desc == "" {
			s.Desc = s.ID
		} else {
			s.Desc = fmt.Sprintf("%s %s", s.ID, s.Desc)
		}
		hash.Write([]byte(s.Desc))
		s.ID = fmt.Sprintf("%040x", hash.Sum(nil))
		if seen[s.ID] {
			log.Fatalf("duplicate sha1: %s", s.ID)
		}
		seen[s.ID] = true
		hash.Reset()
		fmt.Printf("%60a\n", s)
	}
}

const (
	queryNameField = iota

	_ // queryStartField
	_ // queryEndField
	_ // repeatTypeField
	_ // repeatStartField
	_ // repeatEndField
	_ // strandField
	_ // alignment similarity
	_ // alignment positive fraction
	_ // scoreField
	_ // query coverage fraction - not used because we need class name anyway.
	_ // repeat coverage fraction - not used because we need class name anyway.

	numberOfFields
)

func unmangle(mapfile string) {
	table := make(map[string]string)
	sc := seqio.NewScanner(fasta.NewReader(os.Stdin, linear.NewSeq("", nil, alphabet.DNA)))
	for sc.Next() {
		s := sc.Seq().(*linear.Seq)
		id := strings.Fields(s.Desc)[0]
		if id == "" {
			log.Fatalf("no id for sequence %s", s.ID)
		}
		table[s.ID] = id
	}

	f, err := os.Open(mapfile)
	if err != nil {
		log.Fatalf("failed to open map file %q: %v", mapfile, err)
	}
	s := bufio.NewScanner(f)
	for s.Scan() {
		line := s.Text()
		fields := strings.Fields(line)
		if len(fields) != numberOfFields {
			log.Fatalf("unexpected number of fields in line %q", line)
		}
		id := table[fields[0]]
		if id == "" {
			log.Fatalf("no id for map query %s", fields[0])
		}
		fields[0] = id
		for i, f := range fields {
			if i != 0 {
				fmt.Print("\t")
			}
			fmt.Print(f)
		}
		fmt.Println()
	}
}
