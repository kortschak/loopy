// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// dedup-ccs-event breaks gff features from a PB sequencing run and blasr
// alignment passed through the reefer pipeline into
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
	"sort"
	"strings"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
)

var (
	in = flag.String("in", "", "specify input gff file (required)")
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

	names := make(map[string]map[string]struct{})

	sc := featio.NewScanner(gff.NewReader(f))
	for sc.Next() {
		feat := sc.Feat().(*gff.Feature)
		read := feat.FeatAttributes.Get("Read")
		if read == "" {
			continue
		}
		read = strings.Fields(read)[0]
		idx := strings.LastIndex(read, "/")
		e, ok := names[read[:idx]]
		if !ok {
			e = make(map[string]struct{})
			names[read[:idx]] = e
		}
		e[read[idx+1:]] = struct{}{}
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
			s := make([]string, 0, len(coords))
			for c := range coords {
				s = append(s, c)
			}
			sort.Strings(s)
			fmt.Fprintf(nonUnique, "%s\t%v\n", name, s)
		}
	}
}
