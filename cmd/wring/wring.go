// Copyright ©2015 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// wring extracts a set of sequences from a SAM file based on a reefer GFF.
package main

import (
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/hts/sam"

	"github.com/kortschak/loopy/blasr"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Fprintln(os.Stderr, "invalid invocation: must have at least one reads file")
		os.Exit(1)
	}

	extract := make(map[string][2]int)
	sc := featio.NewScanner(gff.NewReader(os.Stdin))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		read := f.FeatAttributes.Get("Read")
		if read == "" {
			continue
		}
		fields := strings.Fields(read)
		name := fields[0]
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			log.Fatalf("failed to parse %q: %v", read, err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			log.Fatalf("failed to parse %q: %v", read, err)
		}
		extract[name] = [2]int{start, end}
	}
	err := sc.Error()
	if err != nil {
		log.Fatalf("error during GFF read: %v", err)
	}

	for _, reads := range os.Args[1:] {
		sf, err := os.Open(reads)
		if err != nil {
			log.Fatalf("failed to open %q: %v", reads, err)
		}
		sr, err := sam.NewReader(blasr.NewPBFixReader(sf))
		if err != nil {
			log.Fatalf("failed to open SAM input %q: %v", reads, err)
		}
		for {
			r, err := sr.Read()
			if err != nil {
				if err != io.EOF {
					log.Fatalf("unexpected error reading SAM: %v", err)
				}
				break
			}

			v, ok := extract[r.Name]
			if !ok {
				continue
			}
			// Currently reefer only expects a single hit per read,
			// so any multiples are due to duplicate read file input.
			// Update this behaviour if we change reefer to look at
			// remapping soft-clipped segments.
			delete(extract, r.Name)

			reverse := r.Flags&sam.Reverse != 0
			rng := fmt.Sprintf("//%d_%d", v[0], v[1])
			if reverse {
				rng += "(-)"
				len := r.Seq.Length
				v[0], v[1] = len-v[1], len-v[0]
			}
			v[0] = feat.OneToZero(v[0])
			s := linear.NewSeq(
				r.Name+rng,
				alphabet.BytesToLetters(r.Seq.Expand())[v[0]:v[1]],
				alphabet.DNA,
			)
			if reverse {
				s.Desc = "(sequence revcomp relative to read)"
			}
			fmt.Printf("%60a\n", s)
		}
		sf.Close()
	}
}
