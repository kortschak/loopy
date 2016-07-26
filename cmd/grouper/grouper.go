// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// grouper reports the genomic extent of a group of reefer features
// where the group has been identified by press or press-global.
package main

import (
	"fmt"
	"log"
	"os"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
)

func main() {
	groups := make(map[string]struct {
		chrom      string
		start, end int
	})

	sc := featio.NewScanner(gff.NewReader(os.Stdin))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		g := f.FeatAttributes.Get("Group")
		if g == "" {
			continue
		}
		grp, ok := groups[g]
		if !ok {
			groups[g] = struct {
				chrom      string
				start, end int
			}{chrom: f.SeqName, start: f.FeatStart, end: f.FeatEnd}
			continue
		}
		if f.FeatStart < grp.start {
			grp.start = f.FeatStart
		}
		if grp.end < f.FeatEnd {
			grp.end = f.FeatEnd
		}
		groups[g] = grp
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during gff read: %v", err)
	}

	for k, v := range groups {
		fmt.Printf("%s\t%d\t%d\t%s\n", v.chrom, v.start, v.end, k)
	}
}
