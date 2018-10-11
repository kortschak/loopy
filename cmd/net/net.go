// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// net performs set operation on reefer pressed events. Input gff feature score
// field must be either not set or set by previous use of net. The coordinate
// systems used for the different inputs is expected to be the same.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
)

var (
	left   = flag.String("a", "", "specify left gff file (required)")
	right  = flag.String("b", "", "specify right gff file (required)")
	thresh = flag.Float64("thresh", 0.90, "specify minumum jaccard similarity for identity between events - must be >= value used by press")
	op     = flag.String("op", "sub", `specify set operation (from "sub" (a\b), "union" (a∪b), "intersect" (a∩b)`)
)

func main() {
	flag.Parse()
	if *left == "" || *right == "" || !validOp(*op) {
		flag.Usage()
		os.Exit(1)
	}

	a, err := events(*left)
	if err != nil {
		log.Fatal(err)
	}
	b, err := events(*right)
	if err != nil {
		log.Fatal(err)
	}

	var c []*gff.Feature
	switch *op {
	case "sub":
		c = sub(a, b, *thresh)
	case "union":
		c = union(a, b, *thresh)
	case "intersect":
		c = intersect(a, b, *thresh)
	}
	w := gff.NewWriter(os.Stdout, 60, true)
	for _, v := range c {
		w.Write(v)
	}
}

func validOp(op string) bool {
	return op == "sub" || op == "union" || op == "intersect"
}

// events returns the maximally extended events from the press gff file given.
func events(file string) (map[int]*gff.Feature, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, fmt.Errorf("failed to open %q: %v", file, err)
	}
	defer f.Close()
	set := make(map[int]*gff.Feature)
	sc := featio.NewScanner(gff.NewReader(f))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		r := strings.TrimRightFunc(f.FeatAttributes.Get("Repeat"), func(r rune) bool {
			return r == ' ' || ('0' <= r && r <= '9')
		})
		g := f.FeatAttributes.Get("Group")
		gid, err := strconv.Atoi(g)
		if err != nil {
			return nil, fmt.Errorf("failed to parse group ID: %v", err)
		}
		p, ok := set[gid]
		if !ok {
			f.FeatAttributes = gff.Attributes{
				{Tag: "Group", Value: g},
				{Tag: "Repeat", Value: r},
			}
			if f.FeatScore == nil {
				f.FeatScore = new(float64)
				*(f.FeatScore) = 1
			}
			set[gid] = f
			continue
		}
		*(p.FeatScore)++
		if f.FeatStart < p.FeatStart {
			p.FeatStart = f.FeatStart
		}
		if f.FeatEnd > p.FeatEnd {
			p.FeatEnd = f.FeatEnd
		}
	}
	if err := sc.Error(); err != nil {
		return nil, fmt.Errorf("error during gff read: %v", err)
	}
	return set, nil
}

// sub returns the result of the set operation a\b. It does this using the
// naive O(n^2) approach rather than using a collection of interval trees
// since len(a) and len(b) are small.
func sub(a, b map[int]*gff.Feature, thresh float64) []*gff.Feature {
	for ka, ea := range a {
		for _, eb := range b {
			if jaccard(ea, eb) >= thresh {
				delete(a, ka)
				break
			}
		}
	}
	c := make([]*gff.Feature, 0, len(a))
	for _, e := range a {
		c = append(c, e)
	}
	return c
}

// union returns the result of the set operation a∪b. It does this using the
// naive O(n^2) approach rather than using a collection of interval trees
// since len(a) and len(b) are small.
func union(a, b map[int]*gff.Feature, thresh float64) []*gff.Feature {
	for ka, ea := range a {
		if ka < 0 {
			// Ignore newly added events from b.
			continue
		}
		for kb, eb := range b {
			if jaccard(ea, eb) >= thresh {
				a[ka].FeatAttributes = gff.Attributes{
					{Tag: "GroupA", Value: fmt.Sprint(ka)},
					{Tag: "GroupB", Value: fmt.Sprint(kb)},
				}
			} else {
				a[ka].FeatAttributes = gff.Attributes{
					{Tag: "GroupA", Value: fmt.Sprint(ka)},
				}
				eb.FeatAttributes = gff.Attributes{
					{Tag: "GroupB", Value: fmt.Sprint(kb)},
				}
				a[-kb] = eb
			}
		}
	}
	c := make([]*gff.Feature, 0, len(a))
	for _, e := range a {
		c = append(c, e)
	}
	return c
}

// intersect returns the result of the set operation a∩b. It does this using the
// naive O(n^2) approach rather than using a collection of interval trees
// since len(a) and len(b) are small.
func intersect(a, b map[int]*gff.Feature, thresh float64) []*gff.Feature {
	var c []*gff.Feature
	for ka, ea := range a {
		for kb, eb := range b {
			if jaccard(ea, eb) >= thresh {
				r := strings.TrimRightFunc(ea.FeatAttributes.Get("Repeat"), func(r rune) bool {
					return r == ' ' || ('0' <= r && r <= '9')
				})
				ea.FeatAttributes = gff.Attributes{
					{Tag: "Group", Value: fmt.Sprint(ka)},
					{Tag: "GroupOther", Value: fmt.Sprint(kb)},
					{Tag: "Repeat", Value: r},
				}
				c = append(c, ea)
			}
		}
	}
	return c
}

func jaccard(a, b *gff.Feature) float64 {
	n := intersection(a, b)
	return float64(n) / (float64(a.Len() + b.Len() - n))
}

func intersection(a, b *gff.Feature) int {
	if a.SeqName != b.SeqName {
		return 0
	}
	return max(0, min(a.FeatEnd, b.FeatEnd)-max(a.FeatStart, b.FeatStart))
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
