// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// ranks reports the repeat types and counts for each group from a GFF on stdin.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
)

var doGrouping = flag.Bool("group", false, "output grouped counts")

func main() {
	flag.Parse()

	var grps []map[string]int
	sc := featio.NewScanner(gff.NewReader(os.Stdin))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		r := f.FeatAttributes.Get("Repeat")
		g := f.FeatAttributes.Get("Group")
		typ := strings.Fields(r)[0]
		if !*doGrouping {
			fmt.Printf("%s\t%s\n", g, typ)
		}
		gid, err := strconv.Atoi(g)
		if err != nil {
			log.Fatalf("failed to parse group id: %v", err)
		}
		grps = add(grps, gid, typ)
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during gff read: %v", err)
	}

	if !*doGrouping {
		return
	}
	for gid, g := range grps {
		if g == nil {
			continue
		}
		fmt.Printf("%d\t", gid)
		m := sortedMap(g)
		for i, t := range m {
			if i != 0 {
				fmt.Print(" ")
			}
			fmt.Printf("%s:%d", t.typ, t.n)
		}
		name := nameHeuristic(m)
		fmt.Printf("\t%s\t%s\n", name, trunc(name, 5))
	}
}

func add(grps []map[string]int, gid int, typ string) []map[string]int {
	switch {
	case gid == len(grps):
		grps = append(grps, make(map[string]int))
	case gid > len(grps) && gid < cap(grps):
		grps = grps[:gid+1]
	case gid > len(grps):
		t := make([]map[string]int, gid+1)
		copy(t, grps)
		grps = t
	}
	if grps[gid] == nil {
		grps[gid] = make(map[string]int)
	}
	grps[gid][typ]++
	return grps
}

type mapElement struct {
	typ string
	n   int
}

type byCount []mapElement

func (m byCount) Len() int { return len(m) }
func (m byCount) Less(i, j int) bool {
	if m[i].n < m[j].n {
		return true
	}
	// Heuristic for sort that longer names are likely to be
	// a tighter definition, so use them in preference.
	return m[i].n == m[j].n && len(m[i].typ) < len(m[j].typ)
}
func (m byCount) Swap(i, j int) { m[i], m[j] = m[j], m[i] }

func sortedMap(g map[string]int) []mapElement {
	m := make([]mapElement, 0, len(g))
	for typ, n := range g {
		m = append(m, mapElement{typ: typ, n: n})
	}
	sort.Sort(sort.Reverse(byCount(m)))
	return m
}

func nameHeuristic(g []mapElement) string {
	if len(g) == 0 {
		return ""
	}

	// Majority rule.
	var n int
	for _, e := range g {
		n += e.n
	}
	r := float64(g[0].n) / float64(n)
	if r > 0.5 || (r == 0.5 && len(g) > 2) {
		return g[0].typ
	}

	// Alu heuristic.
	if isAlu(g[0].typ) {
		return trunc(g[0].typ, 5)
	}

	// Fusion.
	var names []string
	for _, t := range g {
		names = append(names, t.typ)
	}
	return strings.Join(names, "/")
}

func isAlu(t string) bool {
	return strings.HasPrefix(strings.ToLower(t), "alu")
}

func trunc(name string, n int) string {
	return name[:min(5, len(name))]
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
