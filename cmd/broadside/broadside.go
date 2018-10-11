// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// broadside reports the repeat types and counts for each group of events in a trio
// of individuals from a press GFF on stdin. Each of the individuals must be aligned
// to the same reference.
//
// The use of this program makes most sense when the input GFF stream is collection of
// features that are in fil indivdual, but not in the pat or mat individuals. This
// operation can be performed using the net command.
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
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

var (
	fil = flag.String("fil", "", "specify bam and bai files containing filial genome alignments")
	pat = flag.String("pat", "", "specify bam and bai files containing paternal genome alignments")
	mat = flag.String("mat", "", "specify bam and bai files containing maternal genome alignments")
)

func main() {
	flag.Parse()
	if *pat == "" || *mat == "" || *fil == "" {
		flag.Usage()
		os.Exit(1)
	}

	p, err := newCounter(*pat)
	if err != nil {
		log.Fatal(err)
	}
	defer p.Close()
	m, err := newCounter(*mat)
	if err != nil {
		log.Fatal(err)
	}
	defer m.Close()
	f, err := newCounter(*fil)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()

	// Collate each GFF feature on stdin into
	// its group of features.
	var grps []map[string]featGroup
	sc := featio.NewScanner(gff.NewReader(os.Stdin))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		g := f.FeatAttributes.Get("Group")
		gid, err := strconv.Atoi(g)
		if err != nil {
			log.Fatalf("failed to parse group id: %v", err)
		}
		grps = add(grps, gid, f)
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during gff read: %v", err)
	}

	// For each group of features, find the counts of
	// overlapping reads.
	for gid, g := range grps {
		if g == nil {
			continue
		}
		// Iterate over each group's features, counting
		// alignmens that overlap.
		sm := sortedMap(g)
		name, n := nameHeuristic(sm)
		fmt.Printf("%d\t%d\t%s\t", gid, n, name)
		for i, t := range sm {
			if i == 0 {
				nf, err := f.overlapping(t.f)
				if err != nil {
					log.Fatal(err)
				}
				np, err := p.overlapping(t.f)
				if err != nil {
					log.Fatal(err)
				}
				nm, err := m.overlapping(t.f)
				if err != nil {
					log.Fatal(err)
				}
				fmt.Printf("%d\t%d\t%d\n", nf, np, nm)
			}
		}
	}
}

// counter is a BAM/BAI reader that counts mapped reads that overlap
// a GFF feature.
type counter struct {
	f   *os.File
	r   *bam.Reader
	h   *sam.Header
	idx *bam.Index
}

// newCounter returns a counter based on path and path.bai.
func newCounter(path string) (*counter, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("failed to open bam file: %v", err)
	}
	r, err := bam.NewReader(f, 0)
	if err != nil {
		return nil, fmt.Errorf("failed to open bam stream: %v", err)
	}

	ir, err := os.Open(path + ".bai")
	if err != nil {
		return nil, fmt.Errorf("failed to open bai file: %v", err)
	}
	idx, err := bam.ReadIndex(ir)
	if err != nil {
		return nil, fmt.Errorf("failed to open bai data: %v", err)
	}
	ir.Close()

	return &counter{f: f, r: r, h: r.Header(), idx: idx}, nil
}

// overlapping returns the number of mapped BAM reads overlapping f.
func (c *counter) overlapping(f *gff.Feature) (int, error) {
	ref, ok := getReference(c.h.Refs(), f.SeqName)
	if !ok {
		return -1, fmt.Errorf("could not find reference for %q", f.SeqName)
	}
	chunks, err := c.idx.Chunks(ref, max(0, f.FeatStart-1e4), min(ref.Len, f.FeatEnd+1e4))
	if err != nil {
		return -1, fmt.Errorf("failed to get chunks: %v", err)
	}
	it, err := bam.NewIterator(c.r, chunks)
	if err != nil {
		return -1, fmt.Errorf("failed to create iterator: %v", err)
	}
	defer it.Close()

	var n int
	for it.Next() {
		rec := it.Record()
		if rec.Start() < f.FeatStart && f.FeatEnd < rec.End() {
			n++
		}
	}
	return n, nil
}

// getReference returns the sam.Reference with the specified name.
func getReference(refs []*sam.Reference, name string) (ref *sam.Reference, ok bool) {
	for _, r := range refs {
		if r.Name() == name {
			return r, true
		}
	}
	return nil, false
}

// Close closes the bam.Reader held by the counter.
func (c *counter) Close() error {
	err := c.r.Close()
	if err != nil {
		return err
	}
	return c.f.Close()
}

type featGroup struct {
	n int // n is the GFF score of the feature.
	f *gff.Feature
}

// add the feature f to grps[f's group ID][f's repeat typ].
func add(grps []map[string]featGroup, gid int, f *gff.Feature) []map[string]featGroup {
	switch {
	case gid == len(grps):
		grps = append(grps, make(map[string]featGroup))
	case gid > len(grps) && gid < cap(grps):
		grps = grps[:gid+1]
	case gid > len(grps):
		t := make([]map[string]featGroup, gid+1)
		copy(t, grps)
		grps = t
	}
	if grps[gid] == nil {
		grps[gid] = make(map[string]featGroup)
	}
	r := f.FeatAttributes.Get("Repeat")
	typ := strings.Fields(r)[0]
	p, ok := grps[gid][typ]
	if !ok {
		grps[gid][typ] = featGroup{f: f, n: int(*f.FeatScore)}
		return grps
	}
	if f.FeatStart < p.f.FeatStart {
		p.f.FeatStart = f.FeatStart
	}
	if f.FeatEnd > p.f.FeatEnd {
		p.f.FeatEnd = f.FeatEnd
	}

	p.n += int(*f.FeatScore)
	grps[gid][typ] = p
	return grps
}

type mapElement struct {
	typ string
	n   int // n is the GFF score of the feature.
	f   *gff.Feature
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

// return a sort descending slice of the groups in g.
func sortedMap(g map[string]featGroup) []mapElement {
	m := make([]mapElement, 0, len(g))
	for typ, f := range g {
		m = append(m, mapElement{typ: typ, n: f.n, f: f.f})
	}
	sort.Sort(sort.Reverse(byCount(m)))
	return m
}

// make a reasonable guess at what the name of the repeat type
// of the features in g it.
func nameHeuristic(g []mapElement) (name string, n int) {
	if len(g) == 0 {
		return "", 0
	}

	// Majority rule.
	for _, e := range g {
		n += e.n
	}
	r := float64(g[0].n) / float64(n)
	if r > 0.5 || (r == 0.5 && len(g) > 2) {
		return g[0].typ, n
	}

	// Alu heuristic.
	if isAlu(g[0].typ) {
		return trunc(g[0].typ, 5), n
	}

	// Fusion.
	var names []string
	for _, t := range g {
		names = append(names, t.typ)
	}
	return strings.Join(names, "/"), n
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

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
