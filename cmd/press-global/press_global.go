// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// press-global identifies, annotates and counts unique reefer events
// without reference to a set of censor annotations.
//
// The arguments for press-global differ from press in that the input
// on stdin is the set of reefer results, rather than the censor features.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/iterator"
	"gonum.org/v1/gonum/graph/simple"
	"gonum.org/v1/gonum/graph/topo"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/store/interval"
)

var (
	thresh   = flag.Float64("thresh", 0.90, "specify minumum jaccard similarity for identity between events")
	curve    = flag.String("curve", "", "specify the tsv output file for threshold response")
	gffOut   = flag.String("gff", "", "specify the gff output file for remapping")
	deletion = flag.Bool("del", false, "specify that the input are deletions")
)

func main() {
	flag.Parse()

	var v []*gff.Feature
	trees := make(map[string]*interval.IntTree)

	sc := featio.NewScanner(gff.NewReader(os.Stdin))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		fields := strings.Fields(f.FeatAttributes.Get("Read"))
		if len(fields) != 3 {
			log.Fatalf("bad record: %+v", f)
		}
		var err error
		e := *f
		e.FeatStart, err = strconv.Atoi(fields[1])
		if err != nil {
			log.Fatalf("bad record: %+v: %v", f, err)
		}
		e.FeatEnd, err = strconv.Atoi(fields[2])
		if err != nil {
			log.Fatalf("bad record: %+v: %v", f, err)
		}
		b := baseCoordsOf(&e, f, *deletion)
		t, ok := trees[b.SeqName]
		if !ok {
			t = &interval.IntTree{}
			trees[b.SeqName] = t
		}
		t.Insert(gffInterval{id: uintptr(len(v)), Feature: b}, true)
		v = append(v, b)
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during gff read: %v", err)
	}
	for _, t := range trees {
		t.AdjustRanges()
	}

	g := thresholdGraph{WeightedUndirectedGraph: simple.NewWeightedUndirectedGraph(1, 0), thresh: *thresh}
	for i, from := range v {
		if g.Node(int64(i)) == nil {
			g.AddNode(simple.Node(i))
		}
		for _, _to := range trees[from.SeqName].Get(gffInterval{Feature: from}) {
			to := _to.(gffInterval)
			if from == to.Feature {
				continue
			}
			jac := jaccard(from, to.Feature)
			if jac > 0 {
				g.SetWeightedEdge(simple.WeightedEdge{F: simple.Node(i), T: simple.Node(to.id), W: jac})
			}
		}
	}

	cc := topo.ConnectedComponents(g)
	fmt.Printf("number of unique events = %d, total number of nodes = %d\n", len(cc), g.Nodes().Len())
	if *gffOut != "" {
		gf, err := os.Create(*gffOut)
		if err != nil {
			log.Fatalf("failed to create gff file %q: %v", *gffOut, err)
		}
		w := gff.NewWriter(gf, 60, true)
		w.WriteComment("Right coordinates (field 5) and strand (field 7) are hypothetical.")
		for i, c := range cc {
			for _, e := range c {
				f := v[e.ID()]
				f.FeatAttributes = append(f.FeatAttributes, gff.Attribute{Tag: "Group", Value: fmt.Sprint(i)})
				w.Write(f)
			}
		}
		gf.Close()
	}

	if *curve != "" {
		cf, err := os.Create(*curve)
		if err != nil {
			log.Fatalf("failed to create curve file %q: %v", *curve, err)
		}
		fmt.Fprintln(cf, "thresh\treduction")
		for g.thresh = 0.05; g.thresh < 1.04; g.thresh += 0.05 {
			fmt.Fprintf(cf, "%.2f\t%f\n", g.thresh, 1-float64(len(topo.ConnectedComponents(g)))/float64(g.Nodes().Len()))
		}
		cf.Close()
	}
}

func baseCoordsOf(f, ref *gff.Feature, isDeletion bool) *gff.Feature {
	b := *ref
	b.Source = "press/global"
	if isDeletion {
		b.Feature = "deletion"
		delta := f.Len() / 2
		b.FeatStrand = seq.None
		b.FeatStart += delta
		b.FeatEnd -= delta
		return &b
	}
	b.Feature = "insertion"
	b.FeatStrand = f.FeatStrand
	b.FeatStart = ref.FeatStart + f.FeatStart
	b.FeatEnd = ref.FeatStart + f.FeatEnd
	return &b
}

type gffInterval struct {
	id uintptr
	*gff.Feature
}

func (i gffInterval) ID() uintptr { return i.id }
func (i gffInterval) Range() interval.IntRange {
	return interval.IntRange{Start: i.FeatStart, End: i.FeatEnd}
}
func (i gffInterval) Overlap(b interval.IntRange) bool {
	return i.FeatEnd > b.Start && i.FeatStart < b.End
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

// thesholdGraph is an undirected graph where edges must be above
// a given threshold to be returned or traversed.
type thresholdGraph struct {
	*simple.WeightedUndirectedGraph
	thresh float64
}

// From returns all nodes in g that can be reached directly from n.
func (g thresholdGraph) From(n int64) graph.Nodes {
	if g.Node(n) == nil {
		return nil
	}

	var nodes []graph.Node
	for _, to := range graph.NodesOf(g.WeightedUndirectedGraph.From(n)) {
		if g.HasEdgeBetween(n, to.ID()) {
			nodes = append(nodes, to)
		}
	}

	return iterator.NewOrderedNodes(nodes)
}

// HasEdgeBetween returns whether an edge exists between nodes x and y.
func (g thresholdGraph) HasEdgeBetween(x, y int64) bool {
	if !g.WeightedUndirectedGraph.HasEdgeBetween(x, y) {
		return false
	}
	w, _ := g.Weight(x, y)
	return w >= g.thresh
}

// Edge returns the edge from u to v if such an edge exists and nil otherwise.
// The node v must be directly reachable from u as defined by the From method.
func (g thresholdGraph) Edge(u, v int64) graph.Edge {
	return g.EdgeBetween(u, v)
}

// EdgeBetween returns the edge between nodes x and y.
func (g thresholdGraph) EdgeBetween(x, y int64) graph.Edge {
	e := g.WeightedUndirectedGraph.EdgeBetween(x, y)
	if e == nil {
		return nil
	}
	if w, _ := g.Weight(x, y); w < g.thresh {
		return nil
	}
	return e
}
