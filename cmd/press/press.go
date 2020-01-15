// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// press identifies, annotates and counts unique reefer events.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strings"

	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/iterator"
	"gonum.org/v1/gonum/graph/simple"
	"gonum.org/v1/gonum/graph/topo"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
)

var (
	in     = flag.String("in", "", "specify input gff file (required)")
	ref    = flag.String("ref", "", "specify input reference gff file (required)")
	thresh = flag.Float64("thresh", 0.90, "specify minumum jaccard similarity for identity between events")
	curve  = flag.String("curve", "", "specify the tsv output file for threshold response")
	gffOut = flag.String("gff", "", "specify the gff output file for remapping")
)

func main() {
	flag.Parse()
	if *in == "" || *ref == "" {
		flag.Usage()
		os.Exit(1)
	}

	f, err := os.Open(*in)
	if err != nil {
		log.Fatalf("failed to open %q: %v", *in, err)
	}
	events := make(map[string]*gff.Feature)
	got := make(map[string]bool)
	sc := featio.NewScanner(gff.NewReader(f))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		events[strings.TrimSuffix(f.SeqName, "(-)")] = f
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during gff read: %v", err)
	}
	f.Close()

	f, err = os.Open(*ref)
	if err != nil {
		log.Fatalf("failed to open %q: %v", *ref, err)
	}
	var v []*gff.Feature
	sc = featio.NewScanner(gff.NewReader(f))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		fields := strings.Fields(f.FeatAttributes.Get("Read"))
		if len(fields) != 3 {
			log.Fatalf("bad record: %+v", f)
		}
		e, ok := events[fmt.Sprintf("%s//%s_%s", fields[0], fields[1], fields[2])]
		if ok {
			got[fmt.Sprintf("%s//%s_%s", fields[0], fields[1], fields[2])] = true
			v = append(v, baseCoordsOf(e, f))
		}
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during gff read: %v", err)
	}
	f.Close()

	if len(events) != len(v) {
		log.Println("failed to collect all reference features:")
		for k := range events {
			if !got[k] {
				log.Printf("missing: %s", k)
			}
		}
		log.Fatal("terminating")
	}

	g := thresholdGraph{WeightedUndirectedGraph: simple.NewWeightedUndirectedGraph(1, 0), thresh: *thresh}
	// The sets of event are small at this stage,
	// so we do things the naive way rather than
	// setting up a set of interval trees.
	for i := range v[:len(v)-1] {
		for j := range v[i+1:] {
			g.SetWeightedEdge(simple.WeightedEdge{F: simple.Node(i), T: simple.Node(j + i + 1), W: jaccard(v[i], v[j+i+1])})
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

func baseCoordsOf(f, ref *gff.Feature) *gff.Feature {
	b := *ref
	b.Source = "press"
	b.Feature = "insertion"
	b.FeatStrand = f.FeatStrand
	b.FeatAttributes = append(f.FeatAttributes, ref.FeatAttributes...)
	b.FeatStart = ref.FeatStart + f.FeatStart
	b.FeatEnd = ref.FeatStart + f.FeatEnd
	return &b
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
