// Copyright ©2015 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// rinse removes events that are either too close to the end of a read
// or a contig, or map to a site of a repeat of the same class.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/store/interval"
)

var (
	in      = flag.String("in", "", "insertion event gff file")
	mapfile = flag.String("map", "", "read mapping gff file")
	ref     = flag.String("ref", "", "annotation gff file")
	contigs = flag.String("contigs", "", "contig fasta file")
	buf     = flag.Int("buffer", 100, "minimum distance from end of read")
)

func main() {
	flag.Parse()
	if *in == "" || *ref == "" || *mapfile == "" || *contigs == "" {
		flag.Usage()
		os.Exit(0)
	}

	refTrees, err := readAnnotations(*ref)
	if err != nil {
		log.Fatalf("failed to read annotation trees: %v", err)
	}
	mapping, err := readMappings(*mapfile)
	if err != nil {
		log.Fatalf("failed to read mapping file: %v", err)
	}
	contigLength, err := readContigs(*contigs)
	if err != nil {
		log.Fatalf("failed to read contig file: %v", err)
	}

	f, err := os.Open(*in)
	if err != nil {
		log.Fatalf("failed to open %q: %v", *in, err)
	}

	w := gff.NewWriter(os.Stdout, 60, true)

	sc := featio.NewScanner(gff.NewReader(f))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		ok, err := within(*buf, f.SeqName)
		if err != nil {
			log.Fatalf("failed to parse sequence name: %s: %v", f.SeqName, err)
		}
		if !ok {
			log.Printf("too close to read end: excluding %+v", f)
			continue
		}

		repeat := f.FeatAttributes.Get("Repeat")
		if repeat == "" {
			continue
		}
		fields := strings.Fields(repeat)

		contigSide, ok := mapping[strings.Split(f.SeqName, "//")[0]]
		if !ok {
			log.Fatalf("unexpected sequence name in input: %q", f.SeqName)
		}
		if contigSide.FeatStart < *buf || contigLength[contigSide.SeqName]-contigSide.FeatEnd < *buf {
			log.Printf("too close to contig end: excluding %+v", f)
			continue
		}
		t, ok := refTrees[contigSide.SeqName]
		if !ok {
			log.Fatalf("no tree for %v mapped by %v", contigSide.SeqName, f.SeqName)
		}
		var n int
		hits := t.Get(gffInterval{Feature: contigSide})
		for _, h := range hits {
			f := h.(gffInterval)
			repeat := f.FeatAttributes.Get("Repeat")
			if repeat == "" {
				continue
			}
			hitClass := strings.Fields(repeat)[1]
			if fields[1] == hitClass {
				n++
			}
		}
		if n != 0 {
			log.Printf("too many hits: excluding %+v", f)
			for _, h := range hits {
				log.Printf("\t%+v", h.(gffInterval).Feature)
			}
			continue
		}
		w.Write(f)
	}
	err = sc.Error()
	if err != nil {
		log.Fatalf("error during GFF read: %v", err)
	}
}

func within(buffer int, name string) (bool, error) {
	fields := strings.Split(name, "//")
	if len(fields) != 2 {
		return false, fmt.Errorf("too many fields: %q", name)
	}
	readRangeIdx := strings.LastIndex(fields[0], "/")
	if readRangeIdx < 0 {
		return false, fmt.Errorf("no path separator: %q", fields[0])
	}

	readStart, readEnd, err := underscorePair(fields[0][readRangeIdx+1:])
	if err != nil {
		return false, err
	}
	readLen := readEnd - readStart

	featStart, featEnd, err := underscorePair(strings.TrimSuffix(fields[1], "(-)"))
	if err != nil {
		return false, err
	}

	if featStart < buffer {
		return false, nil
	}
	if readLen-featEnd < buffer {
		return false, nil
	}
	return true, nil
}

func underscorePair(s string) (left, right int, err error) {
	fields := strings.Split(s, "_")
	if len(fields) != 2 {
		return 0, 0, fmt.Errorf("too many fields: %q", s)
	}
	left, err = strconv.Atoi(fields[0])
	if err != nil {
		return 0, 0, err
	}
	right, err = strconv.Atoi(fields[1])
	if err != nil {
		return 0, 0, err
	}
	return left, right, nil
}

func readMappings(file string) (map[string]*gff.Feature, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	mapping := make(map[string]*gff.Feature)
	sc := featio.NewScanner(gff.NewReader(f))
	for id := uintptr(1); sc.Next(); id++ {
		f := sc.Feat().(*gff.Feature)
		read := f.FeatAttributes.Get("Read")
		if read == "" {
			continue
		}
		// Currently reefer only expects a single hit per read.
		mapping[strings.Fields(read)[0]] = f
	}
	if err != nil {
		log.Fatalf("error during GFF read: %v", err)
	}
	return mapping, nil
}

func readContigs(file string) (map[string]int, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	lengths := make(map[string]int)
	sc := seqio.NewScanner(fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNA)))
	for sc.Next() {
		s := sc.Seq()
		lengths[s.Name()] = s.Len()
	}
	if err != nil {
		log.Fatalf("error during fasta read: %v", err)
	}
	return lengths, nil
}

func readAnnotations(file string) (map[string]*interval.IntTree, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	trees := make(map[string]*interval.IntTree)
	sc := featio.NewScanner(gff.NewReader(f))
	for id := uintptr(1); sc.Next(); id++ {
		f := sc.Feat().(*gff.Feature)
		t, ok := trees[f.SeqName]
		if !ok {
			t = &interval.IntTree{}
			trees[f.SeqName] = t
		}
		t.Insert(gffInterval{f, id}, true)
	}
	err = sc.Error()
	if err != nil {
		log.Fatalf("error during GFF read: %v", err)
	}
	for _, t := range trees {
		t.AdjustRanges()
	}
	return trees, nil
}

type gffInterval struct {
	*gff.Feature
	id uintptr
}

func (f gffInterval) ID() uintptr { return f.id }
func (f gffInterval) Range() interval.IntRange {
	return interval.IntRange{Start: f.Feature.FeatStart, End: f.Feature.FeatEnd}
}
func (f gffInterval) Overlap(b interval.IntRange) bool {
	// Half-open interval indexing.
	return f.Feature.FeatEnd > b.Start && f.Feature.FeatStart < b.End
}
