// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// catch-global looks for target site duplications flanking reefer event
// output by press-global.
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
)

var in = flag.String("in", "", "specify input gff file (required)")

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
	events := make(map[string][]*gff.Feature)
	fsc := featio.NewScanner(gff.NewReader(f))
	for fsc.Next() {
		f := fsc.Feat().(*gff.Feature)
		fields := strings.Fields(f.FeatAttributes.Get("Read"))
		if len(fields) != 3 {
			log.Fatalf("bad record: %+v", f)
		}
		events[fields[0]] = append(events[fields[0]], f)
	}
	if err := fsc.Error(); err != nil {
		log.Fatalf("error during gff read: %v", err)
	}
	f.Close()

	for _, ref := range flag.Args() {
		f, err = os.Open(ref)
		if err != nil {
			log.Fatalf("failed to open reference %q: %v", ref, err)
		}
		ssc := seqio.NewScanner(fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNA)))
		for ssc.Next() {
			seq := ssc.Seq().(*linear.Seq)
			for _, f := range events[seq.Name()] {
				fields := strings.Fields(f.FeatAttributes.Get("Read"))
				if len(fields) != 3 {
					log.Fatalf("bad record: %+v", f)
				}
				start, err := strconv.Atoi(fields[1])
				if err != nil {
					log.Fatalf("failed to get start coordinate: %v", err)
				}
				end, err := strconv.Atoi(fields[2])
				if err != nil {
					log.Fatalf("failed to get end coordinate: %v", err)
				}
				tmp := *seq
				tmp.ID += fmt.Sprintf("//%d_%d", start, end)
				tmp.Seq = tmp.Seq[start:end]
				fmt.Printf("%60a\n", &tmp)
			}
		}
		if err := ssc.Error(); err != nil {
			log.Fatalf("error during fasta read: %v", err)
		}
		f.Close()
	}
}
