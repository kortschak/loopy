// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// catch looks for target site duplications flanking reefer event
// output by press.
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/biogo/align"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

type mat [3]int

var alnmat = mat{1, -2, -3}

func (v *mat) Set(s string) error {
	fields := strings.Split(s, ",")
	if len(fields) != 3 {
		return fmt.Errorf("invalid number of fields: %q", s)
	}
	var err error
	for i, f := range fields {
		v[i], err = strconv.Atoi(f)
		if err != nil {
			return fmt.Errorf("invalid fields: %v", err)
		}
	}
	return nil
}

func (v *mat) String() string { return fmt.Sprintf("%d,%d,%d", v[0], v[1], v[2]) }

var (
	in       = flag.String("in", "", "input gff file (required)")
	thresh   = flag.Int("thresh", 6, "minimum TSD half alignment length (ungapped)")
	window   = flag.Int("window", 200, "window for TSD search")
	fastaOut = flag.String("fasta-out", "", "write insertions to this file if option not empty")
)

func main() {
	flag.Var(&alnmat, "align", "specify the match, mismatch and gap parameters")
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

	w := gff.NewWriter(os.Stdout, 60, true)
	w.WriteComment("Right coordinates (field 5) and strand (field 7) are hypothetical.")

	var out *os.File
	if *fastaOut != "" {
		out, err = os.Create(*fastaOut)
		if err != nil {
			log.Fatalf("failed to create fasta insertion output file %q: %v", *fastaOut, err)
		}
		defer out.Close()
	}

	hw := *window / 2
	sw := makeTable(alphabet.DNAgapped, alnmat)
	for _, ref := range flag.Args() {
		f, err = os.Open(ref)
		if err != nil {
			log.Fatalf("failed to open reference %q: %v", ref, err)
		}
		ssc := seqio.NewScanner(fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNAgapped)))
	loop:
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

				if out != nil {
					insert := *seq
					if insert.Desc != "" {
						insert.Desc += " "
					}
					insert.Desc += fmt.Sprintf("[%d,%d)", start, end)
					insert.Seq = insert.Seq[start:end]
					fmt.Fprintf(out, "%60a\n", &insert)
				}

				lOff := max(0, start-hw)
				lEnd := min(len(seq.Seq), start+hw)
				rOff := max(0, end-hw)
				rEnd := min(len(seq.Seq), end+hw)

				// Ensure windows don't overlap.
				if lEnd > rOff {
					lEnd = (lEnd + rOff) / 2
					rOff = lEnd
				}

				left := *seq
				left.ID = "prefix"
				left.Seq = left.Seq[lOff:lEnd]
				right := *seq
				right.ID = "postfix"
				right.Seq = right.Seq[rOff:rEnd]

				aln, err := sw.Align(&right, &left)
				if err != nil {
					log.Fatal(err)
				}

				fa := align.Format(&right, &left, aln, '-')
				for _, seg := range fa {
					var n int
					for _, l := range seg.(alphabet.Letters) {
						if l != '-' {
							n++
						}
					}
					if n < *thresh {
						continue loop
					}
				}

				var sc int
				for _, seg := range aln {
					type scorer interface {
						Score() int
					}
					sc += seg.(scorer).Score()
				}
				f.FeatAttributes = append(f.FeatAttributes, gff.Attribute{
					Tag: "TSD", Value: fmt.Sprintf(`%v %d %d %v "%v" %d`,
						fa[0], aln[len(aln)-1].Features()[0].End()+lOff,
						aln[0].Features()[1].Start()+rOff, fa[1],
						aln, sc),
				})
				w.Write(f)
			}
		}
		if err := ssc.Error(); err != nil {
			log.Fatalf("error during fasta read: %v", err)
		}
		f.Close()
	}
}

func makeTable(alpha alphabet.Alphabet, alnmat mat) align.SW {
	match := alnmat[0]
	mismatch := alnmat[1]
	gap := alnmat[2]
	sw := make(align.SW, alpha.Len())
	for i := range sw {
		row := make([]int, alpha.Len())
		for j := range row {
			row[j] = mismatch
		}
		row[i] = match
		sw[i] = row
	}
	for i := range sw {
		sw[0][i] = gap
		sw[i][0] = gap
	}
	return sw
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
