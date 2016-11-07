// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// fathom filters events based on length of element, reading from stdin.
package main

import (
	"flag"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
)

var thresh = flag.Int("thresh", 0, "specify minimum element length")

func main() {
	flag.Parse()

	w := gff.NewWriter(os.Stdout, 60, false)
	sc := featio.NewScanner(gff.NewReader(os.Stdin))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		r := f.FeatAttributes.Get("Repeat")
		fields := strings.Fields(r)
		if len(fields) < 4 {
			log.Fatal("invalid repeat attribute")
		}
		end, err := strconv.Atoi(fields[3])
		if err != nil {
			log.Fatalf("failed to parse end coordinate: %v", err)
		}
		remainder, err := strconv.Atoi(fields[4])
		if err != nil {
			log.Fatalf("failed to parse remains coordinate: %v", err)
		}
		length := end + remainder
		if length < *thresh {
			continue
		}
		w.Write(f)
	}
}
