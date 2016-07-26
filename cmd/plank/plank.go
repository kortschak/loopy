// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// plank drops GFF lines from stdin containing Read attributes in
// the exclude parameter file.
package main

import (
	"bufio"
	"flag"
	"log"
	"os"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
)

var (
	exclude = flag.String("exclude", "", "specify file containing excluded reads")
	retain  = flag.Bool("retain", false, "write excluded reads to stderr")
)

func main() {
	flag.Parse()
	if *exclude == "" {
		flag.Usage()
		os.Exit(1)
	}

	nameSet := make(map[string]struct{})
	f, err := os.Open(*exclude)
	if err != nil {
		log.Fatalf("failed to open exclude file %q: %v", *exclude, err)
	}
	ls := bufio.NewScanner(f)
	for ls.Scan() {
		nameSet[ls.Text()] = struct{}{}
	}
	err = ls.Err()
	if err != nil {
		log.Fatalf("failed to read exclude file: %v", err)
	}

	w := gff.NewWriter(os.Stdout, 60, true)
	var excl *gff.Writer
	if *retain {
		excl = gff.NewWriter(os.Stderr, 60, true)
	}
	sc := featio.NewScanner(gff.NewReader(os.Stdin))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		n := f.FeatAttributes.Get("Read")
		if _, ok := nameSet[n]; ok {
			if excl != nil {
				_, err := excl.Write(f)
				if err != nil {
					log.Fatalf("failed to write feature: %v", err)
				}
			}
			continue
		}
		_, err := w.Write(f)
		if err != nil {
			log.Fatalf("failed to write feature: %v", err)
		}
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during gff read: %v", err)
	}
}
