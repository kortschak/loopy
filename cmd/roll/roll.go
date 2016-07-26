// Copyright ©2016 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// roll outputs a list of read names from a reefer (or later) GFF
// outout on stdin. It calls the roll.
package main

import (
	"fmt"
	"log"
	"os"
	"sort"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
)

func main() {
	nameSet := make(map[string]struct{})
	sc := featio.NewScanner(gff.NewReader(os.Stdin))
	for sc.Next() {
		f := sc.Feat().(*gff.Feature)
		n := f.FeatAttributes.Get("Read")
		if n == "" {
			continue
		}
		nameSet[n] = struct{}{}
	}
	if err := sc.Error(); err != nil {
		log.Fatalf("error during gff read: %v", err)
	}

	names := make([]string, 0, len(nameSet))
	for n := range nameSet {
		names = append(names, n)
	}
	sort.Strings(names)
	for _, n := range names {
		fmt.Println(n)
	}
}
