// Copyright ©2015 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Package blasr provides interaction with the BLASR long read aligner.
package blasr

import (
	"errors"
	"fmt"
	"os/exec"
	"sort"
	"strings"
	"text/template"

	"github.com/biogo/external"
)

var ErrMissingRequired = errors.New("blasr: missing required argument")

// BLASR defines parameters for the blasr aligner.
type BLASR struct {
	// Usage: blasr reads.{bam|fasta|bax.h5|fofn} genome.fasta [-options]
	//
	Cmd string `buildarg:"{{if .}}{{.}}{{else}}blasr{{end}}"` // blasr

	// Input Files:
	Reads  string `buildarg:"{{.}}"` // "reads.{bam|fasta|bax.h5|fofn}"
	Genome string `buildarg:"{{.}}"` // "genome.fasta"

	SuffixArray string `buildarg:"{{if .}}--sa{{split}}{{.}}{{end}}"`   // -sa: suffix array file
	TupleCounts string `buildarg:"{{if .}}--ctab{{split}}{{.}}{{end}}"` // -ctab: table of tuple counts

	// Output file options:
	Aligned string `buildarg:"{{if .}}--out{{split}}{{.}}{{end}}"` // -out: outfile (stdout if empty)

	// SAM output options:
	SAM           bool   `buildarg:"{{if .}}--sam{{end}}"`                    // -sam: write output in SAM format
	Clipping      string `buildarg:"{{if .}}--clipping{{split}}{{.}}{{end}}"` // -clipping: no/hard/subread/soft clipping for SAM
	SAMQV         bool   `buildarg:"{{if .}}--printSAMQV{{end}}"`             // -printSAMQV: quality values to SAM output
	CIGARSeqMatch bool   `buildarg:"{{if .}}--cigarUseSeqMatch{{end}}"`       // -cigarUseSeqMatch: use '=' and 'X' to represent match

	// Non-SAM output options:
	Format int  `buildarg:"{{if .}}--m{{split}}{{.}}{{end}}"` // -m: output format if !SAM
	Header bool `buildarg:"{{if .}}--header{{end}}"`          // -header: output header

	TitleTable string `buildarg:"{{if .}}--unaligned{{split}}{{.}}{{end}}"` // -titleTable: output table of reference sequence titles

	Unaligned string `buildarg:"{{if .}}--unaligned{{split}}{{.}}{{end}}"` // -unaligned: outfile for unaligned reads

	// Alignment options:
	MinSeedLength       int    `buildarg:"{{if .}}--minMatch{{split}}{{.}}{{end}}"`              // -minMatch: minimum seed length
	MaxAligmentLength   int    `buildarg:"{{if .}}--maxMatch{{split}}{{.}}{{end}}"`              // -maxMatch: maximum longest common prefix
	MaxAnchors          int    `buildarg:"{{if .}}--maxAnchorsPerPosition{{split}}{{.}}{{end}}"` // -maxAnchorsPerPosition
	AdvanceExactMatches int    `buildarg:"{{if .}}--advanceExactMatches{{split}}{{.}}{{end}}"`   // -advanceExactMatches
	Candidates          int    `buildarg:"{{if .}}--nCandidates{{split}}{{.}}{{end}}"`           // -nCandidates: candidates to keep for the best alignment
	Concordant          bool   `buildarg:"{{if .}}--concordant{{end}}"`                          // -concordant
	ConcordantTemplate  string `buildarg:"{{if .}}--concordantTemplate{{split}}{{.}}{{end}}"`    // -concordantTemplate
	FastMaxInterval     bool   `buildarg:"{{if .}}--fastMaxInterval{{end}}"`                     // -fastMaxInterval
	HardIntervalCut     bool   `buildarg:"{{if .}}--aggressiveIntervalCut{{end}}"`               // -aggressiveIntervalCut
	FastSDP             bool   `buildarg:"{{if .}}--fastSDP{{end}}"`                             // -fastSDP

	// Dynamic programming and overlap options:
	UseQuality  bool `buildarg:"{{if .}}--useQuality{{end}}"`  // -useQuality
	AffineAlign bool `buildarg:"{{if .}}--affineAlign{{end}}"` // -affineAlign

	// Hit refinement options:
	SDPTupleSize int    `buildarg:"{{if .}}--sdpTupleSize{{split}}{{.}}{{end}}"` // -sdpTupleSize
	ScoreMatrix  string `buildarg:"{{if .}}--scoreMatrix{{split}}{{.}}{{end}}"`  // -scoreMatrix
	AffineOpen   int    `buildarg:"{{if .}}--affineOpen{{split}}{{.}}{{end}}"`   // -affineOpen: gap open penalty
	AffineExtend int    `buildarg:"{{if .}}--affineExtend{{split}}{{.}}{{end}}"` // -affineOpen: gap extend penalty

	// Read and alignment filtering and reporting options:
	BestN              int     `buildarg:"{{if .}}--bestn{{split}}{{.}}{{end}}"`              // -bestn: report the top 'n' alignments
	HitPolicy          string  `buildarg:"{{if .}}--hitPolicy{{split}}{{.}}{{end}}"`          // -hitPolicy: policy to treat multiple hits
	RandomSeed         int     `buildarg:"{{if .}}--randomSeed{{split}}{{.}}{{end}}"`         // -randomSeed: prng seed
	NoResort           bool    `buildarg:"{{if .}}--noSortRefinedAlignments{{end}}"`          // -noSortRefinedAlignments
	AdjacentIndels     bool    `buildarg:"{{if .}}--allowAdjacentIndels{{end}}"`              // -allowAdjacentIndels
	MinReadLength      int     `buildarg:"{{if .}}--minReadLength{{split}}{{.}}{{end}}"`      // -minReadLength
	MinSubreadLength   int     `buildarg:"{{if .}}--minSubreadLength{{split}}{{.}}{{end}}"`   // -minSubreadLength
	MinRawSubreadScore int     `buildarg:"{{if .}}--minRawSubreadScore{{split}}{{.}}{{end}}"` // -minRawSubreadScore
	MaxScore           int     `buildarg:"{{if .}}--maxScore{{split}}{{.}}{{end}}"`           // -maxScore
	MinAlignmentLength int     `buildarg:"{{if .}}--minAlnLength{{split}}{{.}}{{end}}"`       // -minAlnLength
	MinSimilarity      float64 `buildarg:"{{if .}}--minPctSimilarity{{split}}{{.}}{{end}}"`   // -minPctSimilarity: candidates to keep for the best alignment
	MinAccuracy        float64 `buildarg:"{{if .}}--minPctAccuracy{{split}}{{.}}{{end}}"`     // -minPctAccuracy: candidates to keep for the best alignment

	// Subsampling options:
	SubSampleFraction float64 `buildarg:"{{if .}}--subsample{{split}}{{.}}{{end}}"`         // -subsample: fraction of reads to randomly subsample
	HoleNumbers       []int   `buildarg:"{{if .}}--holeNumbers{{split}}{{holes .}}{{end}}"` // -holeNumbers

	// Parallel alignment options:
	Procs  int `buildarg:"{{if .}}--nproc{{split}}{{.}}{{end}}"`  // -nproc: number of processes
	Start  int `buildarg:"{{if .}}--start{{split}}{{.}}{{end}}"`  // -start: index of the first read to begin aligning
	Stride int `buildarg:"{{if .}}--stride{{split}}{{.}}{{end}}"` // -stride: stride over reads
}

// BuildCommand returns an exec.Cmd built from the parameters in b.
func (b BLASR) BuildCommand() (*exec.Cmd, error) {
	if b.Reads == "" || b.Genome == "" {
		return nil, ErrMissingRequired
	}
	cl := external.Must(external.Build(b, template.FuncMap{"holes": holes}))
	return exec.Command(cl[0], cl[1:]...), nil
}

// holes returns a string representation of a list of integers where
// sequential runs are condensed.
func holes(a interface{}) string {
	holes := a.([]int)
	sort.Ints(holes)

	// Make sure the list only contains unique values.
	j := 0
	for i := 1; i < len(holes); i++ {
		if holes[j] >= holes[i] {
			continue
		}
		j++
		if j < i {
			holes[i], holes[j] = holes[j], holes[i]
		}
	}
	holes = holes[:j+1]

	// Format the list into runs where possible.
	var s []string
	for i := 0; i < len(holes); {
		j := i
		for ; j < len(holes) && holes[j]-holes[i] <= j-i; j++ {
		}
		if i == j-1 {
			s = append(s, fmt.Sprint(holes[i]))
		} else {
			s = append(s, fmt.Sprintf("%d-%d", holes[i], holes[j-1]))
		}
		i = j
	}

	return strings.Join(s, ",")
}
