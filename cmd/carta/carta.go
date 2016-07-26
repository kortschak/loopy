// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// carta renders a rings plot of a binned feature distribution on hg19.
package main

import (
	"flag"
	"fmt"
	"image/color"
	"math"
	"os"
	"path/filepath"
	"strings"

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/feat/genome"
	"github.com/biogo/biogo/feat/genome/human/hg19"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/bed"
	"github.com/biogo/graphics/rings"

	"github.com/gonum/plot"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/vg"
	"github.com/gonum/plot/vg/draw"
)

var (
	in     string
	format string

	binLength int
)

const (
	all = iota
	primary
	secondary
)

func init() {
	flag.StringVar(&in, "in", "", "file name of a BED file to be processed.")
	flag.IntVar(&binLength, "length", 1e6, "specifies the density bin length.")
	flag.StringVar(&format, "format", "svg", "specifies the output format of the example: eps, jpg, jpeg, pdf, png, svg, and tiff.")
	help := flag.Bool("help", false, "output this usage message.")
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	if in == "" {
		flag.Usage()
		os.Exit(1)
	}
	for _, s := range []string{"eps", "jpg", "jpeg", "pdf", "png", "svg", "tiff"} {
		if format == s {
			return
		}
	}
	flag.Usage()
	os.Exit(1)
}

var index = map[string]int{}

func init() {
	for i, c := range hg19.Chromosomes {
		index[strings.ToLower(c.Chr)] = i
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func main() {
	bf, err := readBED(in)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	p, err := plot.New()
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	hs, err := tracks(scoreFeatures(bf, binLength, hg19.Chromosomes), 15*vg.Centimeter)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	p.Add(hs...)

	p.HideAxes()

	font, err := vg.MakeFont("Helvetica", 14)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	p.Title.Text = filepath.Base(in)
	p.Title.TextStyle = draw.TextStyle{Color: color.Gray{0}, Font: font}

	err = p.Save(19*vg.Centimeter, 25*vg.Centimeter, filepath.Base(in)+"."+format)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func readBED(in string) ([]*bed.Bed3, error) {
	bf, err := os.Open(in)
	if err != nil {
		return nil, err
	}
	defer bf.Close()

	br, err := bed.NewReader(bf, 3)
	if err != nil {
		return nil, err
	}

	var fs []*bed.Bed3
	sc := featio.NewScanner(br)
	for sc.Next() {
		fs = append(fs, sc.Feat().(*bed.Bed3))
	}
	err = sc.Error()
	if err != nil {
		return nil, err
	}
	return fs, nil
}

func scoreFeatures(b []*bed.Bed3, length int, gen []*genome.Chromosome) []rings.Scorer {
	var n int
	gs := make([][]*feature, len(gen))
	for i, c := range gen {
		bins := make([]*feature, (c.Len()-1)/length+1)
		n += len(bins)
		for j := range bins {
			bins[j] = &feature{
				start: j * length,
				end:   min(c.Len(), (j+1)*length),
				chr:   c,
			}
		}
		gs[i] = bins
	}
	for _, f := range b {
		gs[index[strings.ToLower(f.Chrom)]][(f.Start()+f.End())/2/length].events++
	}

	s := make([]rings.Scorer, 0, n)
	for _, c := range gs {
		for _, b := range c {
			s = append(s, b)
		}
	}
	return s
}

type feature struct {
	start, end int
	name       string
	chr        feat.Feature
	events     int
}

func (f *feature) Start() int             { return f.start }
func (f *feature) End() int               { return f.end }
func (f *feature) Len() int               { return f.end - f.start }
func (f *feature) Name() string           { return f.name }
func (f *feature) Description() string    { return "alignment bin" }
func (f *feature) Location() feat.Feature { return f.chr }
func (f *feature) Scores() []float64 {
	factor := float64(binLength) / float64(f.Len())
	return []float64{float64(f.events) * factor}
}

func tracks(scores []rings.Scorer, diameter vg.Length) ([]plot.Plotter, error) {
	var p []plot.Plotter

	radius := diameter / 2

	// Relative sizes.
	const (
		gap = 0.005

		label = 117. / 110.

		countsInner = 97. / 110.
		countsOuter = 70. / 110.

		karyotypeInner = 100. / 110.
		karyotypeOuter = 1.

		large = 6. / 110.
		small = 2. / 110.
	)

	sty := plotter.DefaultLineStyle
	sty.Width /= 2

	chr := make([]feat.Feature, len(hg19.Chromosomes))
	for i, c := range hg19.Chromosomes {
		chr[i] = c
	}
	hs, err := rings.NewGappedBlocks(
		chr,
		rings.Arc{rings.Complete / 4 * rings.CounterClockwise, rings.Complete * rings.Clockwise},
		radius*karyotypeInner, radius*karyotypeOuter, gap,
	)
	if err != nil {
		return nil, err
	}
	hs.LineStyle = sty

	p = append(p, hs)

	bands := make([]feat.Feature, len(hg19.Bands))
	cens := make([]feat.Feature, 0, len(hg19.Chromosomes))
	for i, b := range hg19.Bands {
		bands[i] = colorBand{b}
		s := b.Start()
		// This condition depends on p -> q sort order in the $karyotype.Bands variable.
		// All standard genome packages follow this, though here the test is more general than
		// actually required since hs is telocentric.
		if b.Band[0] == 'q' && (s == 0 || hg19.Bands[i-1].Band[0] == 'p') {
			cens = append(cens, colorBand{&genome.Band{Band: "cen", Desc: "Band", StartPos: s, EndPos: s, Giemsa: "acen", Chr: b.Location()}})
		}
	}
	b, err := rings.NewBlocks(bands, hs, radius*karyotypeInner, radius*karyotypeOuter)
	if err != nil {
		return nil, fmt.Errorf("bands: %v", err)
	}
	p = append(p, b)
	c, err := rings.NewBlocks(cens, hs, radius*karyotypeInner, radius*karyotypeOuter)
	if err != nil {
		return nil, fmt.Errorf("centromeres: %v", err)
	}
	p = append(p, c)

	font, err := vg.MakeFont("Helvetica", radius*large)
	if err != nil {
		return nil, err
	}
	lb, err := rings.NewLabels(hs, radius*label, rings.NameLabels(hs.Set)...)
	if err != nil {
		return nil, err
	}
	lb.TextStyle = draw.TextStyle{Color: color.Gray16{0}, Font: font}
	p = append(p, lb)

	smallFont, err := vg.MakeFont("Helvetica", radius*small)
	if err != nil {
		return nil, err
	}

	counts := make([]rings.Scorer, len(scores))
	for i, s := range scores {
		counts[i] = s.(*feature)
	}
	ct, err := rings.NewScores(counts, hs, radius*countsInner, radius*countsOuter,
		&rings.Trace{
			LineStyles: func() []draw.LineStyle {
				ls := []draw.LineStyle{sty}
				ls[0].Color = color.Gray16{0}
				return ls
			}(),
			Join: true,
			Axis: &rings.Axis{
				Angle:     rings.Complete / 4,
				Grid:      plotter.DefaultGridLineStyle,
				LineStyle: sty,
				Tick: rings.TickConfig{
					Marker:    plot.DefaultTicks{},
					LineStyle: sty,
					Length:    2,
					Label:     draw.TextStyle{Color: color.Gray16{0}, Font: smallFont},
				},
			},
		},
	)
	if err != nil {
		return nil, err
	}
	p = append(p, ct)

	return p, nil
}

type colorBand struct {
	*genome.Band
}

func (b colorBand) FillColor() color.Color {
	switch b.Giemsa {
	case "acen":
		return color.RGBA{R: 0xff, A: 0xff}
	case "gvar":
		return color.RGBA{R: 0xff, G: 0x8c, A: 0xff}
	case "stalk":
		return color.RGBA{G: 0x8c, A: 0x80}
	case "gneg":
		return color.Gray{0xff}
	case "gpos25":
		return color.Gray{3 * math.MaxUint8 / 4}
	case "gpos33":
		return color.Gray{2 * math.MaxUint8 / 3}
	case "gpos50":
		return color.Gray{math.MaxUint8 / 2}
	case "gpos66":
		return color.Gray{math.MaxUint8 / 3}
	case "gpos75":
		return color.Gray{math.MaxUint8 / 4}
	case "gpos100":
		return color.Gray{0x0}
	default:
		panic("unexpected giemsa value: " + b.Giemsa)
	}
}

func (b colorBand) LineStyle() draw.LineStyle {
	switch b.Giemsa {
	case "acen":
		return draw.LineStyle{Color: color.RGBA{R: 0xff, A: 0xff}, Width: 1}
	case "gneg", "gvar", "stalk", "gpos25", "gpos33", "gpos50", "gpos66", "gpos75", "gpos100":
		return draw.LineStyle{}
	default:
		panic("unexpected giemsa value: " + b.Giemsa)
	}
}
