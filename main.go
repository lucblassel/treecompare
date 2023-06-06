package main

import (
	"encoding/csv"
	"errors"
	"fmt"
	"io/fs"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
	"github.com/jessevdk/go-flags"
	"github.com/schollz/progressbar/v3"
)

func init_edeges(tr *tree.Tree) ([]*tree.Edge, *tree.EdgeIndex) {
	tr.ReinitIndexes()
	edges := tr.Edges()
	index := tree.NewEdgeIndex(uint64(len(edges)*2), 0.75)
	for i, e := range edges {
		index.PutEdgeValue(e, i, e.Length())
	}

	return edges, index
}

type Stats struct {
	Identifier string
	RF         int
	normRF     float64
	wRF        float64
	KF         float64
	Error      error
}

func (self *Stats) to_csv() []string {
	if self.Error != nil {
		return []string{}
	}
	return []string{
		self.Identifier,
		strconv.Itoa(self.RF),
		strconv.FormatFloat(self.normRF, 'f', -1, 64),
		strconv.FormatFloat(self.wRF, 'f', -1, 64),
		strconv.FormatFloat(self.KF, 'f', -1, 64),
	}
}

func get_headers() []string {
	return []string{
		"id", "rf", "norm_rf", "weighted_rf", "branch_score",
	}
}

func compare(refTree *tree.Tree, compTree *tree.Tree, id string, results chan Stats) {
	var err error
	var common, ref, comp []float64

	refEdges, refIndex := init_edeges(refTree)
	compEdges, compIndex := init_edeges(compTree)

	if err = refTree.CompareTipIndexes(compTree); err != nil {
		results <- Stats{id, 0, 0.0, 0.0, 0.0, err}
		return
	}

	for _, compEdge := range compEdges {
		if !compEdge.Right().Tip() {
			refEdge, ok := refIndex.Value(compEdge)
			if !ok { // unique comp edge
				comp = append(comp, compEdge.Length())
			} else { // Edge in common
				refLen := refEdge.Len
				compLen := compEdge.Length()
				if refLen != compLen {
					common = append(common, refLen-compLen)
				}
			}
		}
	}

	for _, refEdge := range refEdges {
		if !refEdge.Right().Tip() {
			if _, ok := compIndex.Value(refEdge); !ok { // unique ref edge
				ref = append(ref, refEdge.Length())
			}
		}
	}

	rf := len(ref) + len(comp)
	n_rf := float64(rf) / float64(rf+2*len(common))
	wrf := 0.0
	kf := 0.0

	for _, diff := range common {
		wrf += math.Abs(diff)
		kf += math.Pow(diff, 2.0)
	}

	for _, container := range [][]float64{ref, comp} {
		for _, length := range container {
			wrf += math.Abs(length)
			kf += math.Pow(length, 2.0)
		}
	}

	results <- Stats{id, rf, n_rf, wrf, math.Sqrt(kf), nil}
}

func read_tree(path string) (tr *tree.Tree, err error) {
	var file *os.File

	if file, err = os.Open(path); err != nil {
		return
	}
	defer file.Close()

	tr, err = newick.NewParser(file).Parse()

	return
}

func run_job(inputs <-chan [2]string, results chan Stats) {
	var t1, t2 *tree.Tree
	var err error
	var id string

	for pair := range inputs {
		id = strings.TrimSuffix(filepath.Base(pair[0]), filepath.Ext(pair[0]))

		// Read trees
		if t1, err = read_tree(pair[0]); err != nil {
			results <- Stats{id, 0, 0., 0., 0., err}
			continue
		}
		if t2, err = read_tree(pair[1]); err != nil {
			results <- Stats{id, 0, 0., 0., 0., err}
			continue
		}

		// compare trees
		compare(t1, t2, id, results)
	}
}

func get_matching_trees(d1, d2 string) (pairs map[string][2]string, err error) {
	var walker []fs.DirEntry
	pairs = make(map[string][2]string)

	if walker, err = os.ReadDir(d1); err != nil {
		return
	}
	for _, entry := range walker {
		ext := filepath.Ext(entry.Name())
		if entry.IsDir() || !(ext == ".nwk" || ext == ".nw" || ext == ".newick") {
			continue
		}

		basename := strings.Split(entry.Name(), ".")[0]
		pairs[basename] = [2]string{filepath.Join(d1, entry.Name()), ""}
	}

	if walker, err = os.ReadDir(d2); err != nil {
		return
	}
	for _, entry := range walker {
		if entry.IsDir() {
			continue
		}

		basename := strings.Split(entry.Name(), ".")[0]
		if pair, ok := pairs[basename]; ok {
			pair[1] = filepath.Join(d2, entry.Name())
			pairs[basename] = pair
		} else {
			err = errors.New(fmt.Sprintf("Could not find matching tree for %s in %s", entry.Name(), d1))
		}
	}

	return
}

func main() {
	var opts struct {
		Output  string `short:"o" long:"output" default:"output.csv" description:"Path to output"`
		Workers int    `short:"w" long:"workers" default:"6" description:"Number of threads to use"`
		// Positional Arguments
		Dirs struct {
			Dir1 string `description:"Path to real tree directory" positional-arg-name:"REAL"`
			Dir2 string `description:"Path to predict tree directory" positional-arg-name:"PRED"`
		} `positional-args:"yes" required:"yes"`
	}

	parser := flags.NewParser(&opts, flags.Default)
	if _, err := parser.Parse(); err != nil {
		os.Exit(1)
	}

	num_workers := opts.Workers

	log.Printf(
		"Comparing %v and %v to %s (on %d threads)\n\n",
		opts.Dirs.Dir1,
		opts.Dirs.Dir2,
		opts.Output,
		opts.Workers,
	)

	// IO Stuff

	pairs, err := get_matching_trees(opts.Dirs.Dir1, opts.Dirs.Dir2)
	if err != nil {
		panic(err)
	}

	file, err := os.Create(opts.Output)
	defer file.Close()
	if err != nil {
		log.Fatal("failed to open the output file: ", err)
	}
	writer := csv.NewWriter(file)
	defer writer.Flush()

	if err := writer.Write(get_headers()); err != nil {
		log.Fatalln("Could not write header to output:", err)
	}

	// Compute distances

	inputs := make(chan [2]string, num_workers)
	results := make(chan Stats, len(pairs))

	for w := 1; w <= num_workers; w++ {
		go run_job(inputs, results)
	}

	bar := progressbar.Default(int64(len(pairs)))

	for _, pair := range pairs {
		inputs <- pair
	}

	close(inputs)

	var result Stats

	for i := 0; i < len(pairs); i++ {
		if result = <-results; result.Error != nil {
			log.Println("Error: ", result.Error)
		} else {
			if err := writer.Write(result.to_csv()); err != nil {
				log.Fatalln("Could not write CSV row:", err)
			}
		}
		bar.Add(1)
	}
}
