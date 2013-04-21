![logo](shapecompare.png "ShapeCompare logo")
ShapeCompare
============
A tool to visualize and compare multiple SHAPE experiment data.

Copyright (C) 2012, 2013  Guilhem DOULCIER

ShapeCompare is a free software for GNU/Linux, Mac OS and Windows made to visualize and compare results of high-throughput Selective 2′-Hydroxyl Acylation analyzed Primer Extension (hSHAPE) experiments. 

Features
------
See http://www.eleves.ens.fr/home/doulcier/shapecompare/ for an overview.

Usage
-----

```sh
    cli.py sequences <folders>...  [-v | -o <path>]
    cli.py (histo|scatter) <folders>...  [--rna <seq_name> | -v | -o <path>]
    cli.py twod <folders>...  --model <model_name> [--rna <seq_name> | -v | -o <path>]
    cli.py (-i|--interactive) <folders>...   [-v] 
    cli.py (-h | --help)
    cli.py --version
    cli.py --license

    Commands:
      sequence <folders>...  Align on a sequence multiple experiments.
      histo <folders>...     Bar-plot of reactivity (sequence order).
      scatter <folders>...   Scatter plot of mean reactivity (sorted).
      twod <folders>...      Mean reactivity on a secondary structure template.

    Nb. <folders>... can either be ".shape" folder(s) or folder(s) containing them.
    (In the last case all .shape sub-folders will be evaluated.)

    Options:
    -v --verbose             Execute in verbose mode.
    -i --interactive         Execute in interactive mode. (Not implemented yet).
    -o --output-file <path>  Save the plot in a file instead of displaying it.
    --rna <seq_name>         RNA sequence to show (must be in the fasta file).
    --model <model_name>     Model of secondary structure to use.
    -h --help                Show this screen.
    --version                Show version.
    --license                Show license information.
```

### Example
All experiments in the data_coreset folder:
`python cli.py histo ../data/shape/data_coreset/ -v `

One file only:
`./cli.py ../data/shape/data_coreset/7SKWT_RTH_MC5_fit.shape/ `

How To
------
How to get ShapeCompare to play nice with your data ? Follow this steps :

1. Get the the code
2. Construct the files corresponding to your favorite RNA : you will need a fasta file with the sequences, and an other with the alignment. (Inspire you from the one I've made for 7SK to get started) 
3. Modify the `options.py` file accordingly.
4. Your data should be in aligned ShapeFinder format (*.shape folder with `procTrace.aln` inside). If the sequence is different from the first of the fasta file (e.g. mutant) put a `seq.txt` file inside containing the sequence name to match.</li>
5. Enjoy ! (Don't forget the `--help` option).


Options 
-------
You can modify several aspect of shapeCompare behavior by editing the `options.py` file. 


Requirement
----------
ShapeCompare is written in python 2.7 and needs the following library :

-    [Numpy](http://www.scipy.org/) and [Scipy](http://www.scipy.org/)
-    [Matplotlib](http://matplotlib.org/)
-    [Docopt](https://github.com/docopt/docopt) (Included)
 
License
------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
