#!/usr/bin/env python
# -*- coding: utf-8 -*- 
""" ShapeCompare - Command line version
A tool to visualize and compare multiple SHAPE experiment data.

Usage:
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

"""

import sys
import glob
import traceback
from docopt import docopt
from numpy import recfromtxt


import options
from processor.shapedata import ShapeConvertor, ShapeData
from processor.model import Model
from processor.fasta import Fasta


__author__ = "Guilhem Doulcier"
__copyright__ = "Copyright 2012, 2013, Guilhem Doulcier"
__license__ = "GPL"
__version__ = str(recfromtxt("../version.txt"))
__email__ = "guilhem.doulcier@ens.fr"
__date__ = "2013"

if __name__ == '__main__':

    args = docopt(__doc__, version=__version__)

    if args["--rna"] == None:
        fasta = Fasta(options.FASTA_FILE)
        args["--rna"] = fasta.sequences.keys()[-1]

    if args["--verbose"] or args["--license"]:
        print('\nShapeCompare v'+__version__+" - Command line mode")
        print("Copyright (C) 2012,2013 Guilhem DOULCIER \n    This program comes with ABSOLUTELY NO WARRANTY.\n    This program is free software: you can redistribute it and/or modify\n    it under the terms of the GNU General Public License as published by\n    the Free Software Foundation, either version 3 of the License, or\n    (at your option) any later version.")
        if args["--license"]:
            sys.exit(2)

    data = []
    sorted_data = {}
    shapeData = {}
    models = {}

    # --------------------------------------------------------------------------
    ### LOADING
    # --------------------------------------------------------------------------
    if args["--verbose"]:   
        print("\nLoading files...")
        
    for arg in args["<folders>"]:
        arg = str(arg)
        if ".shape" in arg: 
            name = arg.split("/")[-2].split(".")[0]
            try:
                data.append(ShapeConvertor(arg,name))
            except Exception, e:
                print("Error in importing: " +name+
                      " (Error code :" + str(sys.exc_info()[0])+")")
                if args["--verbose"]:
                    print traceback.format_exc()
                else:
                    print("Use -v for more information")
    

        elif arg[-1] == "/":
            for f in glob.glob(arg+"*.shape"):
                name = f.split("/")[-1].split(".")[0]
                try:
                    data.append(ShapeConvertor(f,name))
                except Exception, e:
                    print("Error in importing: " +name+
                          " (Error code :" + str(sys.exc_info()[0])+")")
                    if args["--verbose"]:
                        print traceback.format_exc()
                    else:
                        print("Use -v for more information")
                          


    if args["--verbose"]:
        print(str(len(data)) + " experiment files loaded.")


    #Loading models
    for f in glob.glob(options.MODEL_PATH+"*.csv"):
        name = f.split("/")[-1].split("_")[0]
        if name not in models.keys():
            try:
                models[name] = Model(name)
            except Exception, e:
                    print("Error in importing: " +name+" model (Error code :" + str(sys.exc_info()[0])+")")
                    if args["--verbose"]:
                          print traceback.format_exc()      
                    else:
                          print("Use -v for more information")
    

    if args["--verbose"]:
        print(str(len(models)) + " model files loaded.")


    # --------------------------------------------------------------------------
    # SORTING
    # --------------------------------------------------------------------------
    if args["--verbose"]:
        print("\nInfo:")
    for d in data:
        if args["--verbose"]:
            print("-"+ d.name+" is "+d.seq_type+", anchor: "+d.a+" ("+str(d.a_pos)+") matched with "+str(d.match)+"% of bases (offset "+str(d.offset_bases)+", try "+str(d.match_try)+")") 
        if d.seq_type in sorted_data.keys():
            sorted_data[d.seq_type].append(d)
        else:
            sorted_data[d.seq_type] = [d]
    if args["--verbose"]:
        for key, m in models.iteritems():
            print("-"+m.name+" model is "+m.rna)

    # --------------------------------------------------------------------------
    # ANALYSING
    # --------------------------------------------------------------------------
    if args["--verbose"]:
        print("\nData analysis...")
    
    for key, datum in sorted_data.iteritems():
        if args["--verbose"]:
            print("-"+key)
        try:
            shapeData[key] = ShapeData(datum,key)    
        except Exception, e:
            print("\nError in processing: "+key+" (Error:" + str(sys.exc_info()[0])+")")
            if args["--verbose"]:
                print traceback.format_exc()
            else:
                print("Use -v for more information")
    

    # --------------------------------------------------------------------------
    # DISPLAY
    # --------------------------------------------------------------------------

    if args["--interactive"]:
        print "\n Interactive mode not implemented yet.\n"
    else:
        try:
            import plots as plt
            ax = plt.fig.add_subplot(111)
            if args["histo"]:
                plt.histo.reactivity_mean(ax,shapeData[args["--rna"]])
            elif args["scatter"]:
                plt.histo.scatter(ax,shapeData[args["--rna"]])
            elif args["twod"]:
                plt.twod.mean(ax,models[args["--model"]],shapeData[args["--rna"]])
            else:
                plt.sequences.plot(ax,shapeData,models=models)
        except Exception, e:
            print("\nError in plotting: (Error:" + str(sys.exc_info()[0])+")")
            if args["--verbose"]:
                print traceback.format_exc()
            else:
                print("Use -v for more information")
    
        if args["--output-file"]:
            plt.plot.savefig(args["--output-file"])        
        else:
            plt.plot.show()

            
