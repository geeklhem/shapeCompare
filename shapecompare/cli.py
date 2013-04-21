#!/usr/bin/env python
# -*- coding: utf-8 -*- 
""" Command line version of the shapeCompare tool"""

import sys
import glob
import traceback
from numpy import recfromtxt

import options
from processor.shapedata import ShapeConvertor, ShapeData
from processor.model import Model

__author__ = "Guilhem Doulcier"
__copyright__ = "Copyright 2012, 2013, Guilhem Doulcier"
__license__ = "GPL"
__version__ = str(recfromtxt("../version.txt"))
__email__ = "guilhem.doulcier@ens.fr"
__date__ = "2013"

if __name__ == '__main__':
    
    args = sys.argv[1:]
    verbose = False
    if '-v' in args:
        verbose = True
    if len(args) == 0:
        print("\nUsage : cli.py folders (*.shape or a folder containing them) [OPTIONS]...\n")
        print("OPTIONS:")
        print("-v : Verbose mode")
        print("-l : licence informations")
        sys.exit(2)
   
    if verbose:
        print('\nShapeCompare v'+__version__+" - Command line mode")
        print("Copyright (C) 2012,2013 Guilhem DOULCIER \n    This program comes with ABSOLUTELY NO WARRANTY.\n    This program is free software: you can redistribute it and/or modify\n    it under the terms of the GNU General Public License as published by\n    the Free Software Foundation, either version 3 of the License, or\n    (at your option) any later version.")
   
    data = []
    sorted_data = {}
    shapeData = {}
    models = {}

    # --------------------------------------------------------------------------
    ### LOADING
    # --------------------------------------------------------------------------
    if verbose:   
        print("\nLoading files...")
        
    for arg in args:
        arg = str(arg)
        if ".shape" in arg: 
            name = arg.split("/")[-2].split(".")[0]
            try:
                data.append(ShapeConvertor(arg,name))
            except Exception, e:
                print("Error in importing: " +name+" (Error code :" + str(sys.exc_info()[0])+")")
                print traceback.format_exc()

        elif arg[-1] == "/":
            for f in glob.glob(arg+"*.shape"):
                name = f.split("/")[-1].split(".")[0]
                try:
                    data.append(ShapeConvertor(f,name))
                except Exception, e:
                    print("Error in importing: " +name+" (Error code :" + str(sys.exc_info()[0])+")")
                    print traceback.format_exc()



    if verbose:
        print(str(len(data)) + " experiment files loaded.")


    #Loading models
    for f in glob.glob(options.MODEL_PATH+"*.csv"):
        name = f.split("/")[-1].split("_")[0]
        if name not in models.keys():
            try:
                models[name] = Model(name)
            except Exception, e:
                    print("Error in importing: " +name+" model (Error code :" + str(sys.exc_info()[0])+")")

    if verbose:
        print(str(len(models)) + " model files loaded.")


    # --------------------------------------------------------------------------
    # SORTING
    # --------------------------------------------------------------------------
    print("\nInfo:")
    for d in data:
        print("-"+ d.name+" is "+d.seq_type+", anchor: "+d.a+" ("+str(d.a_pos)+") matched with "+str(d.match)+"% of bases (offset "+str(d.offset_bases)+", try "+str(d.match_try)+")") 
        if d.seq_type in sorted_data.keys():
            sorted_data[d.seq_type].append(d)
        else:
            sorted_data[d.seq_type] = [d]
    for key, m in models.iteritems():
        print("-"+m.name+" model is "+m.rna)

    # --------------------------------------------------------------------------
    # ANALYSING
    # --------------------------------------------------------------------------
    if verbose:
        print("\nData analysis...")
    
    for key, datum in sorted_data.iteritems():
        if verbose:
            print("-"+key)
        try:
            pass
            shapeData[key] = ShapeData(datum,key)    
        except Exception, e:
            print("\nError in processing: "+key+" (Error:" + str(sys.exc_info()[0])+")")
            print traceback.format_exc()


    # --------------------------------------------------------------------------
    # DISPLAY
    # --------------------------------------------------------------------------
    import plots as plt
    ax = plt.fig.add_subplot(111)
    plt.sequences.plot(ax,shapeData,models=models)
    plt.plot.show()
    
