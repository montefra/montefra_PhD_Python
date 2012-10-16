#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Out of the input catalogue, write in two separated files the particles of samples
d1 and d2 (as defined in Guo, Zehavi & Zheng 2011).
The criteria for separation are saved in a file"""

import argparse as ap
import my_functions as mf
import numpy as np
import sys  #system stuff

def parse(argv):
  """
  This function accept a list of strings, create and fill a parser istance 
  and return a populated namespace

  Parameters
  ----------
  argv: list of strings
    list to be parsed

  output: namespace
  ---------
  """
  
  description = """Divide the input catalogues in two as in Guo, Zehavi & Zheng 2011"""
  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("catname", action="store", type=ap.FileType('r'), 
      help="File name of the catalogue to split. ra, dec and z must be in the first three columns")
  p.add_argument("sepname", action="store", nargs="?", type=ap.FileType('r'), 
      help="""Name of the file containing the criteria for the separation in two columns.
      The first column must contain 1 or 2, and is used to separate the catalogue;
      the second column must contain 1 or 0, which is used to keep the particle or not.
      If not used, this information are assumed to be in two columns of 'catname'""")

  ir_group = p.add_mutually_exclusive_group()
  ir_group.add_argument("-i", "--insert", action="store", nargs=3, default=["d1", "d2", ".dat"], 
      help="""Output file names created inserting '%(dest)s[0]' or '%(dest)s[1]'
      before '%(dest)s[2]'. The output file has the same structure of the input file if 'sepname'
      given, otherwise only the first three columns are saved.""")
  ir_group.add_argument("-r", "--replace", action="store", nargs=3, 
      help="""Output file names created replacing '%(dest)s[0]' with '%(dest)s[1]' or
      '%(dest)s[2]'. The output file has the same structure of the input file if 'sepname'
      given, otherwise only the first three columns are saved.""")

  p.add_argument('--version', action='version', version='0.1')
  p.add_argument("-v", "--verbose", action="store_true")

  p.add_argument("-c", "--columns", action="store", nargs=2, default=(-3,-2), 
      help="Columns containing the the information to separate the catalogue.")

  p.add_argument("--overwrite", action="store_true", 
      help="If given does not check if the output file exists or not before saving it.")

  p.add_argument("--fmt", default="%7.6e", action="store", nargs='+', help="Format of the output files")

  return p.parse_args(args=argv)

if __name__ == "__main__":   #if it's the main

  args = parse(sys.argv[1:])

  #create the output file names
  if(args.replace == None):   #first file
    ofiled1, skip = mf.insert_src(args.catname.name, args.insert[::2], overwrite=args.overwrite)
    ofiled2, skip = mf.insert_src(args.catname.name, args.insert[1:], overwrite=args.overwrite)
  else:
    ofiled1, skip = mf.replace_src(args.catname.name, args.replace[:2], overwrite=args.overwrite)
    ofiled2, skip = mf.replace_src(args.catname.name, args.replace[::2], overwrite=args.overwrite)

  cat = np.loadtxt(args.catname)   #read the input catalogue
  if(args.sepname != None):  #if a file with the criteria for the separation is given
    d, hasz = np.loadtxt(args.sepname, usecols=args.columns).T
  else:   #if not keep it from the catalogue and resize the catalogue taking out the last two columns
    d, hasz = cat[:,args.columns].T
    cat = cat[:,:3]   #drop the last two columns
  if(args.verbose == True):
    print("Catalogue read")

  d1 = d==1    #index of all the elements in catalogue 1
  d2 = (d==2)&(hasz==1)   #list of all the elements in catalogue 2 with redshift

  if(args.verbose == True):
    print("Saving the output file in '{}' and '{}'".format(ofiled1, ofiled2))
  #print out the file 
  np.savetxt(ofiled1, cat[d1,:], fmt=args.fmt, delimiter="\t")
  np.savetxt(ofiled2, cat[d2,:], fmt=args.fmt, delimiter="\t")

  exit()
