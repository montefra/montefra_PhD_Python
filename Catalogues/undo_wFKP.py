#!/usr/bin/python
# -*- coding: utf-8 -*-
"""The weights w are preplaced by w/w_fpk, with w_fkp = 1/(1+ n(z) * pw )"""

import argparse as ap
import argparse_custom as apc
import my_functions as mf
import numpy as np
import scipy.interpolate as spi
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
  
  description = """The weights w are preplaced by w/w_fpk, with w_fkp = 1/(1+ n(z) * pw )"""
  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("pw", action="store", type=float, help="Value of 'pw' to use for the FKP weight")
  p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'), 
      help="Input file name(s), containing ra and dec in the first two columns")

  p = apc.version_verbose( p, '0.1' )

  p, group = apc.insert_or_replace(p)
  p, group = apc.overwrite_or_skip(p)

  p.add_argument("-n", "--n-col", type=int, default=5, 
      help="""Column containin the number density""")
  p.add_argument("-w", "--w-col", type=int, default=3, 
      help="""Column containin the weights""")

  p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+', 
      help="Format of the output files")

  return p.parse_args(args=argv)

if __name__ == "__main__":   #if it's the main

  args = parse(sys.argv[1:])

  means = []

  for fn in args.ifname:
    if(args.verbose == True):
      print("Process catalogue '{0}'.".format(fn.name))

    #create the output file name and check it
    if(args.replace == None):
      ofile, skip = mf.insert_src(fn.name, args.insert, overwrite=args.overwrite, skip=args.skip)
    else:
      ofile, skip = mf.replace_src(fn.name, args.replace, overwrite=args.overwrite, skip=args.skip)
    if(skip == True):
      print("Skipping")
      continue

    cat = np.loadtxt( fn, usecols=(args.w_col, args.n_col) )  
    if(args.verbose == True):
      print( "Catalogue read" )

    #cat[:, args.w_col] *= 1. + cat[:, args.n_col] * args.pw 
    #mean = np.mean(cat[:, args.w_col])
    cat[:, 0] *= 1. + cat[:, 1] * args.pw 
    mean = np.mean(cat[:, 0])
    means.append( mean ) 
    print( mean, np.std( cat[:, 0] ) )
    #print( mean, np.std( cat[:, args.w_col] ) )

    #np.savetxt( ofile, cat, delimiter="\t", fmt=args.fmt )
    if(args.verbose == True):
      print( "File '{0}' saved".format(ofile) )

  print( "mean of the mean: {0}".format(np.mean(means)) )

  exit()
