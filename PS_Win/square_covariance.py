#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse as ap
import argparse_custom as apc
import my_functions as mf
import numpy as np
import sys

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
  
  description = """Given a covariance matrix in a column of the input file(s), make it square"""

  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'), 
      help="Input file name(s), containing ra and dec in the first two columns")

  p = apc.version_verbose( p, '1' )

  p, group = apc.insert_or_replace1(p)
  p, group = apc.overwrite_or_skip(p)

  p.add_argument("-c", "--column", action="store", type=int, default=4, help="Column number in the input file containing the covariance.")

  return p.parse_args(args=argv)


if __name__ == "__main__":   # if is the main

  args = parse(sys.argv[1:])

  for fn in args.ifname:  #file name loop
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

    cov = np.loadtxt(fn, usecols=[args.column])   #read the covariance

    dimcov = int(np.sqrt(cov.size))  #size of the output covariance matrix

    cov = cov.reshape([dimcov,dimcov])   #reshape to a square

    np.savetxt(ofile, cov, fmt="%8.7e", delimiter="\t")   #save the matrix
    if(args.verbose == True):
      print( "File '{0}' saved".format(ofile) )

  exit()
