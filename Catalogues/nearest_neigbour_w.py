#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Read the las two columns of a file with 'ra dec z sample has_z nearneighbour' and
the corresponding file wit 'x, y, z, w, b, n, nb, L, red', delete all the entries with has_z=0 and 
upweights the nearest neighbour
"""

import argparse as ap
import common_functions as cf
import healpy as H
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
  
  description = """The nearest galaxies to the ones which have no redshift are upweighted
  The second to last column of the 'ifradez' catalogue contains 1 if the redshift has to be kept, 0 if it has to be 
  substituted with the nearest neigbour one, whose index is in the last column. The output file has the same structure
  of 'ifxyz': the objects without redshifts are dropped and the nearest neibour upweighted"""
  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("ifradez", action="store", type=ap.FileType('r'), 
      help="""Input file name, containing 1 or 0 in the second to last ones. -1 or index of the nearest neighbour in the last column.""")
  p.add_argument("ifxyz", action="store", type=ap.FileType('r'), 
      help="""Input file name, containing the catalogue with a column to upweight according to the nearest neighbour correction""")

  p.add_argument('--version', action='version', version='0.1')
  p.add_argument("-v", "--verbose", action="store_true")

  ir_group = p.add_mutually_exclusive_group()
  ir_group.add_argument("-i", "--insert", action="store", nargs=2, default=[".d1d2", ".dat"], 
      help="""Output file name created inserting '%(dest)s[0]' before '%(dest)s[1]' in 'ifxyz'.
      Only the first three columns of the file are saved""")
  ir_group.add_argument("-r", "--replace", action="store", nargs=2, 
      help="""Output file name created replacing '%(dest)s[0]' with '%(dest)s[1]' in 'ifxyz'.
      Only the first three columns of the file are saved. If none of these two options given, 'insert' is assumed""")

  group = p.add_mutually_exclusive_group()
  group.add_argument("--overwrite", action="store_true", 
      help="If given does not check if the output file exists or not before saving it.")
  group.add_argument("--skip", action="store_true", 
      help="Skip already existing output file names. No operation done on the corresponding input file.")

  p.add_argument("-w", "--w-column", action="store", type=int, default=3,
      help="Column containing the weight to increase by 1 for each nearest neighbour without redshift")
  p.add_argument("-s", "--s-columns", action="store", type=int, nargs=2, default=(-2,-1),
      help="Columns containing the information about which ojects has redshift and which is the nearest neighbour")

  p.add_argument("--fmt", default="%7.6e", action="store", nargs='+', help="Format of the output files")

  return p.parse_args(args=argv)



if __name__ == "__main__":   #if it's the main

  args = parse(sys.argv[1:])

  #create the output file name
  if(args.replace == None):
    ofile, skip = mf.insert_src(args.ifxyz.name, args.insert, overwrite=args.overwrite, skip=args.skip)
  else:
    ofile, skip = mf.replace_src(args.ifxyz.name, args.replace, overwrite=args.overwrite, skip=args.skip)
  if(skip == True):
    print("Skipping")
    exit()

  if(args.verbose == True):
    print("Reading catalogues {0} and {1}".format(args.ifradez.name, args.ifxyz.name))
  has_z, nn = np.loadtxt(args.ifradez, dtype='int', usecols=args.s_columns).T
  cat = np.loadtxt( args.ifxyz )
  if(args.verbose == True):
    print("Catalogues read")

  for n in nn[has_z==0]:   #loop though the objects without redshift
    cat[n,args.w_column] += 1   #increase by one the nearest neighbour

  if(args.verbose == True):
    print("Saving the output file in '{}'".format(ofile))
  #print out the file 
  of = open(ofile, 'w')
  np.savetxt(of, cat[has_z!=0,:], fmt=args.fmt, delimiter="\t")
  of.close()

  exit()
