#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Read one or more files and extract a part of the file according to the constraints given for a column"""

import argparse as ap
import argparse_custom as apc
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
  
  description = """  Read one or more files and extract a part of the file according to the constraints given for a column.
  Examples: 
    %(prog)s 3 'col < 4' filename
    %(prog)s 3 '(col < 4) | (col > 8)' filename
    %(prog)s 3 '(col > 8) & (col < 4)' filename
  'col' is the internal name of the column to be checked: use this can in the 'constr' string.
  WARNING: The constraint is evaluated with 'eval' and no check is done upon local or global variables."""

  p = ap.ArgumentParser(description=description, formatter_class=ap.RawDescriptionHelpFormatter)

  p.add_argument("column", action="store", type=int, help="Column to check")
  p.add_argument("constr", action="store", 
      help="""Constraints to be applied to the desired column""")
  p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'), 
      help="Input file name(s), containing ra and dec in the first two columns")

  p = apc.version_verbose( p, '0.1' )

  p, group = apc.insert_or_replace1(p, print_def=True)
  p, group = apc.overwrite_or_skip(p)

  p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+', 
      help="Format of the output files. (default: %(default)s)")

  return p.parse_args(args=argv)

if __name__ == "__main__":   #if it's the main

  args = parse(sys.argv[1:])

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

    cat = np.loadtxt( fn.name )
    if(args.verbose == True):
      print( "Catalogue read" )

    col = cat[:, args.column]
    col = eval( args.constr )

    np.savetxt( ofile, cat[col, :], delimiter="\t", fmt=args.fmt )
    if(args.verbose == True):
      print( "File '{0}' saved".format(ofile) )


  exit()
