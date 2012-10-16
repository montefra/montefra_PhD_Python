#!/usr/bin/python
# -*- coding: utf-8 -*-
"""read a file with k, P(k), P(k)+shotnoise, rescale the shot noise by the given fraction and save P'(k) instead of P(k)"""

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
  
  description = """ Given a file (list) containing power spectra with structure k, P(k), P(K)+shotnoise,
  and a fraction 'f', saves a file with k, P(k)-f*shotnoise """

  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("fraction", action="store", type=float, help="Fraction for rescaling the shot noise")
  p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'), 
      help="Input file name(s), containing k, P(k) and P(k)+shotnoise (maybe something else).")

  p = apc.version_verbose( p, '0.1' )

  p, group = apc.insert_or_replace1(p, print_def=True)
  p, group = apc.overwrite_or_skip(p)

  p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, 
      help="Format of the output files. (default: %(default)s)")

  return p.parse_args(args=argv)

if __name__ == "__main__":   #if it's the main

  args = parse(sys.argv[1:])

  for fn in args.ifname:  #file name loop
    if(args.verbose == True):
      print("Process file '{0}'.".format(fn.name))

    #create the output file name and check it
    if(args.replace == None):
      ofile, skip = mf.insert_src(fn.name, args.insert, overwrite=args.overwrite, skip=args.skip)
    else:
      ofile, skip = mf.replace_src(fn.name, args.replace, overwrite=args.overwrite, skip=args.skip)
    if(skip == True):
      print("Skipping")
      continue

    pk = np.loadtxt( fn )  #load the file
    noise = (pk[:,2]-pk[:,1]).mean()   #get the shot noise
    pk[:,1] -= args.fraction*noise  #decrease the power spectrum by the fraction of the shot noise

    np.savetxt( ofile, pk, fmt=args.fmt, delimiter="\t" )
  #endfor file name loop

  exit()
