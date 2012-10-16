#!/usr/bin/python
# -*- coding: utf-8 -*-

import glob
import numpy as np
import optparse as op


def options(p):
  """
  This function accept the option parser istance and fill it with options.

  The version must be set when creating the optparse istance

  Parameters
  ----------
  p: option parser istance

  output: tuple with the options and the arguments
  ---------
  """

  p.set_usage("""
  %prog [options] file_names 
  Cut the content of the input file names (also with shell-like wild card)
  according to the given minimum and maximum range in the first column of the file.
  The output file name is determined from the input one through replacement""")

  p.add_option("--min", action="store", type=float, dest="kmin", default=0.001, help="Minimum value of k used to cut. [Default: %default]")
  p.add_option("--max", action="store", type=float, dest="kmax", default=0.3, help="Maximum value of k used to cut. [Default: %default]")

  p.add_option("-r", "--replace", action="store", nargs=2, dest="repl", default=[".dat", ".cut.dat"], help="Replace repl[0] with repl[1] in the input file names, in order to get the output file names. [Default: %default]")

  return p.parse_args()


if __name__ == "__main__":   # if is the main

  (opt, args) = options(op.OptionParser(version="%prog version 0.1"))   #create the object optparse

  for a in args:   #go through the arguments
    for f in glob.iglob(a):    #get all the files names

      ps = np.loadtxt(f)   #read the input file
      
      np.savetxt(f.replace(opt.repl[0],opt.repl[1]), np.compress((ps[:,0]>opt.kmin)&(ps[:,0]<opt.kmax), ps, axis=0), fmt='%8.7e', delimiter='\t')

  exit()

