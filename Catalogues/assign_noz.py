#!/usr/bin/python
# -*- coding: utf-8 -*-
#Assign n(z) read from a file

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
  %prog [options] n_z_filename catalogues_filenames
  Given a list of catalogues ('catalogues_filenames') which contain the redshift in one 
  column [the last by default] and a file containing n(z) ('n_z_filename'), 
  it add n(z) in the catalogues and save them.
  The n(z) computed for each file are by default averaged and plotted on the screen.
  Through the options it is possible to save them file, as well as decide if plot or save 
  the average or the single values.
  """)

  p.add_option("-v", "--verbose", action="store_true", dest="verb", help="Verbose")

  p.add_option("-z", "--z-column", action="store", dest="zcol", type=int, default="-1", help="Column in the catalogue which contains the redshift. [Default: %default]")
  p.add_option("-n", "--nz-column", action="store", dest="nzcol", type=int, default="5", help="Column in the catalogue in which n(z) is to be saved. [Default: %default]")
  p.add_option("--zn-file-col", action="store", dest="zn_f_col", type=int, nargs=2, default=[0,1], 
      help="""Columns in 'n_z_filename' containing z and n(z). [Default: %default]""")

  p.add_option("-o", "--outputsub", action="store", dest="outsub", nargs=2, default=["", ".out"],
               help="""Modify the input file name and use it as output file. 
	       If outsub[0]="", output[1] added to the end of the input file name,
	       if outsub[1]="", output[0] added at the beginning of the input file name,
	       else substitute outsub[0] with output[1]. [Default: %default]""")
  p.add_option("-f", "--force-overwrite", action="store_true", default=False, dest="force", help="It forces to overwrite the input file")

  return p.parse_args()

if __name__ == "__main__":   #if it's the main

  opt, args = options(op.OptionParser(version="%prog version 0.1"))   #create the object optparse

  #read n(z)
  z, noz = np.loadtxt(args[0], usecols=opt.zn_f_col).T

  #loop through the catalogues and add a columns with n(z)
  for fn in args[1:]:
    if(opt.verb):
      print("Processing file:", fn)
    cat = np.loadtxt(fn)

    cat[:,opt.nzcol] = np.interp(cat[:,opt.zcol], z, noz, left=0, right=0)
  
    #create the output file name
    if(opt.force == True):
      ofile = fn
    else:
      if(opt.outsub[0] == ""):
	ofile = fn+opt.outsub[1]
      elif(opt.outsub[1] == ""):
	ofile = opt.outsub[0]+fn
      else:
	ofile = fn.replace(opt.outsub[0], opt.outsub[1])
      if(ofile == fn):   #if the substitution didn't work add outsub[1] at the end of the file name
	ofile = fn+".out"

    np.savetxt(ofile, cat, fmt="%8.7e", delimiter='\t')

  exit()
