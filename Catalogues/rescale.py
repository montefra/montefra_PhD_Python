#!/usr/bin/python
# -*- coding: utf-8 -*-
"""rescale the x, y and z positions of object in a catalogue by a given amount """

#import constants as const
import glob   #allows search of files with unix wildcards
import numpy as np
import optparse as op  #import optsparse: allows nice command line option handling
import sys


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
	      Subtract from the positions x,y and z of the object in catalogues 
	      the desired values, in order to rescale the positions.
	      They must be in the first three columns of the input file.
	      'file_names' can be a list and accept also shell wildcards.""")

  p.add_option("-v", "--verbose", action="store_true", dest="verb", help="Verbose")

  p.add_option("-r", "--rescale", action="store", type=float, nargs=3, dest="resc",
               help="""If set use those values to rescale the x,y,z otherwise compute absolute 
minimum from all the input files and use that value""")

  p.add_option("--comments", action="store", dest="comm", default="#",
	       help="The character used to indicate the start of a comment. [Default: %default]")
  p.add_option("--delimiter", action="store", dest="deli", 
               help="The string used to separate values. By default, this is any whitespace.")
  p.add_option("--skiprow", action="store", dest="skip", type=int, default=0, help="Skip the first 'skiprows' lines. [Default: %default]")

  p.add_option("-o", "--outputsub", action="store", dest="outsub", nargs=2, default=["", "out"],
               help="""Modify the input file name and use it as output file. If outsub[0]="", 
output[1] added to the end of the input file name, if outsub[1]="", output[0] added at the beginning of the input file name, else substitute outsub[0] with output[1]. [Default: %default]""")
  p.add_option("-f", "--force-overwrite", action="store_true", default=False, dest="force", help="It forces to overwrite the input file")

  return p.parse_args()

if __name__ == "__main__":   # if is the main

  opt, args = options(op.OptionParser(version="%prog version 0.1"))

  #check the parameters
  if(opt.outsub[0] == opt.outsub[1]):
    print("The input and output file must have different names")
    sys.exit(10)

  #search for the absolute minimum if not given
  if(opt.resc == None):
    print("Finding the absolute maximum and minimum")
    absmin = np.array([+np.inf,]*3)
    for arg in args:   #loop through the arguments
      for fn in glob.iglob(arg):   #loop to the file names associated to each argument. This allows for shell wildcard
	cat = np.loadtxt(fn, comments=opt.comm, delimiter=opt.deli, skiprows=opt.skip, usecols=(0,1,2))
	pmin = np.amin(cat, axis=0)
	absmin = np.amin([absmin,pmin],axis=0) #get absolute maximum and minimum

    print("absolute minimum:", absmin)
    opt.resc=absmin

  #rescale
  for arg in args:   #loop through the arguments
    for fn in glob.iglob(arg):   #loop to the file names associated to each argument. This allows for shell wildcard
      if(opt.verb):
	print("Processing file:", fn)
      cat = np.loadtxt(fn, comments=opt.comm, delimiter=opt.deli, skiprows=opt.skip)

      cat[:,:3] = cat[:,:3] - opt.resc + 0.01

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

      np.savetxt(ofile, cat, fmt="%8.7e", delimiter='\t')  #save the new file

  exit()
