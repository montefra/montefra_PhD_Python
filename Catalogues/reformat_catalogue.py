#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Converts files with a given structure 
to files with structure
x-xmin	y-ymin	z-zmin	weight*(fiber*comp)=1	bias=1	n(z)	n_b(z)=n(z)	M=1	redfhift=1
"""

import glob   #allows search of files with unix wildcards
import numpy as np
import optparse as op  #import optsparse: allows nice command line option handling

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
	      Reorganize the columns of a given catalogue in order to match my power spectrum code
	      'file_names' can be a list and accept also shell wildcards. The name of the output file 
	      name is derived from the input one and the modification can be set with the options
	      The output file has the structure
	      x	y	z	w	bias	n(z)	n_b(z)	M	redshift  (not used now)
	      x	y	z	w(FKP)	bias	n(z)	w_sys	M	redshift
	      (The exact repositioning of the columns from the input to the output can be adjusted 
	      modifying function 'reorganise_cols')""")
 
  p.add_option("-v", "--verbose", action="store_true", dest="verb", help="Verbose")

  p.add_option("--delimiter", action="store", dest="deli", 
               help="The string used to separate values. By default, this is any whitespace.")
  p.add_option("--skiprow", action="store", dest="skip", type=int, default=0, help="Skip the first 'skiprows' lines. [Default: %default]")

  p.add_option("-r", "--rescale", action="store", type=float, nargs=3, dest="resc",
               help="""If set use those values to rescale the x,y,z otherwise no riscaling done""")

  p.add_option("-o", "--outputsub", action="store", dest="outsub", nargs=2, default=["", ".out"],
               help="""Modify the input file name and use it as output file. 
	       If outsub[0]="", output[1] added to the end of the input file name,
	       if outsub[1]="", output[0] added at the beginning of the input file name,
	       else substitute outsub[0] with output[1]. [Default: %default]""")
  p.add_option("-f", "--force-overwrite", action="store_true", default=False, dest="force", help="It forces to overwrite the input file")

  return p.parse_args()

def reorganise_cols(cat, resc=None):
  """Takes catalogue 'cat' and returns the catalogue to be printed out.
  Parameters
  cat: 2D numpy array
    input catalogue
  resc: list of 3 float (optional)
    if given and of dimension 3 the values are used to rescale the positions of the particles in the catalogues. not checked

  Output
  ocat: 2D numpy array
    output catalogue
  """

  ocat = np.ones((cat.shape[0], 9))   #create the output catalogue

  #all the catalogues I've used have the positionx x,y,z as first three columns
  ocat[:,:3] = cat[:,:3]
  if(resc != None):   #if rescaling required, rescaling applied
    ocat[:,:3] = ocat[:,:3] - resc + 0.01

  #will catalogues adapted by ariel
  #ocat = np.hstack((ocat, np.ones_like(ocat[:,:2])))   #add two columns of one
  #ocat = np.vstack((ocat.T, cat[:,8], np.ones_like(ocat[:,2]))).T   #add weights and a column of one
  #ocat = np.vstack((ocat.T, cat[:,5], cat[:,5])).T    #add twice n(z)
  #ocat = np.hstack((ocat, np.ones_like(ocat[:,:2])))   #add two columns of one

  #ariel own catalogues
  #ocat[:,3] = cat[:,9] *cat[:,11]    #copy the weights
  #ocat[:,4] = cat[:,5]
  #ocat[:,5] = cat[:,6]*cat[:,8]   #n(z)*compl
  #ocat[:,6] = cat[:,9] #*cat[:,11]    #copy the weights
  #ocat[:,6] = cat[:,6]*cat[:,7]   #copy the weights
  #ocat[:,6] = cat[:,7]    #copy the weights
  #ocat[:,8] = cat[:,4]    #copy the redshift

  #reconstructed mock catalogues with redshifts in the 5th column
  ocat[:,3] = cat[:,3]    #FKP weights
  ocat[:,8] = cat[:,4]    #redshift

  #reconstructed catalogue, save fpk weights in the 8th column and weiths set to 1
  #ocat[:,7] = cat[:,3]    #fkp weights
  #ocat[:,8] = cat[:,8]    #redshift


  return ocat

if __name__ == "__main__":   # if is the main

  opt, args = options(op.OptionParser(version="%prog version 1.0"))

  absmin, absmax= np.array([+np.inf,]*3), np.array([-np.inf,]*3)  #initialise the max and minimum 
  for arg in args:   #loop through the arguments
    for fn in glob.iglob(arg):   #loop to the file names associated to each argument. This allows for shell wildcard

      if(opt.verb):
	print("Processing file:", fn)
      ida = np.loadtxt(fn, skiprows=opt.skip)   #read the file

      oda = reorganise_cols(ida, resc=opt.resc)   #reforamt the catalogue

      #find the maximum and minimum in each catalogue
      pmin, pmax = np.amin(oda[:,:3], axis=0), np.amax(oda[:,:3], axis=0)
      if(opt.verb):
	print(fn+". maximum and minimum:", pmax, pmin)
      absmin, absmax = np.amin([absmin,pmin],axis=0), np.amax([absmax,pmax],axis=0)   #get absolute maximum and minimum

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

      np.savetxt(ofile, oda, fmt="%8.7e", delimiter='\t')


  print("Absolute maximum and minimum:", absmax, absmin)
  print("Difference:", absmax-absmin)

exit()





