#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Converts catalogues with ra, dec, redshift into positions x y z in Mpc/h assuming a cosmology

The output file has the structure
x	y	z	w	bias	n(x,y,z)	n_b(z)	M	redfhist
(The columns that are not present in the input file are filled with 1, except the redshift

"""

#import constants as const
import common_functions as cf
import cosmologydir as c
import my_functions as mf
import numpy as np
import scipy.interpolate as spi
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
  
  import argparse as ap  #import optsparse: allows nice command line option handling
  import argparse_custom as apc

  description = """ Convert the file from ra, dec, redshift into cartesian coordinates assuming a cosmology.
    The name of the output file name is derived from the input one and the modification 
    can be set with the options. ra, dec and redshift are by default assumed to be in the first three columns,
    and all the following columns are by default copied after x, y and z. If the number of columns in 
    the input file (or in the ones read) is more than 8 the exceding ones are cut.
    Redshift is copied in the last column. 
    The output file has the structure
    x	y	z	w	bias	n(x,y,z)	n_b(z)	M	redshift
    (The columns that are not present in the input file are filled with 1, except the redshift)"""
  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'), 
      help="Input file name(s), containing ra and dec in the first two columns")

  p = apc.version_verbose( p, '1' )

  p, group = apc.insert_or_replace1(p)
  p, group = apc.overwrite_or_skip(p)

  description="""Cosmology to use to convert ra, dec and z in distances. 
  h0=1 is equivalent to have the distance in Mpc/h with any other value of h0"""
  p, cosmo = apc.cosmology_group( p, description=description, h0_def=1. )
  cosmo.add_argument("--zrange", action="store", nargs=2, type=float,
      help="""Lower and upper limit for the redshift. If this option is given
      the distance is computed only ones and then interpolated in the values in the files""")
  cosmo.add_argument("--nbins", action="store", type=int, default='500', 
      help='Number of bins in redshift.')

  p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+', 
      help="Format of the output files")

  p.add_argument("--usecols", action="store", nargs="+", type=int,
      help="""Read the selected columns. By default read all the columns.If thi option is used, 
      make sure the the first three columns read in are ra, dec and redshift.""") 

  return p.parse_args(args=argv)


def rdz2xyz(rdz, dis):
  """Get an array with ra, dec and redshift and return an array with x, y and z in Mpc/h 

  ----------
  Parameters
  rdz: 2D numpy array
    array with ra, dec, redshift
  dis: function
    function that returns the distance at given redshifts

  Output
  xyz: 2D numpy array
    comoving coordinates: x, y, z
  """
  xyz = np.empty_like(rdz)  #create the output array
  rdz[:,0] *= np.pi/180     #convert ra from deg to rad
  rdz[:,1] = np.pi/2. - rdz[:,1] * np.pi/180     #convert dec from deg to rad
  rdz[:,2] = dis( rdz[:,2] )       #transform redshift in distance in Mpc/h
  xyz[:,0] = rdz[:,2] * np.sin(rdz[:,1]) * np.cos(rdz[:,0])   #fill the x coordinates
  xyz[:,1] = rdz[:,2] * np.sin(rdz[:,1]) * np.sin(rdz[:,0])   #fill the y coordinates
  xyz[:,2] = rdz[:,2] * np.cos(rdz[:,1])    #fill the z coordinates

  return xyz

if __name__ == "__main__":   # if is the main

  args = parse(sys.argv[1:])

  cos = cf.set_cosmology( args.om, args.ok, args.wde, h0=args.h0 ) #set the cosmology
  d = c.distance.Distance( cos )   #create the distance object
  if( args.zrange == None ):  #If no interpolation required
    dis = d.comoving_distance_z #reference the function that compute the distance
  else:   #otherwise compute the comoving distance 
    z = np.linspace( args.zrange[0], args.zrange[1], num=args.nbins ) #find the redshifts
    D_c = d.comoving_distance_z(z)    #get the distance
    dis = spi.UnivariateSpline( z, D_c )   #create the spline function

  absmin, absmax= np.array([+np.inf,]*3), np.array([-np.inf,]*3)  #initialise the max and minimum 

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

    cat = np.loadtxt( fn.name, usecols=args.usecols )  #read the input catalogu

    out = np.ones((cat.shape[0], 9))   #create the output catalogue

    out[:,8] = cat[:,2]   #save the redshift in the last column of the output file
    if(cat.shape[1]<=8):  #if the total number of columns is less than 8
      out[ :, 3:cat.shape[1] ] = cat[ :, 3:]   #all are copied
    else:
      out[ :, 3:8 ] = cat[ :, 3:8 ]   #only the first 5 columns after ra, dec, redshift are copied

    out[:,:3] = rdz2xyz(np.copy(cat[:,:3]), dis)   #convert ra, dec, red in x,y,z in Mpc/h
    #find the absolute max and minimum
    pmin, pmax = np.amin(out[:,:3], axis=0), np.amax(out[:,:3], axis=0)
    absmin, absmax = np.amin([absmin,pmin],axis=0), np.amax([absmax,pmax],axis=0)   #get absolute maximum and minimum

    np.savetxt(ofile, out, fmt=args.fmt, delimiter='\t')


  print("Absolute maximum and minimum:", absmax, absmin)
  print("Difference:", absmax-absmin)


  exit()
