#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Given a set of catalogues with x, y and z in the first three columns,
compute the redshift corresponding to the objects and write the output 
in an extra column or in an existing one"""

import argparse as ap
import argparse_custom as apc
import common_functions as cf
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
  
  description = """Given a set of catalogues with x, y and z in the first three columns,
    compute the redshift corresponding to the objects and write the output 
    in an extra column or in an existing one"""
  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'), 
      help="Input file name(s), containing ra and dec in the first two columns")

  p = apc.version_verbose( p, '1' )

  p, group = apc.insert_or_replace1(p)
  p, group = apc.overwrite_or_skip(p)

  p.add_argument("-z", "--z-col", type=int, help="""If given the redshifts are substituted in 
      '%(dest)s' column. Otherwise an extra column is appended in the output file""")

  p.add_argument("-f", "--zdist-file", type=ap.FileType('r'), 
      help="""File containing a table with resdshift and comoving distance.
      If not given, this table is computed using the cosmological parameters
      that can be set with the options.""")

  description="""If the file with the table of z and D_c(z) is not given, 
  this is compute using the following parameters"""
  p, cosmo = apc.cosmology_group( p, description=description )
  cosmo.add_argument("--zrange", type=float, nargs=2, action="store", default=[0.4, 0.8],
      help="Maximum and minimum redshift where to compute the comoving distance")
  cosmo.add_argument("--nbins", action="store", type=int, default='500', 
      help='Number of bins in redshift.')

  p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+', 
      help="Format of the output files")

  return p.parse_args(args=argv)

if __name__ == "__main__":   #if it's the main

  args = parse(sys.argv[1:])

  if( args.zdist_file != None ):  #read the file with D_c(z)
    z, D_c = np.loadtxt( args.zdist_file ).T
  else:   #or compute it
    z = np.linspace( args.zrange[0], args.zrange[1], num=args.nbins )
    import cosmologydir as c
    cos = cf.set_cosmology( args.om, args.ok, args.wde ) #set the cosmology
    #set the distance object
    dis = c.distance.Distance(cos)
    #compute the table with D_c(z) 
    if(args.h0 == None):  
      D_c = dis.comoving_distance_zh(z)   #Mpc/h
    else:    
      D_c = dis.comoving_distance_z(z)    #Mpc

  #initialise the smooth spline fit to optain z at given D_c
  z_Dc = spi.UnivariateSpline( D_c, z )

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
    
    cat = np.loadtxt( fn )  
    if(args.verbose == True):
      print( "Catalogue read" )

    #compute comoving distance of objects from the center
    dist = np.sqrt( cat[:,0]**2 + cat[:,1]**2 + cat[:,2]**2 ) 
    #get the redshift
    z = z_Dc( dist ) 
    if(args.verbose == True):
      print( "Redshift computed" )

    if( args.z_col == None ):  #if no column given z appended to the catalogue
      cat = np.vstack( (cat.T, z) ).T 
    else:  #otherwise substitute z in the desired column
      cat[:, args.z_col ] = z

    np.savetxt( ofile, cat, delimiter="\t", fmt=args.fmt )
    if(args.verbose == True):
      print( "File '{0}' saved".format(ofile) )

  #end file name loop

  exit()
