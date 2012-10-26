#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Converts catalogues with ra, dec, redshift into positions x y z in Mpc/h assuming a cosmology

The output file has the structure
x	y	z	w	bias	n(x,y,z)	n_b(z)	M	redfhist
(The columns that are not present in the input file are filled with 1, except the redshift

"""

#import constants as const
import my_functions as mf
import numpy as np

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
  
  import argparse as ap 
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
  p = ap.ArgumentParser(description=description,
      formatter_class=ap.ArgumentDefaultsHelpFormatter)

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
      help="""Read the selected columns. By default read all the columns.If thi
      option is used, make sure the the first three columns read in are ra, dec
      and redshift.""") 

  description = """Parameters related to the parallel computation"""
  p, parallel = apc.parallel_group( p, description=description )

  return p.parse_args(args=argv)  #end def parse(argv)

def comoving_distance( om, ok, wde, h0, zrange=None, nbins=None):
  """Compute the comoving distance given the cosmology and returns a function
  reference. If *zrange* is not *None*, the the comoving distance is evaluated
  in *nbins* in *zrange* and then a UnivariateSpline is returned
  Parameters
  ----------
  om, ok, wde, h0: floats
    omega matter, omega curvature, dark energy equation of state and reduced
    hubble parameter (h0=1 is equivalent of getting the distance in Mpc/h) 
  zrange: 2 element list
    lower and upper limit in the redshift range
  nbins: int
    number of bins in redshift 
  output
  ------
  dis: function
    function that evaluates the comoving distance at given redshift(s)

  Examples
  --------
    d = distance( 0.27, 0, -1, 1)
    z = np.linspace(0, 1, num=50)
    comdis = d(z)
  """

  import cosmologydir as c
  cos = c.Setcosmology( om=om, ok=ok, wde=wde, h=h0 ) #set the cosmology
  d = c.Distance( cos )   #create the distance object
  if( args.zrange == None ):  #If no interpolation required
    dis = d.comoving_distance_z #reference the function that compute the distance
  else:   #otherwise compute the comoving distance 
    import scipy.interpolate as spi  #scipy interpolation routines
    z = np.linspace( zrange[0], zrange[1], num=nbins ) #find the redshifts
    D_c = d.comoving_distance_z(z)    #get the distance
    dis = spi.UnivariateSpline( z, D_c )   #create the spline function
  return dis   #end def comoving_distance( ... )

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

  return xyz   #end def rdz2xyz(rdz, dis):

def convert_save(f, distance, **kwargs ):
  """
  Read file *f*, converts ra, dec, z into cartesian coordinates, computing the
  comoving distance at redshift z stored in *distance*, and save to a new file
  Parameters
  ----------
  f: file object or string
    file containing ra, dec and z
  distance: function
    function that evaluates the comoving distance at given redshift(s)
  kwargs: keyword arguments
  output
  ------
  max, min: lists
    maximum and minimum values of x, y and z
  If kwargs['skip'] == True and the output file name already exists, a *None*
  is returned

  accepted kwargs that affects the function
  +verbose: verbose mode [True|False] 
  +replace: replace string *replace[0]* with *replace[1]* in f.name
  +insert: insert string *insert[0]* before *insert[1]* in f.name
  +skip: existing file names skipped [True|False]
  +overwrite: existing file names overwritten [True|False]
  +usecols: columns to read from the input files. the first three must be ra,
    dec and redshift

  """
  if( type(f) == file ):  #if f is a file object
    fname = f.name  #get the file name
  else:  #it's alread the file name
    fname = f  #if

  if(kwargs['verbose'] == True):
    print("Process catalogue '{0}'.".format(fname))

  #create the output file name and check it
  if(kwargs['replace'] == None):
    ofile, skip = mf.insert_src(fname, kwargs['insert'],
	overwrite=kwargs['overwrite'], skip=kwargs['skip'])
  else:
    ofile, skip = mf.replace_src(fname, kwargs['replace'],
	overwrite=kwargs['overwrite'], skip=kwargs['skip'])
  if(skip == True):
    print("Skipping file '{0}'".format(fname))
    return None

  cat = np.loadtxt( f, usecols=kwargs['usecols'] )  #read the input catalogu

  out = np.ones((cat.shape[0], 9))   #create the output catalogue

  out[:,8] = cat[:,2]   #save the redshift in the last column of the output file
  if(cat.shape[1]<=8):  #if the total number of columns is less than 8
    out[ :, 3:cat.shape[1] ] = cat[ :, 3:]   #all are copied
  else:#only the first 5 columns after ra, dec, redshift are copied
    out[ :, 3:8 ] = cat[ :, 3:8 ]   

  out[:,:3] = rdz2xyz(np.copy(cat[:,:3]), distance)   #convert ra, dec, red in x,y,z in Mpc/h

  #save the converted catalogue
  np.savetxt(ofile, out, fmt=kwargs['fmt'], delimiter='\t')

  #return the max and minimum
  return np.amax(out[:,:3], axis=0), np.amin(out[:,:3], axis=0)
  #end def convert_save(f, distance, **kwargs ):


if __name__ == "__main__":   # if is the main

  import sys
  args = parse(sys.argv[1:])

  #compute the comoving distance for the given cosmology
  dis = comoving_distance( args.om, args.ok, args.wde, args.h0,
      zrange=args.zrange, nbins=args.nbins)

  #if parallel computation required, check that Ipython.parallel.Client 
  #is in installed and that the ipycluster has been started
  if args.parallel :

    import ipython_parallel as IPp
    import os
    #the absolute path and file name of this script
    path, fname = os.path.split( os.path.abspath(sys.argv[0]) )
    function_name = 'rdz2xyz'  #name of the function to import
    #command to run on all the engines
    imports = [ 'import numpy as np', 'import my_functions as mf',
	#add the script directory to the python path
	'import sys', 'sys.path.append("{0}")'.format(path),     
	#import the desired function in the namespace
	'from {0} import {1}'.format( os.path.splitext(fname)[0], function_name) ]  
    args.parallel, lview = IPp.start_load_balanced_view( to_execute=imports )

  #run the script in serial mode
  if( args.parallel == False ):  #if: parallel
    #initialise the list of maxima and minima in the output file
    maxi, mini = [], []
    for fn in args.ifname:  #file name loop
      #convert the coordinates and return maxima and minima
      temp = convert_save(fn, dis, **vars(args) ) 
      if( temp != None ):
	maxi.append( temp[0] )
	mini.append( temp[1] )
  #run the script using the IPython parallel environment 
  else:    #if: parallel
    engines_id = lview.client.ids  #get the id of the engines_id
    initstatus = lview.queue_status()  #get the initial status

    #submit the jobs and save the list of jobs
    runs = [ lview.apply( convert_save, os.path.abspath(fn.name), dis, **vars(args) ) 
	for fn in args.ifname ]

    if args.verbose :   #if some info is required
      IPp.advancement_jobs( lview, runs, engines_id, update=args.update,
	  init_status=initstatus )
    else:   #if no info at all is wanted
      lview.wait( jobs=runs )  #wait for the end

    #get the maxima and minima from the computations excluding the None
    maxi = [r.result[0] for r in runs if r is not None]
    mini = [r.result[1] for r in runs if r is not None]
  #end if: parallel

  #compute absolute maximum and minimum
  absmax = np.max( maxi, axis=0)
  absmin = np.min( mini, axis=0)

  maxstring_lenght = len("Difference" )
  string_template = "{:<{}}  "+"  ".join(["{:^9}",]*3)
  float_template = "{:<{}}: "+", ".join(["{:>+9.3f}",]*3)

  print string_template.format( " ", maxstring_lenght, "x", "y", "z" )
  print float_template.format( "Maximun", maxstring_lenght, *absmax )
  print float_template.format( "Minimum", maxstring_lenght, *absmin )
  print float_template.format( "Difference", maxstring_lenght, *absmax-absmin )

  exit()
