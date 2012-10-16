#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Separate (a) given catalogue(s) into d1 and d2 sample as in Guo, Zehavi & Zheng 2011"""

import argparse as ap
import common_functions as cf
import healpy as H
import itertools as it
import my_functions as mf
import numpy as np
import sys  #system stuff

headersh = "#ra	dec	z	\n"  #add the header
headerout = "#ra	dec	z	sample	has_z	nearneighbour\n"  #add the header

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
  
  description = """Separate the input catalogues in two as in Guo, Zehavi & Zheng 2011"""
  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'), 
      help="Input file name(s), containing ra and dec in the first two columns")

  p.add_argument('--version', action='version', version='0.1')
  p.add_argument("-v", "--verbose", action="store_true")

  ir_group = p.add_mutually_exclusive_group()
  ir_group.add_argument("-i", "--insert", action="store", nargs=2, default=[".d1d2", ".dat"], 
      help="""Output file name created inserting '%(dest)s[0]' before '%(dest)s[1]' in the input file name.
      The output file is saved appending three columns to the input one:
      first: object belongs to sample d1 ("1") or d2 ("2");
      second: if '0' the object is to be dropped to simulate fiber collision.
      third: if object in d2, the index of the nearest neighbour in d1, otherwise -1.""")
  ir_group.add_argument("-r", "--replace", action="store", nargs=2, 
      help="""Output file name created replacing '%(dest)s[0]' with '%(dest)s[1]' in the input file name.
      The output file is saved appending three columns to the input one:
      first: object belongs to sample d1 ("1") or d2 ("2");
      second: if '0' the object is to be dropped to simulate fiber collision.
      third: if object in d2, the index of the nearest neighbour in d1, otherwise -1.
      If none of these two options given, 'insert' is assumed""")

  group = p.add_mutually_exclusive_group()
  group.add_argument("--overwrite", action="store_true", 
      help="If given does not check if the output file exists or not before saving it.")
  group.add_argument("--skip", action="store_true", 
      help="Skip already existing output file names. No operation done on the corresponding input file.")

  shuffle = p.add_argument_group(description="Shuffle the catalogue(s) before processing it(them)")
  shuffle.add_argument('--shuffle', action="store", nargs=2, 
      help="""Shuffle the input catalogue(s) before proceding with the separation.
      The shuffles catalogue(s) is saved in a file, whose name is created inserting '%(dest)s[0]' 
      before '%(dest)s[1]'. The output file names are then created on top of this file.""")
  sh_group = shuffle.add_mutually_exclusive_group()
  sh_group.add_argument("--shuffle-overwrite", action="store_true", dest="soverwrite",
      help="If given does not check if the output file exists or not before saving it.")
  sh_group.add_argument("--shuffle-skip", action="store_true", dest="sskip",
      help="Skip already existing output file names. No operation done on the corresponding input file.")
  shuffle.add_argument("--shuffle-fmt", default="%7.6e", action="store", dest="sfmt", nargs='+', help="Format of the shuffled files")

  p.add_argument("-p", "--polygons", action="store", 
      default="/data01/montefra/BOSS/Catalogs/Geometry/boss_geometry_2011_06_10.ply", type=ap.FileType('r'), 
      help="Mangle file containing the survey geometry.")
  p.add_argument("-s", "--poly2sect", action="store", 
      default="/data01/montefra/BOSS/Catalogs/Geometry/boss_geometry_100611_secntiles.dat", type=ap.FileType('r'), 
      help="""File containing in the first column the number of the sector to which the polygon 
      identified by the line number belongs and in the second the number of tiles covering the secotr.""")
  p.add_argument('--temp', action="store", type=ap.FileType('w'), default="temp_polyid.txt", 
      help="""Temporary file to store the output of 'polyid'. Make sure to change it when running 
      this program more than once in parallel from the same directory.""")

  compl = p.add_argument_group(description = """File containing the completeness due to fiber collision per sector, 
      indicated by the line number. If not given, the mean value of the fraction of collided galaxies with redshift 
      in sectors with 1, 2, 3 or 4 tiles is used [0.05, 0.76, 0.88, 0.88].""")
  compl.add_argument('--fc-file', action="store", type=ap.FileType('r'), help="File name.")
  compl.add_argument("--fc-column", action="store", type=int, default=0, help="Column containing the completeness")

  p.add_argument("-n", "--nside", action="store", type=cf.pow2, default=256, 
    help="Healpix parameter: related with the number of pixels. Must be power of 2")

  p.add_argument("-f", "--fiber-collision", default=62., action="store", type=float, dest="fdiam", 
    help="Diameter of the fiber for fiber collision in arcseconds.")

  p.add_argument("--header", action="store_true", help="""If selected the the following headers will be 
      written '{0}' for the shuffled file, '{1}' for the output one""".format(headersh, headerout) )
  p.add_argument("--fmt", default="%7.6e", action="store", nargs='+', help="Format of the output files")

  p.add_argument("--israndom", action="store_true", 
      help="If given, the input files are considered random catalogue and are simply downsampled sector by sector")

  return p.parse_args(args=argv)

def separate_d1d2(nside, pixels, vector, diameter, nest=False):
  """Separate the input catalogue stored in vector according to the distance 
  between objects as in Guo, Zehavi & Zheng 2011

  Parameters
  ----------
  nside : int, scalar or array-like
    The healpix nside parameter, must be a power of 2
  pixels : array of int
    The healpix pixel numbers. Scalar if all input are scalar, array otherwise.
    Must have size equal to the first dimention of vector
  vector: array
    Healpix 3D position vector. Must have the first dimention equal to the size of pixels
  diameter: float
    minimum distance between objects in radians
  nest: bool (optional)
    if True, assume NESTED pixel ordering, otherwise, RING pixel ordering

  output
  ------
  d: int or array of int
    1 if the distance of the object is larger than 'diameter' from 
    the other objects with d==1 preceding it in the array, 2 otherwise
  nnid: int or array of int
    for object with d==2 contains the index of the nearest neighbour, for d==1, is -1
  """

  #initialize an empty list with a number of elements given by the maximum
  #number in pixels and fill it with the index of the objects in each pixel
  pointers = [[] for i in xrange(pixels.max()+1)]
  #insert in position 'p' of 'pointers' the indeces of the object that belong to pixel 'p'
  for i,p in enumerate(pixels):
    pointers[p].append(i)
  
  d = 2*np.ones_like(pixels)  #create an array of 2
  nnid = -1*np.ones_like(pixels)  #create an array of 2

  d[0]=1  #the first galaxy goes automatically to d1

  for i in xrange(1, d.size):  #loop through the rows of the input catalogue

    #find all the pixels that intersect a disk of radius '2*diameter' around the object
    overlaps = H.query_disc(nside, vector[i,:], 2*diameter, inclusive=True, nest=nest)
    overlaps = overlaps[overlaps<=pixels.max()]   # exclude pixels with id larger than the maximum one

    #find all the indeces of the objects that are in the pixels that overlap the disc
    indeces = [y for x in overlaps for y in pointers[x]]
    indeces.remove(i)  #remove the object itself from the list

    collide = 0   #number of object in sample d1 nearer than the given diameter from the one considered
    dist = np.infty   #set the distance to infitity to search for the nearest neighbour

    for gal in indeces:  #loop though the nearby objects
      if( d[gal]==1 ):   #if it is in sample 1, go ahead
	distance = H.rotator.angdist(vector[i,:], vector[gal,:])   #find the distance
	if( distance < diameter ):  #check if is less than the required one
	  collide += 1  #if the distance is less than the required, add 1
	  if( distance < dist ):   #find the nearest neighbour
	    dist = distance   #save the distance
	    nnid[i] = gal       #and the index
    #finished looping through the nearby objects

    #print(collide)
    if( collide == 0 ): d[i] = 1  #if the object does not collide with any other one, move it to sample 1
  #finished looping through the rows of the input catalogue

  return d, nnid  #returns d and nnid

def downsample(sectid, frac_redshift, subsample=None): 
  """In each sector randomly select 'frac_redshift' objects and 
  returns 1 for them. For the other returns 0.

  Parameters
  ----------
  sectid: 1d-array
    sector id of each object. Must have the same size of 'frac_redshift'
  frac_redshift: 1d-array
    fraction of object to reject due to fiber collision for sector for each object. 
    Must have the same size of 'sectid'
  subsample: 1d-array (optional)
    array of indeces of sectid and frac_redshift where to do the downsampling.
    If given, for all the objects not in subsample 1 is returned.
    If not given the full input array is downsampled.

  output
  ------
  ds: 1d-array
    array of 1 or 0 if the object must be kept or discarded
  """
  tosave = np.ones_like(sectid)   #create an arry of ones
  if(subsample == None):   #if no subsample given
    subsample = np.arange(sectid.size)
  subsectid = sectid[subsample]   #extract the subsample
  subfrac_redshift = frac_redshift[subsample]
  #find the unique sectors id of the objects and their position in the array
  un_sectid, indeces = np.unique(subsectid, return_index=True)
  downsample_sect = subfrac_redshift[indeces]  #get fraction for the downsampling
  for sid, dss in it.izip(un_sectid, downsample_sect):
    gal_index = (subsectid==sid).nonzero()[0]    #find all the galaxies in the sector 
    gal_index = gal_index[ np.random.rand(gal_index.size)>dss ] #random selection of the galaxies to discard
    tosave[subsample[gal_index]] = 0    #set to zero objects to discard
  return tosave  #return the array of 0 and 1

if __name__ == "__main__":   #if it's the main

  args = parse(sys.argv[1:])

  args.fdiam = cf.dec2rad( cf.arcsec2deg(args.fdiam) )  #convert the fiber diameter in radians

  #read the file containing the sector to which each polygon belong and the number of tiles per sector
  poly2sect, frac_redshift = cf.read_poly2sect(args.poly2sect, args.verbose)
  if( args.fc_file == None ): #if there is no file with the completeness per sector use the default approximate values
    #fraction of the collided galaxies in CMASS which have redshifts 
    #in sectors with 1, 2 or 3 tiles. From Guo, Zehavi & Zheng 2011
    frac_d2_observed = np.array([0.05, 0.76, 0.88, 0.88])
    frac_redshift = frac_d2_observed[frac_redshift-1]    #number of observed redshifts per polygon
  else:  #read the file containing the percentage of observed redshifts per sector
    frac_redshift = np.loadtxt( args.fc_file, usecols=(int(args.fc_column),) )

  for fn in args.ifname:   #loop through the input file names
 
    if(args.verbose == True):
      print("Process catalogue '{0}'.".format(fn.name))
    cat = np.loadtxt(fn)
    if(args.verbose == True):
      print("Catalogue read")
    if(args.shuffle != None):  #if it is required to shuffle the catalogue
      np.random.shuffle(cat)  #shuffle the catalogue
      ofilesh, skip = mf.insert_src(fn.name, args.shuffle, overwrite=args.soverwrite, skip=args.sskip)   #create the file name for the shuffled catalogue
      if(skip == True):
	print("Skipping")
	continue
      fn = open(ofilesh, "w")   #open the shuffled catalogue output file name
      if(args.header == True):
	of.write(headersh)  #add the header
      np.savetxt(fn, cat, fmt=args.sfmt, delimiter="\t")   #save the file
      fn.close()

    #create the output file name
    if(args.replace == None):
      ofile, skip = mf.insert_src(fn.name, args.insert, overwrite=args.overwrite, skip=args.skip)
    else:
      ofile, skip = mf.replace_src(fn.name, args.replace, overwrite=args.overwrite, skip=args.skip)
    if(skip == True):
      print("Skipping")
      continue

    if(args.israndom == False):   #if the files are actual catalogues
      theta, phi = cf.radec2rad(cat[:,0], cat[:,1])  #convert ra and dec to radiant coordinates
      pixels = H.ang2pix(args.nside, theta, phi, nest=True)  #get the pixel number of the coordinates
      vec = H.ang2vec(theta, phi)   #convert the angles to 3D vectors: used in query_disk
      if(args.verbose == True):
	print("Healpix called")
	print("Beginning the separation")

      d, nnid = separate_d1d2(args.nside, pixels, vec, args.fdiam, nest=True)  #do the separation

      if(args.verbose == True):
	print("Sample separated.\nDownsampling sample d2 according to the observed fiber collision")
    #end if(israndom )
    
    #Downsample sample d2
    #obtain the polygon ID with mangle
    polyid = cf.read_polygons(args.polygons.name, fn.name, args.temp.name)
    sectid = poly2sect[ polyid ]   #get the sector where the galaxy in d2 lays and the number of tiles in that sector
    if( args.fc_file == None ): #if there is no file with the completeness per sector use the default approximate values
      sect_frac_d2 = frac_redshift[ polyid ]  #get the fractions of galaxies with measured redshift per sector
    else:
      sect_frac_d2 = frac_redshift[ sectid ]  #get the fractions of galaxies with measured redshift per sector

    #do the downsampling of sample d2
    if(args.israndom == False):   #if the files are actual catalogues
      zobserved = downsample(sectid, sect_frac_d2, subsample=(d==2).nonzero()[0])  #downsample only d2
    else:   #if the files are randoms
      if(args.verbose == True):
	print("Downsampling the random catalogue")
      zobserved = downsample(sectid, sect_frac_d2)  #downsample everything
      d = 2*np.ones_like(sectid)  #create an array of 2
      nnid = -1*np.ones_like(sectid)  #create an array of -1

    if(args.verbose == True):
      print("Saving the output file in '{}'".format(ofile))
    #print out the file 
    of = open(ofile, 'w')
    if(args.header == True):
      of.write(header)  #add the header
    np.savetxt(of, np.vstack( (cat.T, d, zobserved, nnid) ).T, fmt=args.fmt, delimiter="\t")
    of.close()
  
  #finished looping through the file names

  exit()

