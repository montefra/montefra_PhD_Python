#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Read a file with column with 1 or 0 
if the object is observed or not and compute the completeness per sector"""

import argparse as ap
import common_functions as cf
import itertools as it
import my_functions as mf
import numpy as np
import sys  #system stuff

header = "#Fraction of galaxies with redshift per sector (indicated by the line number)\n"

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
  
  description = """Read a file with a column with 1 or 0 if the object is observed or 
    not and compute the completeness per sector"""
  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("ofile", action="store", help="""Output file name, if only one 
    'ifname' give or if the mean is required, output file root otherwise.""")

  p.add_argument("ifnames", action="store", nargs='+', type=ap.FileType('r'), 
    help="""Input file name(s), containing ra, dec and redshift in the first three columns 
    and with a column of 1 or 0""")

  p.add_argument('--version', action='version', version='0.1')
  p.add_argument("-v", "--verbose", action="store_true")

  p.add_argument("--outend", action="store", default="", help="String to append to the output file name(s)")

  p.add_argument("-c", "--column", action="store", default=-2, help="Column containing the 1 and 0")

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

  p.add_argument("--mean", action="store_true", default=False, help="""If more than one input file,
      compute the mean of the completeness and save only it""")

  p.add_argument("--header", action="store_true", help="""If selected the the following header will be 
      written '{0}'""".format(header) )
  p.add_argument("--fmt", default="%7.6e", action="store", nargs='+', help="Format of the output files")

  group = p.add_mutually_exclusive_group()
  group.add_argument("--overwrite", action="store_true", default=False,
      help="If given does not check if the output file exists or not before saving it.")
  group.add_argument("--skip", action="store_true", default=False,
      help="Skip already existing output file names. No operation done on the corresponding input file.")

  return p.parse_args(args=argv)


if __name__ == "__main__":   #if it's the main

  args = parse(sys.argv[1:])

  #read the file containing the sector to which each polygon belong and the number of tiles per sector
  poly2sect, n_tiles = cf.read_poly2sect(args.poly2sect, args.verbose)
  #create an array of zeros with a number of elements equal to the maximum sector id number 
  completeness = []  #list of arrays

  for fn in args.ifnames:  #loop through the input file names
    if(args.verbose == True):
      print("Read catalogue '{0}'.".format(fn.name))
    hasz = np.loadtxt(fn, usecols=[args.column]) #read the columns containing the 0 and 1

    #obtain the polygon ID with mangle
    polyid = cf.read_polygons(args.polygons.name, fn.name, args.temp.name)
    sectid = poly2sect[ polyid ]   #get the sector where the galaxies lie 

    un_sectid = np.unique(sectid)  #get the unique sector id

    #create an array of zeros with a number of elements equal to the maximum sector id number 
    compl = np.zeros(poly2sect.max()+1)  

    for sid in un_sectid:  #loop through the unique sector ids
      gal_index = (sectid==sid).nonzero()[0]    #find all the galaxies in the sector 
      s_hasz = hasz[ gal_index ]   #hasz in the sector
      #number of objects with redshifts over the total number
      compl[ sid ] = (s_hasz==1).sum() / float(s_hasz.size)
    #end loop through the unique sector ids

    completeness.append(compl)  #append the completeness for each file
  #end loop through the input file names

  if(args.mean == True):
    if(args.verbose == True):
      print("Compute the mean among the files")
    completeness = [np.mean(completeness, axis=0), ]

  for i, compl in enumerate(completeness):  #loop over the completeness list
    if( len(completeness) == 1):  #if only one file to save
      outfile = args.ofile+args.outend   #output file name
    else:  #otherwise
      outfile = "{0}.{1:03d}.{2}".format(args.ofile, i+1, args.outend)   #output file name

    toskip = mf.check_file_exist(outfile)   #check if the output file exists
    if( toskip == True and args.skip == True ):  #if the file is to be skipped go to the next step of the loop
      print( "Skipping not to overwrite file {}.".format(outfile) )
      continue  
    if( toskip == True and args.overwrite == False): #if overwriting is not required, exit the program
      print( "The file {} already exists.".format(outfile) )
      sys.exit(10)

    if(args.verbose == True):
      print("Saving the completeness to file {}".format(outfile))

    of = open(outfile, 'w')  #create the output file name
    if(args.header == True):
      of.write(header)  #add the header
    np.savetxt(of, compl, fmt=args.fmt, delimiter="\t")
    of.close()
  #end loop over the completeness list

  exit()

