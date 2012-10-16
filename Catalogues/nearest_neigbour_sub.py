#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Read a catalogue with ra, dec, z and substitute to the collided galaxies the position of the nearest neighbour
Now is assumed that the last column is 1 if the object has redhisft and 0 if it doesn't
"""

import argparse as ap
import common_functions as cf
import healpy as H
import my_functions as mf
import numpy as np
import sys  #system stuff

header = "#ra	dec	z	\n"

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
  
  description = """Assign to the galaxies which have no redshift, the redshift or also the position of the closest neigbour.
  The second to last column of the input catalogue contains 1 if the redshift has to be kept, 0 if it has to be 
  substituted with the nearest neigbour one, whose index is in the last column. 
  The output contains only the first three columns (ra, dec and z)"""
  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'), 
      help="""Input file name(s), containing ra, dec and z in the first three columns 
      and 1 or 0 in the second to last ones. -1 or index of the nearest neighbour in the last column.""")

  p.add_argument('--version', action='version', version='0.1')
  p.add_argument("-v", "--verbose", action="store_true")

  ir_group = p.add_mutually_exclusive_group()
  ir_group.add_argument("-i", "--insert", action="store", nargs=2, default=[".d1d2", ".dat"], 
      help="""Output file name created inserting '%(dest)s[0]' before '%(dest)s[1]'.
      Only the first three columns of the file are saved""")
  ir_group.add_argument("-r", "--replace", action="store", nargs=2, 
      help="""Output file name created replacing '%(dest)s[0]' with '%(dest)s[1]'.
      Only the first three columns of the file are saved. If none of these two options given, 'insert' is assumed""")

  group = p.add_mutually_exclusive_group()
  group.add_argument("--overwrite", action="store_true", 
      help="If given does not check if the output file exists or not before saving it.")
  group.add_argument("--skip", action="store_true", 
      help="Skip already existing output file names. No operation done on the corresponding input file.")

  p.add_argument("-a", "--radec", action="store_true", default=False,
      help="Assign to the galaxies with no redshifts also ra and dec of the nearest one")

  p.add_argument("--header", action="store_true", help="""If selected the the following header will be 
      written '{0}'""".format(header) )
  p.add_argument("--fmt", default="%7.6e", action="store", nargs='+', help="Format of the output files")

  return p.parse_args(args=argv)



if __name__ == "__main__":   #if it's the main

  args = parse(sys.argv[1:])

  for fn in args.ifname:   #loop through the input file names
    
    #create the output file name
    if(args.replace == None):
      ofile, skip = mf.insert_src(fn.name, args.insert, overwrite=args.overwrite, skip=args.skip)
    else:
      ofile, skip = mf.replace_src(fn.name, args.replace, overwrite=args.overwrite, skip=args.skip)
    if(skip == True):
      print("Skipping")
      continue

    if(args.verbose == True):
      print("Process catalogue '{0}'.".format(fn.name))
    cat = np.loadtxt(fn)
    if(args.verbose == True):
      print("Catalogue read")

    indnoz = (cat[:,-2]==0).nonzero()[0]   #indeces of the objects without redshift
    indnn = np.array( cat[indnoz, -1], dtype=int )    #indeces of the nearest neighbours

    if(args.verbose == True):
      print("Subsituting redshifts (and positions)")
    if(args.radec == False):  #substitue the redshift of the object with the nearest neighbour one
      cat[indnoz,2] = cat[indnn,2] 
    else:   #subsitute ra, dec and z with the nearest neighbour ones
      cat[indnoz,:3] = cat[indnn,:3] 

    if(args.verbose == True):
      print("Saving the output file in '{}'".format(ofile))
    #print out the file 
    of = open(ofile, 'w')
    if(args.header == True):
      of.write(header)  #add the header
    np.savetxt(of, cat[:,:3], fmt=args.fmt, delimiter="\t")
    of.close()

  #input file names loop finished

  exit()
