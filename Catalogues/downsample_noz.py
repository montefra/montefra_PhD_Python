#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Given a file with the fraction for downsample as function of redshift and one or more
catalogues, downsample them according to the given fractions"""

import argparse as ap
import argparse_custom as apc
import itertools as it
import my_functions as mf
import numpy as np
import scipy.interpolate as spip
import sys  #system stuff

header = "#ra	dec	z\n"

def strORint(string):
  """Check if string is an integer or one of these strings
  'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic' and return it with the correct type.
  Otherwise through an error.
  
  Parameters
  ----------
  string: string
    input string to check
    
  output
  ------
  value: int or string
  """

  stchoise = ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'] 
  try:
    value = int(string)  #check if it is an int
  except ValueError:   #if not check if it's a string among the possible choises
    if( string in stchoise ):
      value = string
    else:
      msg = "'{}' is neither an integer nor a string among those: {}".format(string, stchoise)
      raise ap.ArgumentTypeError(msg)
  return( value )

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
  
  description = """Given a file with the fraction for downsample as function of redshift and 
      one or more catalogues, downsample them according to the given fractions"""
  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("downsample", action="store", type=ap.FileType('r'), 
      help="""File with the fraction to downsample as function of redshift.""")
  p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'), 
      help="Input file name(s), containing z in one of the columns")

  p.add_argument('--version', action='version', version='0.1')
  p.add_argument("-v", "--verbose", action="store_true")

  ir_group = p.add_mutually_exclusive_group()
  ir_group.add_argument("-i", "--insert", action="store", nargs=2, default=[".downsampled", ".dat"], 
      help="""Output file name created inserting '%(dest)s[0]' before '%(dest)s[1]' in the input file name.""")
  ir_group.add_argument("-r", "--replace", action="store", nargs=2, 
      help="""Output file name created replacing '%(dest)s[0]' with '%(dest)s[1]' in the input file name.
      If none of these two options given, 'insert' is assumed""")

  group = p.add_mutually_exclusive_group()
  group.add_argument("--overwrite", action="store_true", 
      help="If given does not check if the output file exists or not before saving it.")
  group.add_argument("--skip", action="store_true", 
      help="Skip already existing output file names. No operation done on the corresponding input file.")

  p.add_argument("-z", "--z-column", action="store", type=int, default=2,
      help="Column containing the redshift")

  p.add_argument("-k", "--kind", action="store", type=strORint, default="linear",
      help="""Specifies the kind of interpolation of the downsample file as 
      a string (‘linear’,’nearest’, ‘zero’, ‘slinear’, ‘quadratic, ‘cubic’) 
      or as an integer specifying the order of the spline interpolator to use. Default is ‘linear’.""")

  p.add_argument("--header", action="store_true", help="""If selected the the following header will be 
      written '{0}'""".format(header) )
  p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+', help="Format of the output files")

  return p.parse_args(args=argv)

if __name__ == "__main__":   #if it's the main

  args = parse(sys.argv[1:])

  #read the fraction of according to which the downsample has to be done
  z, fdown = np.loadtxt(args.downsample).T  

  #function with the interpolation. 0 is returned outside x and y range
  func = spip.interp1d(z, fdown, kind=args.kind, fill_value=0., bounds_error=False)

  for fn in args.ifname:  #file name loop
    if(args.verbose == True):
      print("Process catalogue '{0}'.".format(fn.name))

    #create the output file name
    if(args.replace == None):
      ofile, skip = mf.insert_src(fn.name, args.insert, overwrite=args.overwrite, skip=args.skip)
    else:
      ofile, skip = mf.replace_src(fn.name, args.replace, overwrite=args.overwrite, skip=args.skip)
    if(skip == True):
      print("Skipping")
      continue

    cat = np.loadtxt(fn)  #read the catalogue
    if(args.verbose == True):
      print("File read. Downsampling")

    z = cat[:, args.z_column]  #get the redshifts
    keep = np.random.rand(z.size) <= func(z)  #random select objects to keep
    #end zbins loop
    if(args.verbose == True):
      print("Saving the downsampled catalogue into file {}".format(ofile))
    cat = cat[ keep, : ]  #downsample the catalogue
    of = open(ofile, 'w')
    if(args.header == True):
      of.write(header)  #add the header
    np.savetxt(of, cat, fmt=args.fmt, delimiter="\t")   #save the file
    of.close()

  #end file name loop

  exit()

    


