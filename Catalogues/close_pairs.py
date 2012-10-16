#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Given a mangle mask and a (list of) catalogue(s), it gets back the 
polygon id. Then compute the near neigbour pairs and assign the redshift of one
to the other in order to simulate the fiber collision 'nearest neigbour' correction"""

import astrometry.libkd.spherematch as match
import itertools as it
#import mangle as m
import numpy as np
import optparse as op  #import optsparse: allows nice command line option handling
import os
import subprocess as sp

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

  p.set_usage(""" %prog [options] catalogue(s) 
      'catalogue(s)' are one or more catalogues each having ra and dec in the first two columns (this can be adjsuted in the options)

      Given a mangle mask and a (list of) catalogue(s), it gets back the 
      polygon id. Then compute the near neigbour pairs and assign the redshift of one
      to the other in order to simulate the fiber collision 'nearest neigbour' correction
      """)

  p.add_option("-v", "--verbose", action="store_true", dest="verb", help="Verbose")

  p.add_option("--ra-dec", dest="rd", action="store", nargs=2, type=int, default=[0,1], help="Columns containing ra and dec in 'catalogue(s)'. [Default= %default]")
  #p.add_option("--redshift", dest="z", action="store", nargs=1, type=int, default=[2], help="Columns containing the redshift in 'catalogue(s)'")

  p.add_option("-p", "--polygons", action="store", dest="polyfile", default="/data01/montefra/BOSS/Catalogs/Geometry/boss_geometry_2011_06_10.ply", help="Mangle file containing the survey geometry. [Default: %default]")
  p.add_option("-s", "--poly2sect", action="store", dest="poly2sect", default="/data01/montefra/BOSS/Catalogs/Geometry/boss_geometry_100611_secntiles.dat", help="File containing in the first column the number of the sector to which the polygon identified by the line number (starting at z) belongs  [Default: %default]")
  p.add_option("-c", "--mask", action="store", dest="mask", default="/data01/montefra/BOSS/Catalogs/Completeness/mask-CMASS-full-220911-details.txt", help="File containing the sector number (first column), North-South flag (second column, 2=North, 0=South, 1=North_empty),  the fraction of close-pairs observed out of those in target file (second to last) and the fraction of galaxies that need to be removed to remove close pairs (last column). [Default: %default]")
  p.add_option('--temp', action="store", dest="temp", default="temp_polyid.txt", help="Temporary file to store the output of 'polyid'. Make sure to change it when running this program more than once in parallel from the same directory. [Default: %default]")

  p.add_option("-f", "--fibre_collision", action="store", dest="fdiam", type=float, default=62., help="Diameter of the fibers for fiber collision in arcseconds. [Default: %default]")

  p.add_option("-o", "--outputsub", action="store", dest="outsub", nargs=2, default=["", ".out"],
               help="""Modify the input file name and use it as output file. 
	       If outsub[0]="", output[1] added to the end of the input file name,
	       if outsub[1]="", output[0] added at the beginning of the input file name,
	       else substitute outsub[0] with output[1]. [Default: %default]""")

  p.add_option("--fmt", action="store", dest="fmt", default="%+10.9e", help="Format of the output file. [Default: %default]")
  p.add_option("--delimiter", action="store", dest="delim", default="\t", help="Delimiter in the output file. [Default: '%default']")

  return p.parse_args()

if __name__ == "__main__":   # if is the main
 
  opt, args = options(op.OptionParser(version="%prog version 0.1"))   #create the object optparse

  #assign the file name
  catalogues = args

  #initialise mangle
  #mng = m.Mangle(opt.polyfile)
  #if(opt.verb == True):
  #  print("Mangle initialised")
  #read the polygon to sector file
  poly2sect = np.loadtxt(opt.poly2sect, usecols=[0,])
  #read the mask file
  maskdet = np.loadtxt(opt.mask, usecols=(0,1,-2,-1))
  if(opt.verb == True):
    print("Sector and completeness files read")

  #convert the fiber dimeter in degrees and add the column of containing the redshift in the catalogues to the ones of ra and dec
  opt.fdiam /= 3600.
  #opt.rd.extend(opt.z)
  test_loss = []  #used to check how much the sum of the weights differs from the number of input particles

  for fn in catalogues:    #loop through the catalogues
    if(opt.verb == True):
      print("\nFile: {0}".format(fn))
    cat = np.loadtxt(fn)  #read the catalogue
    ra, dec = cat[:,opt.rd].T  #extract ra, dec and z
    w = np.ones_like(ra)   #create a column for the weights

    if(opt.verb == True):
      print("Obtaining the sector id")

    cmd = ['polyid', '-q', opt.polyfile, fn, opt.temp]
    exe = sp.Popen(cmd, stdout=sp.PIPE).wait()
    #polyid = mng.get_polyids(ra, dec)   #get the sector where the galaxy lays: first find the polygon with mangle and then the sector
    polyid = np.loadtxt(cmd[-1],usecols=(2,), skiprows=1, dtype=int)
    sectid = poly2sect[ polyid ]   #get the sector where the galaxy lays: first find the polygon with mangle and then the sector
    unique_id = np.unique(sectid)   #find all the unique sector ids (ordered)
    if(opt.verb == True):
      print("Done")
    if(unique_id[0] == -1):
      os.exit("There is something strange with the poligons and sectors in file {0}".format(fn))

    #loop over all the sectors, find which galaxies belong to each, what is the fraction of 
    #galaxies in close pairs lost, and assign the nearest neighbour redshift
    if(opt.verb == True):
      print("Finding close pairs and put 'nearest-neigbour'-like weights according to the observed data")
    for uid in unique_id:   
      cp_obs = maskdet[maskdet[:,0]==uid,-2]   #find the fraction of close pairs observed out of those in target
      if(cp_obs == 1):   #cp_obs==1 all of the galaxies in close pairs observed and no redshift need to be substitute
	continue
      else:
	gal_index = (sectid==uid).nonzero()    #find all the galaxies in the sector 
	ras, decs = ra[gal_index], dec[gal_index]   #ra and dec of the galaxies in the sector

	#find pairs separated by less than 62". The position in ras and decs of the pairs are in m1,m2. d12 is the distance
	m1, m2, d12 = match.match_radec(ras, decs, ras, decs, opt.fdiam, notself=True)

	colliding = m1.size/2.  #number of colliding pairs
	if(colliding != 0):
	  #discard the doubles. 'match_radec' returns pairs such that m1[i]=m2[j] and m1[j]=m2[i] for i!=j
	  discard = m1<m2

	  #throw random number for every close pairs and, if it is > cp_obs, assign to one the redshift of the other
	  loss_pairs = (np.random.rand(colliding) > cp_obs)
	  which_galaxy = np.random.randint(2, size=loss_pairs.size)  #decide which of the two galaxies get the redshift of the other

	  m12 = np.vstack( (m1[discard], m2[discard]) ).T  #get all the unique pairs and put them together
	  for i, lp, wg in it.izip(it.count(), loss_pairs, which_galaxy):   #increas the weight in galaxies "observed" and decrease in the others
	    if(lp == True):
	      w[m12[i,wg]] += 1
	      w[m12[i,(wg+1)%2]] = 0   #simply put to 0 all the lost pairs (but might make disappear some galaxy that might be there) 
	      #w[m12[i,(wg+1)%2]] -= 1   #decrease by 1, but some catalogue can be discarded twice. Do this (*) to take the negative w

    #copy the revised redshifts back into the catalogue
    #w[w<0] = 0   #(*) force all the weights to be >=0

    #create the output file name avoiding to overvrite it
    if(opt.outsub[0] == ""):
      ofile = fn+opt.outsub[1]
    elif(opt.outsub[1] == ""):
      ofile = opt.outsub[0]+fn
    else:
      ofile = fn.replace(opt.outsub[0], opt.outsub[1])
    if(ofile == fn):   #if the substitution didn't work add outsub[1] at the end of the file name
      ofile = fn+".out"

    fmt = [opt.fmt,]*cat.shape[1]
    fmt.append('%d')
    np.savetxt(ofile, np.vstack( [cat.T,w.T] ).T, fmt=opt.fmt, delimiter=opt.delim)

    test_loss.append(w.sum()/float(w.size))

  print( 'sum(w)/n_object: mean={0:4.3f}, std={1:4.3f}. Should be 1'.format(np.mean(test_loss), np.std(test_loss)) )
  sp.Popen(['rm', opt.temp]).wait()

  exit()

