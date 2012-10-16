#!/usr/bin/python
# -*- coding: utf-8 -*-
#given a list of files read the column with the redshift and computes
# and plots or save the n(z) in [h/Mpc]^2

import argparse_custom as apc
import matplotlib.pyplot as plt
import numpy as np
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
  import argparse as ap

  description = """ %prog [options] area filenames
  Given a list of file which contain the redshift in one column [the third by default]
  it computes the value of N(z) or n(z) if the area in deg^2 given .
  The N(z) or n(z) computed for each file are plotted and showed on the screen.
  Through the options it is possible to save them file, as well as decide if plot or save 
  the average or the single values. The plot are done for fast checks."""
  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("catname", action="store", nargs="+", type=ap.FileType('r'), 
      help="File name of the catalogue containing the redshfit.")
  
  p = apc.version_verbose(p, '0.1')

  p.add_argument("-c", "--column", action="store", type=int, default="2", 
      help="Column in the input file names which contains the redshift.")
  p.add_argument("-w", "--weight", action="store", type=int, 
      help="If used, column in the input file names which contain the weight.")

#  p.add_argument("--skip", action="store", type=int, default=0, 
#      help="Skip 'skip' rows at the beginning of the file.")

  p.add_argument("-n", "--nbins", action="store", type=int, default='50', 
      help="Number of bins per histogram.")
  p.add_argument("-r", "--range", action="store", nargs=2, type=float, 
      help="""Upper and lower range of the bins. If not given it is fixed to 
      the extremes of the first file given""")

  p.add_argument("--mean", action="store_true", help="""If more than one file is given, 
      the mean is computed before plotting or saving the n(z).""")

  description="""If the area of the survey is given, the effective 
      volume is computed and the histogram of n(z) is returned. The cosmological parameters can be adjusted """
  p, cosmo = apc.cosmology_group(p, description=description )
  cosmo.add_argument("-a", "--area", type=float, action="store", help="Area of the survey in steradians")

  p.add_argument("--legend", action="store_true", default=False, 
      help="Legend with the input file name drawn")

  p.add_argument("-s", "--save", nargs="+", action=apc.required_length(1,2),
      help="""If given, the N(z) or n(z) are saved in oname[0](+oname[1]), 
      if only one input file is given or the mean is required, or oname[0]+counter(+oname[1]), otherwise. 
      This function can be called with one or two arguments.""")

  p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+', help="Format of the output files")

  return p.parse_args(args=argv)

def nofz(z_edges, hists, area, om, ok, wde, h0=None):
  """
  Compute the effective volume in (Mpc/h)^3 or Mpc^3 for each redshift bin 
  and return hists[i]/effective_area

  Parameters
  ----------
  z_edges: 1d array
    edges of the the istogram intervals: used to obtain the effective volume in the interval
  hists: list
    list of 1d arrays containing the histograms. The must all have the same edges
  area: float
    area of the footprint in deg^2
  om: float
    Omega_Matter
  ok: float
    Omega_curvature
  wde: float
    dark energy equation of state parameter
  h0: float (optional)
    reduced hubble parameter. If given the effective volume computed in unit of Mpc^3

  output
  ------
  nofz: same format as 'hists'
    list of 1d array containing the input histograms divided by the effective volume per redshif bin
  """

  import cosmologydir as c
  import common_functions as cf

  cos = cf.set_cosmology(om, ok, wde, h0=h0)
  #set the distance object
  dis = c.distance.Distance(cos)

  #compute the effective volumes
  if(h0 == None):
    eff_vol = dis.effective_volume_sr_zh(z_edges) * area
  else:
    eff_vol = dis.effective_volume_sr_z(z_edges) * area

  nofzh = []   #initialize the output array
  for h in hists:  #compute the number density in [Mpc/h]^3 or Mpc^3
    nofzh.append(h/eff_vol)

  return nofzh

def plot_noz(ax, z, hists, std=None, leg=None):
  """
  Plot on axes 'ax' the histograms

  Parameters
  ----------
  ax: matplotlib Axes or AxesSubplot object
    axes where to plot the histograms
  z: 1d array
    absissa of the plot
  hists: list of 1d arrays
    histograms to be plotted. Same dimension as z
  std: 1d arrays (optional)
    standard deviation. If given 'ax.errorbars' is used instead of 'plot'. Same dimension as z
  leg: list of strings (optional)
    strings for the legend

  output
  ------
  ax: matplotlib Axes or AxesSubplot object
    return the axes with the linea objects added
  """
  if( leg == None):  #if the legend is not given create a list of empty strings
    leg = ["",]*len(hists)
    isleg = False  #no legend
  else:
    isleg = True  #legend requred
  if(std == None):   #if no std to be plotted
    for h,l in zip(hists, leg):
      ax.plot(z, h, label=l)
  else:
    for h in hists:
      ax.errorbar(z, h, std)

  if( isleg ): #if the legend is required
    ax.legend(loc=0)   #draw it
  return ax


if __name__ == "__main__":   #if it's the main

  args = parse(sys.argv[1:])

  hists=[]  #initialise the list containing all the histograms

  # read the first file and fix range if not given
  if(args.verbose):
    print("Processing file:", args.catname[0].name)
  if( args.weight == None ):
    z, w = np.loadtxt(args.catname[0], usecols=[args.column,], ), None #skiprows=args.skip), None
  else:
    z, w = np.loadtxt(args.catname[0], usecols=[args.column, args.weight], ).T # skiprows=args.skip).T
  hist, z_edges = np.histogram(z, bins=args.nbins, range=args.range, weights=w)   #create the first histogram
  hists.append(hist)  #append to the list
  if(args.range == None):  #if the range is not given, save it
    args.range = [z_edges[0], z_edges[-1]]

  for fn in args.catname[1:]:   #loop through the remaining file names
    if(args.verbose):
      print("Processing file:", fn.name)
    z = np.loadtxt(fn, usecols=[args.column,], skiprows=args.skip)
    hist, bin_edges = np.histogram(z, bins=args.nbins, range=args.range)   #create the histograms
    hists.append(hist)  #append to the list

  if(args.mean == True):   #if the mean is computed
    sigmahist = np.std(hists, axis=0)
    hists = [np.mean(hists, axis=0), ]

  if(args.area != None):   #if the area is give compute n(z)
    if(args.verbose):
      print("Computing the effective volume for the redshift bins")

    if(args.mean == True):   #if the mean is computed
      hists.append(sigmahist)   #pass mean and std together to nofz
      hists = nofz(z_edges, hists, args.area, args.om, args.ok, args.wde, args.h0)
      sigmahist = hists.pop(1)   #separate back mean and std
    else:    #if no mean requred loop through the histograms
      hists = nofz(z_edges, hists, args.area, args.om, args.ok, args.wde, args.h0)

  #save or plot the histogram(s)
  zbins = (z_edges[1:]+z_edges[:-1])/2.   #get the mean of the histogram bins
  if(args.save == None):    #plot stuff
    if(args.verbose):
      print("Plot N(z) or n(z).")

    ax = plt.figure().add_subplot(111)   #Axes object
    if(args.mean == True):   #if the mean is computed
      ax = plot_noz(ax, zbins, hists, std=sigmahist)
    else:
      if( args.legend == True ):
	leg = [fn.name for fn in args.catname ]  #legend of file names
      else:
	leg = None
      ax = plot_noz(ax, zbins, hists, leg=leg)

    ax.set_xlabel('z')   #set the axes label 
    if(args.area == None):  #N(z)
      ax.set_ylabel('N(z)')
    else:
      if(args.h0 == None): #n(z)
	ax.set_ylabel('n(z) [(Mpc/h)$^3$]')
      else:
	ax.set_ylabel('n(z) [Mpc$^3$]')
    plt.tight_layout()
    plt.show()
  #end plot if
  else:   #print the n(z)
    if(args.verbose):
      print("Save N(z) or n(z) in files.")
    if(len(args.save) == 1):
      args.save.append("")   #append an empty string if only one argument is given for the output 

    if(args.mean == True ):   #save single file
      np.savetxt( ''.join(args.save), np.vstack([zbins, hists[0], sigmahist]).T, delimiter="\t", fmt=args.fmt )
    elif( len(args.catname)==1 ):
      np.savetxt( ''.join(args.save), np.vstack([zbins, hists[0]]).T, delimiter="\t", fmt=args.fmt )
    else:    #save multiple files with a counter
      for i, h in enumerate(hists):
	np.savetxt( "{0:03d}".format(i).join(args.save), np.vstack([zbins, h]).T, delimiter="\t", fmt=args.fmt )

  #end print if

  exit()

