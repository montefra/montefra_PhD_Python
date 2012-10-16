# -*- coding: utf-8 -*-
"""It contains functions to prepare the contour plots from the 
MCMC chains produced by COSMOMC"""

import numpy as np  #numpy
import smooth as sm  #import the gaussian kernel density estimate function 
import sys

def get_paramnames(file_roots, ext=".paramnames", verbose=False):
  """Reads the paramnames files of the cosmomc
  
  Parameters
  ----------
  file_roots: list
    list of file root of the chain to be plotted
  ext: string (optional)
    extention of the parameters file
  verbose: bool (optional)
    if True print more output
  output: list
    list of strings containing the parameters names in human readeble and latex syntax
  """
  for i in file_roots:    # goes through the paramnames files to get the list of parameters
    parafile = i+ext
    if(verbose == True):
      print "Reading file %s" % parafile
    parf = open(parafile, 'r')
    try:   # if the list of string from paramnames already exist concatenate otherwise create
      paramnames
    except NameError:
      paramnames = parf.readlines()
    else:
      if(paramnames != parf.readlines()):
        print "The parameter name file '%s' do not correspond with the previous one" %parafile
	sys.exit(10)
    parf.close()
  return paramnames   #return the list o parameters names


def get_chains(file_roots, cols, ext="_total.txt", verbose=False):
  """Reads the paramnames files of the cosmomc

  Parameters
  ----------
  file_roots: list
    list of file root of the chain to be plotted
  cols: list
    list containing the columns to read from the files
  ext: string (optional)
    extention of the chain files
  verbose: bool (optional)
    if True print more output
  output: list of array
    list of array containing all the input chains
  """

  chains = []   #save the arrays of the chains
  for i in file_roots:     #loop through the total chains
    chainfile = i+ext
    if(verbose == True):
      print "Reading file %s" % chainfile
    chains.append( np.loadtxt(chainfile, usecols=cols) )   #load the chain file

  return chains


def hist2D(a, xr=None, yr=None, bins=30, smooth=False):
  """
  Given a list of arrays a containing len(a) indipendent chains concatenated along axis=0
  and whose extrama are stored in 'row' returns a [bins*(bins*(len(row)-1))] array 
  containing all the 2D histograms
  Smoothing is applied is required

  Parameters
  ----------
  a: list of arrays
    each element must have 2 or 3 columns. If has 3 columns the first one will be considered 
    as containing the weights.
  xr: list (optional)
    x_range of the histogram: if not given will be determined from 'a' itself
  yr: list (optional)
    y_range of the histogram
  bins: int (optional)
    number of bins for the 2D histogram. This numer will be doubled when smothing
  smooth: bool (optional)
    if False numpy.histogram2D used, otherwise fast_kde (gaussian kernel density estimate)

  output: list of arrays, list of arrays, list of arrays
    hist2D: list of arrays containing all the 2D histograms
    xe: list of arrays containing the bin means along the first dimension
    ye: list of arrays containing the bin means along the second dimension
  """

  h2D, xe, ye = [], [], []  #lists of outputfiles
  for c in a:  #loop over the input list of arrays
    if(c.shape[1] == 2):  #if a has 2 columns add 1 in the first one
      c = insert(c, 0, 1., axis=1)
    elif(c.shape[1] == 3):  #if a has 3 columns do nothing
      pass
    else:
      print "Too many or too few columns"
      sys.exit(20)

    txr, tyr = xr, yr
    if(xr == None):  #if not given get the x and y ranges
      txr = [np.amin(c[:,1]), np.amax(c[:,1])]
    if(yr == None):
      tyr = [np.amin(c[:,2]), np.amax(c[:,2])]

    if smooth==True:
      h2D.append(np.flipud(sm.fast_kde(c[:,1], c[:,2], gridsize=(bins,bins), extents=(txr[0],txr[1],tyr[0],tyr[1]), weights=c[:,0])))  #create a smooth 2D histogram
      txe, xst = np.linspace(txr[0], txr[1], num=bins, endpoint=False, retstep=True)
      tye, yst = np.linspace(tyr[0], tyr[1], num=bins, endpoint=False, retstep=True)
      txe, tye = txe+0.5*xst, tye+0.5*yst
      xe.append(txe)  #arrays containing the mean of x and y axis of the 2d histogram
      ye.append(tye)  #arrays containing the mean of x and y axis of the 2d histogram
    else:
      th2D, txe, tye = np.histogram2d(c[:,1], c[:,2], bins=bins, range=[txr, tyr], weights=c[:,0])  #create the 2D histogram
      h2D.append(th2D.T)  #append the histrogram
      txe, tye = (txe[:-1] + txe[1:])/2., (tye[:-1] + tye[1:])/2.
      xe.append(txe)  #list of arrays containing the mean of x and y axis of the 2d histogra
      ye.append(tye)

  return h2D, xe, ye


def h2D2conflev(a, levels, ind=None):
  """
  This function gets 2D histograms and a list of confidence levels
  and returns an 2D array containing the amplitudes (relative to the maxima)
  of the histogram that correspond to the given levels for the histograms in a

  a: array like
    array containing squared 2D histograms. The codes do not check wherever those are squared or not
  levels: 1D array
    1D array containing the required confidence levels
  ind: list (optional)
    The n elements of 'ind' are considered as the index extrema for the n-1 indipendent
    histogramns stored in 'a' along axis=0

  output: 2D array like
    amplitude, relative to the maximum, of the 2D histograms corresponding to the confidence
    levels 'levels'
  """

  levels.sort()   #sort so there's no problem
  n_levels = levels.size   #number of levels

  output = []  #create a list containing all the 
  for c in a:
    [biny, binx] = c.shape
    toutput = np.zeros(n_levels+1)   #allocate the temporary output
    sc = np.sum(c)  #sum of all the elements into the 2d histogram  
    toutput[0] = np.amax(c)
    tlev = sc*levels

    sortc = np.argsort(c.flatten(), kind='mergesort')[::-1]   #flatten the histogram. sort its arguments and invert the order

    tot = 0.  #initialize the sum
    l = 0  #integer to go through the levels
    for ind in sortc:  #go through the indeces
      tc=c[ind/binx,ind%binx]
      tot += tc    #compute the sum of the amplitude starting for decreasing amplitude of the 2D histogram
      if(tot == tlev[l]):
        toutput[l+1] = tc
        l += 1
        if(l >= n_levels):
          break
      elif(tot > tlev[l]):
        toutput[l+1] = ((tc-tco)*tlev[l] + tco*tot - tc*(tot-tc)) / tc
        l += 1
        if(l >= n_levels):
          break
      tco=tc  #save the value of amplitude of the histogram for the next loop if interpolation needed

    output.append(toutput[::-1])   #append the levels of the 2D histogram

  return output

