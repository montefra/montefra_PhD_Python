#!/usr/bin/python
# -*- coding: utf-8 -*-
#merge window functions in a unique file 

import glob
import numpy as np
import optparse as op
import smooth

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

  p.set_usage("""
  %prog [options] output_filename input_filenames
  Given a list of input file names (also shell-like wildcard) containing the window functions,
  it unify them and save them in 'output_filename'.
  """)

  p.add_option("-f", "--fraction-kN", action="store", type=float, dest="fkN", default=0.65, help="All the modes larger than fkN*k[-1] will be discarded. [Default: %default]")

  p.add_option("-k", "--k-average", action="store", type=float, dest="kav", help="If given, average the window function for k>'kav' and substitute this value to the true window function. at these scales noise dominates")

  p.add_option("-s", "--smooth", action="store", dest="smooth", nargs=2, help="Smooth the signal for k>'kav' with the window 'smooth[0]' and with a width of 'smooth[1]'. Choices for the window are:'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

  return p.parse_args()

def unify_win(wins, noise=0, fkN=0.65):
  """Unify the window functions in a unique object.
  The windows have to be ordered from larger to smaller scale coverage.
  For each window function all the modes larger than fkN*wins[i][0,-1] 
  (fkN of the Nyquist frequency) will be discarder. When two window functions 
  overlap the one on larger scale will be considered.

  Parameters
  ----------
  wins: list of array
    list of window functions. Each element with at least 2 columns: k, W(k)
  noise: float (optional)
    shot noise to subtract from W(k)
  fkN: float (optional)
    fraction of the maximum wave number to discard

  output: nD array
    array with the same number of columns of the input ones, containing the merged window function

  The input wins needs to be ordered from the more large scales to the more small scales """
  k_min = 0.
  win = []   #output total window
  for w in wins:
    w[:,1] -= noise   #subtract the shot noise
    k = w[:,0]    #get the k
    imin, imax = np.sum(k<k_min), np.sum(k<(fkN*k[-1]))   #get the first and the last bins
    win.extend(w[imin:imax,:])    #add the part of interest of the window function 
    k_min = k[imax]  #substitute the k_min for the next loop with the largest k of this loop

  return np.array(win)  #return the full window

if __name__ == "__main__":   # if is the main

  opt, args = options(op.OptionParser(version="%prog version 1"))   #create the object optparse

  wins = []
  for a in args[1:]:   #go through the various input file names
    for fn in glob.iglob(a):
      wins.append(np.loadtxt(fn))
  if(len(wins)==0):
    raise IOError("No valid input file name given among {0}".format(", ".join(args[1:])))

  order = np.empty(len(wins))  #check the order of the window functions through wins[i][0,0]
  for i,w in enumerate(wins):
    order[i] = w[0,0]  #copy the first k of every window 
  wins = [wins[i] for i in np.argsort(order)]   #list of window function ordered from smaller to larger wins[i][0,0]

  win = unify_win(wins, fkN=opt.fkN)   #unify window functions

  if(opt.kav != None):  #substitute every value for k>kav with the mean of hte window for those k
    kindex = (win[:,0]>opt.kav).sum()
    if(opt.smooth == None):
      win[-kindex:,1] = np.mean(win[-kindex:,1])
    else:
      print("understand how to implement it properly")
      win[-kindex:,1] = smooth.smooth1D(win[-kindex:,1], window=opt.smooth[0], window_len=int(opt.smooth[1]))[:kindex]

  np.savetxt(args[0], win, fmt='%7.6e', delimiter='\t')

  exit()

