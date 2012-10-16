#!/usr/bin/python
# -*- coding: utf-8 -*-
#returns the parameters of a model obtain with a chi^2 minimisation
#inspired by the package "leastsq.py" 
#from 'http://sites.google.com/site/applegatearchive/software/python-model-fitting'

import leastsq as lchi
import numpy as np
import optparse as op
import ps_func as pf
import scipy.optimize as spo


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
  %prog [options] measured kmin kmax
  It computes the best model parameters given a measured power spectrum
  and a model (defined in function 'model(x, params)') using a chi^2 minimisation for
  kmin < k < kmax.
  'measured' must contain at least two columns x, f(x). If it has third column,
  and no variance or covariance file is given (option: not set for now), this is
  assumed to be the variance.
  The convolution with a window function before performing the chi^2 is done if
  Wij and kj are given.
  """)

  p.add_option("-l", "--linear", action="store", dest="flin", default="/data01/montefra/SDSSDR7/Theory/lasdamas_matterpower_z0.278.dat", help="File name of the linear power spectrum. [Default: %default]")
  p.add_option("-m", "--mode-coupling", action="store", dest="f1loop", default="/data01/montefra/SDSSDR7/Theory/lasdamas_1loop_log_z0.278.dat", help="File name of the linear power spectrum. [Default: %default]")
  p.add_option("-g", "--guess", action="store", dest="guess", type=float, default=[2., 0.2, 0.7], help="Initial guess of the model parameters: b, k_star, A_MC. [Default: %default]")

  p.add_option("-v", "--variance", dest='isvar', action='store_true', default=False, help='If set the errors used. If "cov" and "covcol" are not set, the variance is assumed to be in the third column of "measured".')
  p.add_option('-c', '--covariance', dest='cov', action='store', help='File name containing the covarince matrix. The variance used is the square root of the diagonal elements.')
  p.add_option('--cov-column', dest='covcol', action='store', type=int, help="If a number given, the covariance matrix is assumed to be in the 'covcol' column of the file: only that column is reshaped to a square.")

  p.add_option("-w", "--window", action="store", dest="win", nargs=2, help="If given, the two strings must be the window matrix 'W_ij' and the values of k_j.")
  p.add_option("--window-extra", action="store", dest="win_extr", nargs=2, help="If given, the two strings must be W0j and G20i for the integral constraint")

  return p.parse_args()


def residuals(params, func, x, y, theo, win, errs):
  """
  This function returns the difference between a measured and a model power spectrum

  Parameters
  ----------
  params: list
    list of model marameters
  func: function
    function that returns the model power spectrum given x, params, theo, win
  x: array
    value of k of the measured power spectrum
  y: array
    value of P(k) of the measured power spectrum
  theo: list
    list of arrays used to build the model power spectrum
  win: list
    list of arrays containing the files to convolve the model with the window 
    function wij, kj, W0j, G02i. If the latter two are given the integral 
    contraint is also used with the convolution
  errs: array
    errors in the measurement

  ---------
  output: array
    different between the measured and the model power spectrum
  """
  predicted = func(x, params, theo, win)
  if errs is None:
    return y - predicted
  else:
    return np.divide(y - predicted, errs)

if __name__ == "__main__":   #if it's the main

  opt, args = options(op.OptionParser(version="%prog version 0.1"))   #create the object optparse
 
  #read the measured power spectrum
  measured = np.loadtxt(args[0])
  #check the range of scales
  krange = (measured[:,0] >= float(args[1])) & (measured[:,0] <= float(args[2]))
  measured = measured[krange,:]

  if(opt.isvar == False):
    errs = None
  else:
    if(opt.cov == None):
      errs = measured[:,2]
    else:
      if(opt.covcol == None):
	errs = np.sqrt( np.diag(np.loadtxt(opt.cov))[krange] )   #if covariance is already square in the file
      else:  #if the covariance is on one column of the file, reshape it and get the variance
	errs = np.loadtxt(opt.cov, usecols=[opt.covcol,])
	size = int(np.sqrt(errs.size))
	errs = np.sqrt( np.diag(errs.reshape([size,size]) )[krange] )

  #read the theory power spectra to build the model and save in a list
  theo= [np.loadtxt(opt.flin), np.loadtxt(opt.f1loop)]
  theo[0] = np.vstack( (theo[1][:,0], np.interp( theo[1][:,0], theo[0][:,0], theo[0][:,1] )) ).T  #interpolate the values of the linear power spectrum in the k values of the 1loop
  
  #if the window matrix files are given read them and save in a list
  win = [None, ]*4
  if(opt.win):
    win[0] = np.loadtxt(opt.win[0])[krange,:]  #Wij  in the range of 'i' required
    win[1] = np.loadtxt(opt.win[1])   #kj
  if(opt.win_extr):
    win[2] = np.loadtxt(opt.win_extr[0])   #W0j
    win[3] = np.loadtxt(opt.win_extr[1])   #G02i

  params, covar, info, mesg, ier = spo.leastsq(residuals, opt.guess,
      args = (pf.model, measured[:,0], measured[:,1], theo, win, errs), full_output = True)
      
  print(ier, mesg)
  print("best chi2 paramaters: ", params)

  exit()
