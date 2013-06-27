# -*- coding: utf-8 -*-
"""
This module implements some function with non available within standard python packages
"""

import itertools as it
import my_functions as mf
import numpy as np    # import numpy
import scipy.stats as sps   #statistical funcitons from scipy
import sys

def append_weights(a, w):
  """
  Given the 2D array 'a' and a 1d array of weights 'w' extend a such
  that every row of a is repeated n times according to the values is w.
  The repetition is done appending the new values along axis=0

  Parameters
  ----------
  a: array like
    1D or 2D array
  w: array like
    1D array containg the values the weights of the elements in the columns of 'a'.
    Must have the same dimention as the 0th of 'a'.

  output: array like
    contains 'a' with the rows replicated according to the correponding weights
  """

  if(w.shape[0] != a.shape[0]):   #check that they are fine
    print "ERROR: the lengh of 'w' does not correspond to the one of the 0th column of a"
    sys.exit(100)

  return a.repeat(np.array(w, dtype='i'), axis=0)  #convert w to int and repeat elements of a according to the values in w


def stddev(a, weights=None, axis=None):
  """
  Compute the standard deviation of the array a along a given axis

  Parameters
  ----------
  a: array like
    contains the values whose standard deviations are desired
  weights: array like, optional
    array containg the values the weights of the elements in the columns of 'a'.
    Must have the same dimention as the 0th of 'a' or the same dimention of 'a'.
  axis: integer, optional
    Axis over which the sum is taken. By default axis is None, and all elements are summed.

  output: array like
    standard deviations
  """

  if(weights != None):  #if the weights are given used hand made recipe (appending the weights can create an array too large for numpy.std)
    if(a.shape != weights.shape):
      if(a.shape[0] != weights.shape[0] and len(weights.shape)==1):
        print "ERROR: the code needs either weights.shape=a.shape or weights.shape[0] = a.shape[0]"
	sys.exit(110)
      else:
        weights = np.column_stack([weights,]*a.shape[1])   #so weights.shape=a.shape
    swx2 = np.sum(weights*a*a, axis=axis)
    sw = np.sum(weights, axis=axis)
    sw2 = np.sum(weights*weights, axis=axis)
    swx = np.sum(weights*a, axis=axis)
    std = (swx2*sw - swx*swx) / (sw*sw - sw2)
    if( std.ndim == 0 ):   #if std is a single element convert to a 1d array in order to continue
      std = np.array( [std,] )
    for i,s,w,x in it.izip(it.count(),mf.convert2array(std), mf.convert2array(sw), mf.convert2array(swx)):
      if(abs(s*w/x) > 1e-09):
	std[i] = np.sqrt(s)
      else:
	std[i] = 0.
  else:
    std = np.std(a, axis=axis)  #if no weights used

  return std

def cov(a, weights=None):
  """
  Compute the covariance: each column is a parameters and each row a realisation 

  Parameters
  ----------
  a: array like
    contains the values whose standard deviations are desired
  weights: array like, optional
    1D array containg the values the weights of the elements in the columns of 'a'.
    Must have the same dimention as the 0th of 'a' or the same dimention of 'a'.

  output: array like
    covariance
  """

  if(weights != None):
    if(a.shape != weights.shape):
      if(a.shape[0] != weights.shape[0] and len(weights.shape)==1):
        print "ERROR: the code needs either weights.shape=a.shape or weights.shape[0] = a.shape[0]"
	sys.exit(110)
      else:
        weights = np.column_stack([weights,]*a.shape[1])   #so weights.shape=a.shape
    


  return cov


def percentile(a, weights=None, perc=np.array([16.,84.])):
  """Compute the percentile along the 0th axis

  Parameters
  ----------
  a: array like
    (1D or 2D) containing the values whose percentile are desired
  weights: array like, optional
    1D array containg the values the weights of the elements in the columns of 'a'.
    Must have the same dimention as the 0th of 'a'. If 'None', weights=1
  perc: 1D numpy array
    list of floats containing the required percentiles
    they will be ordered from smaller to bigger

  output: array like
    array of size len(perc) * a.shape[1]+1, the 'i'th rows containing the values 
    corresponding to the 'i'th the percentile [WARNING: the percentiles are ordered]
    The 0th column contains the percentile values
    If the percentile falls between two values, linear interpolation between the nearest values will be performed
  """

  if(weights != None):  #if the weights are given
    a = append_weights(a, weights)

  perc.sort()  #sort the percentile scores

  output = np.zeros((perc.shape[0],a.shape[1]+1))  #create the output array: output[:,1:] correspond to the columns of 'a', the rows are of the size of 'perc' ordered from the lower to the upper score
  output[:,0] = perc   #output[:,0] contains informations about the percentile scores

  for i,p in enumerate(perc):   #go through the various percentiles
    output[i,1:] = sps.scoreatpercentile(a, p)  #and compute the corresponding score

  return output

def gaussfit(x,y, mu_zero=False):
  """Given x and y computes the gaussian fit, 
  transforming the gaussian into a polinomial
  "p(x) = p[0] * x**deg + ... + p[deg]"
   and using numpy.polyfit

  Parameters
  ----------
  x: 1-D array
    x coordinates
  y: 1-D array
    y = f(x)
  mu_zero: bool (optional)
    if true force the mean of the gaussian to 0

  output: A, mu, sigma
    such that f(x) = A*exp(-(x-mu)^2/(2*sigma^2))
  """

  xf, yf = x[y>0.], y[y>0.]  #take away the negative values of y
  deg = 2  #number of degree of the polinomial when mu, sigma and A are to fit
  if( mu_zero == True ):  #if mu forced to 0
    deg = 1  #only linear fit required y=p[0] * t + p[1]
    xf = xf**2    #with t = x**2
  p = np.polyfit(xf, np.log(yf), deg)   #execute the polinomial fit for a parabola from ln(y) = f(x) 

  sigma = np.sqrt(-1./(2.*p[0]))     # sigma^2
  if( mu_zero == True ):  #if mu forced to 0
    mu = 0.        # x mean
    A = p[1] #ln(amplitude)
  else:
    mu = p[1] * sigma*sigma        # x mean
    A = p[2] + mu*mu/(2.*sigma*sigma)  #ln(amplitude)
  
  return np.exp(A), mu, sigma

def GR_criterion(chains):
  """
  Compute the Geldman and Rubin criterion from a list of MCMC chains

  Parameters
  ----------
  chains: list
    list of arrays, each one containing a single MCMC chain with the structure: weight, likelihood, parameters

  optput: 1D array
    list of the R values for the parameters
  """
  n_step = np.empty(len(chains))   #it contains the number of step for each chain
  for i,c in enumerate(chains):   #count the number of steps
    n_step[i] = np.sum(c[:,0])
  n_step = np.amin(n_step)   #get the number of steps in the shortest chain

  ch_mean, ch_var = np.empty((len(chains), chains[0].shape[1]-2)), np.empty((len(chains), chains[0].shape[1]-2))  #it will contain the various in-chain mean and variances^2 for the parameters
  for i,c in enumerate(chains):
    tc = append_weights(c[:,2:], c[:,0])[-n_step:,:]  #replicate weights and cut all the chains such that they are all the same length
    ch_mean[i,:] = np.average(tc, axis=0)     #compute the in-chain average
    ch_var[i,:] = np.power(stddev(tc, axis=0), 2.)    #compute the in-chain square of the variances

  W = np.average(ch_var, axis=0)            #mean of the square of the variances 
  B = np.power(stddev(ch_mean, axis=0), 2.) #quare variance of the mean of the chains
  var = (1.-1./n_step)*W + B     #estimated variance
  old_err_state = np.seterr(divide='raise')   #change error handling
  try:
    R = np.divide(var,W)      #square of the scale reduction factor
  except FloatingPointError:   #avoid the printing of warnings
    pass

  return np.sqrt(R)

