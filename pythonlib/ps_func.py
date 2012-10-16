#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Module containing some function useful when dealing with power spectra and window functions"""

import itertools as it
import my_statistic as ms
import numpy as np
import scipy.integrate as spi

def shell_volume(k, dk=None):
  """Approximated formula for the volume in fourier space for small bin of width dk around k
  k^2*delta_k/2*pi^2

  Parameters
  ----------
  k: 1D array
    k values of the bins 
  dk: float (optional)
    width of the bins, if not given computed from k

  output: 1D array with the same shape of k
    volume of a spherical shell centered in k and of width dk
  """
  if(dk==None):
    dk = k[1]-k[0]

  return (k*k*dk) / 2.*np.pi*np.pi


def meansigma_k_shell(k, sigma):
  """Performs the average of function 'sigma' in spherical shells between k[i] and k[i+i]:
  N int_k[i]^k[i+1] dk' k'^2 sigma^2(k'). N=18*pi^2/(k[i+i]^3 - k[i]^3)^2

  Parameters
  ----------
  k: 1D array
    N+1 values of extremes of the N bins in which the mean is performed
  sigma: 2D array
    array containing the value of sigma to be integrated with k' and sigma(k') in the first 2 columns

  output: 1D array
    average of sigma in the k-shels
  """

  N = 18.*np.pi*np.pi/np.power(k[1:]*k[1:]*k[1:] - k[:-1]*k[:-1]*k[:-1] , 2.)   #normalisation factor
  ms = np.empty_like(k[1:])  #create the output array with 1 element less than k
  ks, k2s2 = sigma[:,0], sigma[:,0]*sigma[:,0]*sigma[:,1]*sigma[:,1]
  func = lambda x,kk,s: np.interp(x, kk, s, left=0, right=0)   #integration function
  for i,km, kp in it.izip(it.count(), k[:-1], k[1:]):  #integrate in the bins
    ms[i] = spi.quad(func, km, kp, args=(ks,k2s2))[0]
  
  avs = N*ms    #spherical average of sigma^2

  return avs


def mymodel(plin, p1loop, k, kstar=0.15, AMC=1., b=1.):
  """Given a linear and a 1loop power spectrum, the function interpolates them in k
  and the return the model power spectrum: 
  P(k) = b^2(exp(-(k/kstar)^2)*Plin(k) + AMC*P1loop(k))
  If k exceed the limits of the input power spectra, interpolation padd with 0

  Parameters
  ----------
  plin, p1loop: 2D array
    linear and 1loop power spectra with k and P(k) in the first and second column
  k: 1D array
    values of k where to evaluate the model
  kstar, AMC, b: floats (optional)
    model parameters

  output: (k.size,2) array
    model power spectrum
  """
  iplin = np.interp(k, plin[:,0], plin[:,1], left=0., right=0.)     #interpolated linear P(k)
  ip1loop = np.interp(k, p1loop[:,0], p1loop[:,1], left=0, right=0)     #interpolated 1loop P(k)
  model = b * (np.exp(-np.power(k/kstar,2.))*iplin + AMC*ip1loop)   #model power spectrum

  return np.vstack((k,model)).T

def bias(measure, model, var=None, kmin=None, kmax=None):
  """given the the power spectrum and the reference one, returns the square of the bias
  
  Parameters
  ----------
  measure: 2D array
    measured power spectrum: measure[:,0]:k; measure[:,1]: P(k)
  model: 1D array
    model for whom the bias must be computed, interpolated in measure[:,0]
  var: 1D array (optional)
    variance of 'measure'
  kmin: double (optional)
    minimum k for the fit
  kmax: double (optional)
    maximum k for the fit
  
  ----------
  output: float
    bias to associate with the model
  """

  if(var==None):  #set the variance to 1
    var = np.ones_like(measure[:,0])
  else:
    if(var.size != measure.shape[0]):
      raise SystemExit("The dimention of the variance must agree with the number of rows in 'measure'")
  k,pk = measure[:,0], measure[:,1]
  if(kmin == None):   #set the min and max k if not given
    kmin = k[0]
  if(kmax == None):
    kmax = k[-1]

  imin, imax = (k<kmin).sum(), (k<kmax).sum()   #minimum and maximum index of the fitting range
  sigma = var[imin:imax]*var[imin:imax]   #variance square
  dt_s2 = np.sum(pk[imin:imax]*model[imin:imax]/sigma)     #sum data*theory/sigma^2
  t2_s2 = np.sum(model[imin:imax]*model[imin:imax]/sigma)  #sum theory^2/sigma^2

  return dt_s2/t2_s2   #return the bias

def convolve(ps, wij, kj, w0j=None, g20i=None):
  """Given a power spectrum and the four quantities related to the window matrix,
  returns the convolved power spectrum.
  if w0j and/or g20i are not given the correction due to integral constraint is not done
  
  Parameters
  ----------
  ps: 2D array
    power spectrum to convolve: ps[:,0]=k; ps[:,1]=P(k)
  wij: 2D array
    window matrix
  kj: 1D array
    values of k of the 'j' dimention of the window matrix
  w0j: 1D array  (optional)
    window matrix at k_i=0
  G20i: 1D array  (optional)
    window function evaluated in k_i=0 and for the k of the 'i' dimention of the window matrix

  ----------
  output: 1D array
    convolved ps
  """

  intps = np.interp(kj, ps[:,0], ps[:,1], right=0., left=0.)   #interpolate the power spectrum
  conv = np.dot(wij, intps)     #do matrix multiplication
  if(w0j != None and g20i != None):
    conv -= np.sum(w0j*intps)*g20i[1:]/g20i[0] #and integral constraint

  return conv

def norm(win, gfit=0):
  """Given an array of at least two columns, computed the area under the curve discribed by the first two columns.

  Parameters
  ----------
  win: 2-D array
    array function to integrate

  gfit: int (optional)
    if gfin>0, win[:gfin,:2] used for a gaussian fit and then the result is used to integrate between 0 and win[0,0]

  ----------
  output: float
    normalisation of the window function
  """
  N = spi.trapz(win[:,1], x=win[:,0]) 
  if(gfit > 0):
    a,m,s = ms.gausfit(win[:int(gfit),0], win[:int(gfit),1])  #gaussian fit of the very large scale window function
    N += spi.quad(lambda x: a*np.exp(-np.power(x-m,2.)/(2.*s*s)), 0, win[0,0])[0]

  return N

def kaiser_boost(beta):
  """Given beta=f/b, with f=d ln D/d ln a, b the linear bias and D the linear growth factor,
  return the linear kaiser boost: S = 1 + 2/3 beta + 1/5 beta^2
  Parameters
  ----------
  beta: float or iterable
    f/b
  output: 
    S = 1 + 2/3 beta + 1/5 beta^2
  """
  return 1. + 2.*beta/3. + 1.*beta*beta/5.

def realspace_b(b2_S, fgrowth):
  """Given the 'bias' in redshift space b_s^2 and the logarithmic derivative 
  of the linear growth factor D with respect to the scale factor a, extract the 
  linear real space bias
  Parameters
  ----------
  b2_S: float
    redshift space amplitude w.r.t. dark matter in real space
  fgrowth: float
    logarithmic derivative of the linear growth factor with respect to the scale factor
  output: float
    real space linear bias
  """

  return np.sqrt(b2_S - 4.*fgrowth*fgrowth/45.) - fgrowth/3.

def model(x, params, theo, win=[None,]*4):
  """
  This function returns the model given set of parameters, the theory power spectra
  and the files of the window function

  Parameters
  ----------
  x: array
    value of k of the measured power spectrum
  params: list
    list of model marameters: bias, kstar, amc, ...
  theo: list
    list of power spectra used to build the model power spectrum: plin, p1loop
    they must be evaluated in the same position
  win: list (optional)
    list of arrays containing the files to convolve the model with the window 
    function wij, kj, W0j, G02i. If the latter two are given the integral 
    contraint is also used with the convolution

  ---------
  output: array
    model
  """

  k = theo[0][:,0] #k where to evaluate the model
  mod = (np.exp(-k*k/(params[1]*params[1])) * theo[0][:,1] + params[2] * theo[1][:,1] )  #model
  mod *= params[0]*params[0]  #multiply by the bias

  if(win[1] == None):  #if the window matrix is not given interpolate the theory at 'x'
    mod = np.interp(x, k, mod)
  else:   # convolve with the window function and cut it at the required values of x
    mod = convolve( np.vstack((k,mod)).T, win[0], win[1], w0j=win[2], g20i=win[3] )

  return mod 

def rsd( k, f, b, sigma_v):
  """Spherical averaged redshift space distortion on the monopole of the power spectrum.
  Parameters
  ----------
  k: float or array
    wavenumbers where to evaluate the function
  f: float
    growth rate
  b: float
    bias
  sigma_v: float
    variance of the velocity field
  output
  ------
  RSD: same as k
    redshift space distortion effect
  """
    
  beta = f/b
  x = k * f * sigma_v
  RSD = np.arctan(x)/x * ( (beta+x*x) **2 - 4.*beta*beta) + 2.*beta*beta
  RSD += (beta-x*x)**2 / (x*x+1.)
  RSD /= 2.* x**4

  return RSD
