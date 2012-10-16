#!/usr/bin/python

import itertools as it
import my_functions as mf
import numpy as np
import scipy.constants as spc   #physical constants
import scipy.integrate as spi
import setcosmology

ckms = spc.c * 1e-3   #speed of light in km/s

class Distance(object):
  """Class containing cosmological distances and derived quantities
  scosmo is a 'Setcosmology' object
  """

  def __init__(self, scosmo):
    """Initialization of the class
    'scosmo' is a 'Setcosmology' object
    """
    self.c = scosmo   #assigne the cosmology object to self.c

  def comoving_distance_a(self, a):
    """Compute the comoving distance xi(a) = c*\int_a^1 da'/(a'^2H(a')) in Mpc
    Parameters
    ----------
    a: float, list, tuble, numpy array
      scale factor normalised to 1 at present time.

    output: float
      comoving distance evaluated in a
    """
    a = mf.convert2array(a)
    cdis = np.empty_like(a)   #initialise the array containing the comoving distance
    a2H = lambda t: 1./(np.power(t,2.) * self.c.H_a(t))     #define the integrand as a local function. Integration in logarith
    for i,aa in enumerate(a):
      cdis[i] = spi.quad(a2H, aa, 1.)[0]
    return ckms * cdis   #return the comoving distance
  def comoving_distance_ah(self, a):
    """Compute the comoving distance xi(a) = c*\int_a^1 da'/(a'^2H(a')) in Mpc/h
    Parameters
    ----------
    a: float, list, tuble, numpy array
      scale factor normalised to 1 at present time.

    output: float
      comoving distance evaluated in a
    """
    return self.comoving_distance_a(a)*self.c.h0
  def comoving_distance_z(self, z):
    """Compute the comoving distance xi(z) = c*\int_0^z dz'/(H(z')) in Mpc
    Parameters
    ----------
    z: float, list, tuble, numpy array
      redshift

    output: float
      comoving distance evaluated in z
    """
    z = mf.convert2array(z)
    cdis = np.empty_like(z)   #initialise the array containing the comoving distance
    oneoH = lambda t: 1./ self.c.H_z(t)     #define the integrand as a local function. Integration in logarith
    for i,zz in enumerate(z):
      cdis[i] = spi.quad(oneoH, 0., zz)[0]
    return ckms * cdis   #return the comoving distance
  def comoving_distance_zh(self, z):
    """Compute the comoving distance xi(z) = c*\int_0^z dz'/(H(z')) in Mpc/h
    Parameters
    ----------
    z: float, list, tuble, numpy array
      redshift

    output: float
      comoving distance evaluated in z
    """
    return self.comoving_distance_z(z)*self.c.h0

  def angular_distance_a(self, a):
    """Compute the angular diameter distance d_a=a*r(xi(a)) in Mpc with 
    r(xi) = xi, if Omega_k=0
    r(xi) = sin(sqrt(-Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k<0
    r(xi) = sinh(sqrt(Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k>0
    Parameters
    ----------
    a: float
      scale factor normalised to 1 at present time.

    output: float
      angular diameter distance evaluated in a
    """
    a = mf.convert2array(a)
    xi = self.comoving_distance_a(a) / ckms   #compute the comoving distance without c
    if(self.c.ok == 0.):   #if flat
      r = xi
    elif(self.c.ok < 0.):    #if closed
      r = np.sin(np.sqrt(-self.c.ok) * self.c.H0 * xi) / (self.c.H0 * np.sqrt(np.abs(self.c.ok)))
    else:   #if open
      r = np.sinh(np.sqrt(self.c.ok) * self.c.H0 * xi) / (self.c.H0 * np.sqrt(np.abs(self.c.ok)))
    return a * ckms *r   #angular diameter distance
  def angular_distance_ah(self, a):
    """Compute the angular diameter distance d_a=a*r(xi(a)) in Mpc/h with 
    r(xi) = xi, if Omega_k=0
    r(xi) = sin(sqrt(-Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k<0
    r(xi) = sinh(sqrt(Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k>0
    Parameters
    ----------
    a: float
      scale factor normalised to 1 at present time.

    output: float
      angular diameter distance evaluated in a
    """
    return self.angular_distance_a(a)*self.c.h0
  def angular_distance_z(self, z):
    """Compute the angular diameter distance d_z=r(xi(z))/(1+z) in Mpc with 
    r(xi) = xi, if Omega_k=0
    r(xi) = sin(sqrt(-Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k<0
    r(xi) = sinh(sqrt(Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k>0
    Parameters
    ----------
    z: float, list, tuble, numpy array
      redshift

    output: float
      angular diameter distance evaluated in z
    """
    a = 1./(1+mf.convert2array(z))
    return self.angular_distance_a(a)
  def angular_distance_zh(self, z):
    """Compute the angular diameter distance d_z=r(xi(z))/(1+z) in Mpc/h with 
    r(xi) = xi, if Omega_k=0
    r(xi) = sin(sqrt(-Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k<0
    r(xi) = sinh(sqrt(Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k>0
    Parameters
    ----------
    z: float, list, tuble, numpy array
      redshift

    output: float
      angular diameter distance evaluated in z
    """
    return self.angular_distance_z(z)*self.c.h0

  def luminosity_distance_a(self, a):
    """Compute the luminosity distance d_a=r(xi(a))/a in Mpc with 
    r(xi) = xi, if Omega_k=0
    r(xi) = sin(sqrt(-Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k<0
    r(xi) = sinh(sqrt(Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k>0
    Parameters
    ----------
    a: float
      scale factor normalised to 1 at present time.

    output: float
      luminosity distance evaluated in a
    """
    a = mf.convert2array(a)
    return self.angular_distance_a(a)/(a*a)
  def luminosity_distance_ah(self, a):
    """Compute the luminosity distance d_a=r(xi(a))/a in Mpc/h with 
    r(xi) = xi, if Omega_k=0
    r(xi) = sin(sqrt(-Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k<0
    r(xi) = sinh(sqrt(Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k>0
    Parameters
    ----------
    a: float
      scale factor normalised to 1 at present time.

    output: float
      luminosity distance evaluated in a
    """
    return self.luminosity_distance_a(a)*self.c.h0
  def luminosity_distance_z(self, z):
    """Compute the luminosity distance d_z=r(xi(z))*(1+z) in Mpc with 
    r(xi) = xi, if Omega_k=0
    r(xi) = sin(sqrt(-Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k<0
    r(xi) = sinh(sqrt(Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k>0
    Parameters
    ----------
    z: float, list, tuble, numpy array
      redshift

    output: float
      luminosity distance evaluated in z
    """
    a = 1./(1+mf.convert2array(z))
    return self.angular_distance_a(a)/(a*a)
  def luminosity_distance_zh(self, z):
    """Compute the luminosity distance d_z=r(xi(z))*(1+z) in Mpc/h with 
    r(xi) = xi, if Omega_k=0
    r(xi) = sin(sqrt(-Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k<0
    r(xi) = sinh(sqrt(Omega_k)H_0*xi)/(H_0*sqrt(|Omega_k|))   if Omega_k>0
    Parameters
    ----------
    z: float, list, tuble, numpy array
      redshift

    output: float
      luminosity distance evaluated in z
    """
    return self.luminosity_distance_z(z)*self.c.h0

  def effective_distance_a(self, a):
    """Compute the effective distance d_V(a) = [d_a^2(a) *c*(1./a-1)/ (H(a))]^(1/3) in Mpc
    Parameters
    ----------
    a: float
      scale factor normalised to 1 at present time.

    output: float
      effective distance evaluated in a
    """
    a = mf.convert2array(a)
    z = 1./a - 1.   #get the redshift
    d_a = self.angular_distance_a(a)   #get the angular diameter distance
    dv3 = d_a*d_a*z*ckms / self.c.H_a(a)    #cube of the effective distance
    return np.power(dv3, 1./3.)      #return the effective distance
  def effective_distance_ah(self, a):
    """Compute the effective distance d_V(a) = [d_a^2(a) *c*(1./a-1)/ (H(a))]^(1/3) in Mpc/h
    Parameters
    ----------
    a: float
      scale factor normalised to 1 at present time.

    output: float
      effective distance evaluated in a
    """
    a = mf.convert2array(a)
    z = 1./a - 1.   #get the redshift
    d_a = self.angular_distance_ah(a)   #get the angular diameter distance
    dv3 = d_a*d_a*z*ckms / (np.sqrt(self.c.E_a(a))*100.)    #cube of the effective distance
    return np.power(dv3, 1./3.)      #return the effective distance
  def effective_distance_z(self, z):
    """Compute the effective distance d_V(z) = [d_a^2(z) *c*z/ (H(z))]^(1/3) in Mpc
    Parameters
    ----------
    z: float, list, tuble, numpy array
      redshift

    output: float
      effective distance evaluated in z
    """
    z = mf.convert2array(z)
    d_a = self.angular_distance_z(z)   #get the angular diameter distance
    dv3 = d_a*d_a*z*ckms / self.c.H_z(z)    #cube of the effective distance
    return np.power(dv3, 1./3.)      #return the effective distance
  def effective_distance_zh(self, z):
    """Compute the effective distance d_V(z) = [d_a^2(z) *c*z/ (H(z))]^(1/3) in Mpc/h
    Parameters
    ----------
    z: float, list, tuble, numpy array
      redshift

    output: float
      effective distance evaluated in z
    """
    z = mf.convert2array(z)
    d_a = self.angular_distance_zh(z)   #get the angular diameter distance
    dv3 = d_a*d_a*z*ckms / (np.sqrt(self.c.E_z(z))*100.)    #cube of the effective distance
    return np.power(dv3, 1./3.)      #return the effective distance

  def effective_volume_sr_a(self, a):
    """Compute the effective volume in [Mpc]^3/deg^2
    \int_af^ai dV(a') da'
    with dV(a) = c*d_a^2/(a^4*H(a))
    if len(a)==1 the effective volume the integration is done 
    in the interval a'=[a,1]; otherwise the values of a given 
    are interpreted as the extrema of integration interval and 
    the effective volume of the len(a)-1 intervals is returned
    Parameters
    ----------
    a: float, list, tuble, numpy array
      scale factor normalised to 1 at present time.

    output: float
      effective distance evaluated in a
    """
    a = mf.convert2array(a)
    if(a.size == 1):  #if only one z given, add 0
      a = np.r_[a,0.]
    eff_vol = np.empty(a.size-1)  #initialise the output array
    integrand = lambda t: self.angular_distance_a(t)*self.angular_distance_a(t)/(t*t*t*t * self.c.H_a(t))     #define the integrand as a local function.
    for i, ai, af in it.izip(it.count(), a[:-1], a[1:]):  #do the integration in each redshift bin between z[i] and z[i+1]
      eff_vol[i] = spi.quad(integrand, ai, af)[0]
    return ckms*eff_vol
  def effective_volume_sr_ah(self, a):
    """Compute the effective volume in [Mpc/h]^3/deg^2
    \int_af^ai dV(a') da'
    with dV(a) = c*d_a^2/(a^4*H(a))
    if len(a)==1 the effective volume the integration is done 
    in the interval a'=[a,1]; otherwise the values of a given 
    are interpreted as the extrema of integration interval and 
    the effective volume of the len(a)-1 intervals is returned
    Parameters
    ----------
    a: float, list, tuble, numpy array
      scale factor normalised to 1 at present time.

    output: float
      effective distance evaluated in a
    """
    a = mf.convert2array(a)
    if(a.size == 1):  #if only one z given, add 0
      a = np.r_[a,0.]
    eff_vol = np.empty(a.size-1)  #initialise the output array
    integrand = lambda t: self.angular_distance_ah(t)*self.angular_distance_ah(t) / (t*t*t*t * np.sqrt(self.c.E_a(t)))     #define the integrand as a local function.
    for i, ai, af in it.izip(it.count(), a[:-1], a[1:]):  #do the integration in each redshift bin between z[i] and z[i+1]
      eff_vol[i] = spi.quad(integrand, ai, af)[0]
    return ckms*eff_vol/100.
  def effective_volume_sr_z(self, z):
    """Compute the effective volume in [Mpc]^3/deg^2
    \int_zi^zf dV(z') dz'
    with dV(z) = c*(1+z)^2*d_a^2/H(z)
    if len(z)==1 the effective volume the integration is done 
    in the interval z'=[0,z]; otherwise the values of z given 
    are interpreted as the extrema of integration interval and 
    the effective volume of the len(z)-1 intervals is returned
    Parameters
    ----------
    z: float, list, tuble, numpy array
      redshift

    output: float
      effective distance evaluated in z
    """
    z = mf.convert2array(z)
    if(z.size == 1):  #if only one z given, add 0
      z = np.r_[0.,z]
    eff_vol = np.empty(z.size-1)  #initialise the output array
    integrand = lambda t: (1+t)*(1+t)*self.angular_distance_z(t)*self.angular_distance_z(t)/self.c.H_z(t)     #define the integrand as a local function.
    for i, zi, zf in it.izip(it.count(), z[:-1], z[1:]):  #do the integration in each redshift bin between z[i] and z[i+1]
      eff_vol[i] = spi.quad(integrand, zi, zf)[0]
    return ckms*eff_vol
  def effective_volume_sr_zh(self, z):
    """Compute the effective volume in [Mpc/h]^3/deg^2
    \int_zi^zf dV(z') dz'
    with dV(z) = c*(1+z)^2*d_a^2/H(z)
    if len(z)==1 the effective volume the integration is done 
    in the interval z'=[0,z]; otherwise the values of z given 
    are interpreted as the extrema of integration of each redshift
    interval and the effective volume of the len(z)-1 intervals is
    returned
    Parameters
    ----------
    z: float, list, tuble, numpy array
      redshift

    output: float
      effective distance evaluated in z
    """
    z = mf.convert2array(z)
    if(z.size == 1):  #if only one z given, add 0
      z = np.r_[0.,z]
    eff_vol = np.empty(z.size-1)  #initialise the output array
    #define the integrand as a local function.
    integrand = lambda t: (1+t)*(1+t) * self.angular_distance_zh(t)*self.angular_distance_zh(t) / np.sqrt(self.c.E_z(t))
    for i, zi, zf in it.izip(it.count(), z[:-1], z[1:]):  #do the integration in each redshift bin between z[i] and z[i+1]
      eff_vol[i] = spi.quad(integrand, zi, zf)[0]
    return ckms*eff_vol/100.
