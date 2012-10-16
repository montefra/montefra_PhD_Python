#!/usr/bin/python

import itertools as it
import my_functions as mf
import numpy as np
import scipy.constants as spc   #physical constants
import scipy.integrate as spi
import sys 


class Cosmology(object):
  """Class containing cosmological functions. Version 1.
  "cosmo" is a dictionary with the cosmological parameters wanted
  """

  def __init__(self, **cosmo):
    """Initialization of the class. 
    Keyword accepted:
    "om", "omegam", "OmegaM" or "Omegam":
      Omega matter [default: 0.25]
    "ok", "omegak" or "Omegak":
      Omega curvature [default: 0]. Convention Omega_k>0: close; Omega_k>0: open
    "ob", "omegab", "Omegab":
      Omega baryons [default: 0.04]
    "ora", "omegar", "Omegar":
      Omega radiation [default: 0]
    "wl", "wL", "wde", "wDE", "w0":
      Dark energy equation of state at a=1 (z=0) [default: -1]
    "wa":
      derivative of wDE(a) for model wDE(a) = w0 + wa(1-a) [default: 0]
    "h", "h0", "H", "H0":
      Hubble parameter ("H", "H0") or reduced Hubble parameters ("h", "h0") [default: h=0.74
    """
    keys = sorted(cosmo.keys())   #get the keywords
    self._set_default()
    for k in keys:    #assigne the keywords. If one is not give, default set
      if k in ("om", "omegam", "OmegaM", "Omegam"):
        self.om = cosmo[k]
      elif k in ("ok", "omegak", "Omegak"):
        self.ok = cosmo[k]
      elif k in ("ob", "omegab", "Omegab"):
        self.ob = cosmo[k]
      elif k in ("ora", "omegar", "Omegar"):
        self.ora = cosmo[k]
      elif k in ("wl", "wL", "wde", "wDE", "w0"):
        self.w0 = cosmo[k]
      elif k in ("wa"):
        self.wa = cosmo[k]
      elif k in ("h", "h0", "H", "H0"):
        if k in ("h", "h0"):
	  self.h0 = cosmo[k]
	else:
	  self.h0 = cosmo[k]/100.
      elif k in ("odm", "omegadm", "OmegaDM", "Omegadm", "ol", "ode", "omegal", "omegaL", "omegaDE", "Omegal", "Omegade", "OmegaL", "OmegaDE"):
        print "The parameter '%s' is derived from the other parameters given" %k
      else:
        print "The parameter '%s' is of no interests now" %k

    self._set_derived()
        
  def _set_default(self):
    """Set the default cosmological parameters"""
    self.om = 0.25   #omega matter
    self.ok = 0.     #omega curvature
    self.ob = 0.04   #omega baryon
    self.ora = 0.    #omega radiation
    self.w0 = -1.    #dark energy equation of state 
    self.wa = 0.     #scale factor defivative of the DE equation of state (w(a) = wde + wa (1-a))
    self.h0 = 0.74   #reduced hubble parameter

  def _set_derived(self):
    """Set derived parameters"""
    self.ode = 1 - self.ok - self.om - self.ora   #Omega dark energy: Om+Ode+Or+Ok = -1
    self.odm = self.om - self.ob       #Omega dark matter: Odm+Ob=Om
    self.H0 = self.h0*100.             #Hubble parameter

  def get_cosmology(self):
    "return all the cosmological parameters"
    return dict(om=self.get_om(), odm=self.get_odm(), ob=self.get_ob(), ok=self.get_ok(), ode=self.get_ode(), ora=self.get_ora(), w0=self.get_w0(), wa=self.get_wa(), h0=self.get_h0(), H0=self.get_H0())

  def get_om(self):
    "return Omega matter"
    return self.om
  def get_odm(self):
    "return Omega dark matter"
    return self.odm
  def get_ob(self):
    "return Omega baryon"
    return self.ob
  def get_ok(self):
    "return Omega curvature"
    return self.ok
  def get_ode(self):
    "return Omega dark energy"
    return self.ode
  def get_ora(self):
    "return Omega radiation"
    return self.ora
  def get_w0(self):
    "return dark energy equation of state at a=1 (z=0)"
    return self.w0
  def get_wa(self):
    "return the scale factor derivative of the dark energy equation of state w(a) = w0 + wa(1-a)"
    return self.wa
  def get_h0(self):
    "return the reduced Hubble parameter"
    return self.h0
  def get_H0(self):
    "return the Hubble parameter"
    return self.H0

  def print_cosmology(self):
    """print the cosmological parameters in a nice format"""
    paramnames = dict(om = "Omega_{matter}", 
                  odm = "Omega_{dark matter}", 
		  ob = "Omega_{baryons}", 
		  ok = "Omega_{curvature}", 
		  ode = "Omega_{dark energy}", 
		  ora = "Omega_{radiation}", 
		  w0 = "w_{dark energy}(a=1)", 
		  wa = "d(w_de)/da(a=1)", 
		  h0 = "reduced Hubble parameter h_0", 
		  H0 = "Hubble parameter H_0 [km/(s Mpc)]")
    paramvalues = self.get_cosmology()
    for k in paramnames.keys():
      print "{0:40} {1:4.3e}".format(paramnames[k]+":", paramvalues[k])

  def _de_ev(self,a):
    """Given a set of cosmological parameters computes the 
    evolution factor of dark energy assuming wDE = w0 + wa(1-a) 
    Parameters
    ----------
    a: numpy array
      scale factor normalised to 1 at present time
    output: np.array with a.shape
      a^(-3(1+w0+wa))exp(-3wa(1-a))
    """
    return np.power(a, -3.*(1.+self.w0+self.wa)) * np.exp(-3.*self.wa*(1.-a)) #scale factor evolution of Omega_DE in a wDE = w0 + wa(1-a) model

  def E_a(self, a):
    """Given a set of cosmological parameters computes H^2(a)/H_0^2
    Parameters
    ----------
    a: float, list, tuble, numpy array
      scale factor normalised to 1 at present time
    output: np.array with a.shape
      E(a)
    """
    a = mf.convert2array(a)
    return self.om/(a*a*a) + self.ora/(a*a*a*a) + self.ode*self._de_ev(a) + self.ok/(a*a)

  def _dE_da(self, a):
    """Given a set of cosmological parameters computes the derivative of E(a) w.r.t. a
    Parameters
    ----------
    a: numpy array
      scale factor normalised to 1 at present time
    output: np.array with a.shape
      dE(a)/da
    """
    return -3.*self.om/(a*a*a*a) - 4.*self.ora/(a*a*a*a*a) + 3.*(self.wa-(1+self.w0+self.wa)/a)*self.ode*self._de_ev(a) - 2.*self.ok/(a*a*a)

  def H_a(self, a):
    """Given a set of cosmological parameters returns H(a)
    Parameters
    ----------
    a: float, list, tuble, numpy array
      scale factor normalised to 1 at present time
    output: np.array with a.shape
      H(a)
    """
    a = mf.convert2array(a)
    return self.H0*np.sqrt(self.E_a(a))

  def _integrand_growth_lin(self,a):
    """Given a set of cosmological parameters returns 1/(a**3*H(a)**3)
    Parameters
    ----------
    a: float, list, tuble, numpy array
      scale factor normalised to 1 at present time.
    output: np.array with a.shape
      a**3*H(a)**3
    """
    a3H3 = a*a*a * np.power(self.E_a(a), 1.5)
    return 1./a3H3
  def _integrand_growth_log(self,t):
    """Given a set of cosmological parameters returns 1/(exp(t*2)*H(exp(t))**3)
    where t = ln(a) in order to do a logarithmic spaced integration
    Parameters
    ----------
    t: float, list, tuble, numpy array
      ln(a), where a is the scale factor normalised to 1 at present time.
    output: np.array with a.shape
      a**2*H(a)**3
    """
    a = np.exp(t)   #get back the scale factor
    a2H3 = a*a * np.power(self.E_a(a), 1.5)
    return 1./a2H3

  def _growth_integral_a(self,a):
    """Given a set of cosmological parameters returns H(a)/(H0*a) int_0^a da'/(a'H(a')/H0)^3
    Parameters
    ----------
    a: float
      scale factor normalised to 1 at present time.
    output: 
      scale factor dependent part of the growth factor definition
    """
    return np.sqrt(self.E_a(a))*spi.quad(self._integrand_growth_log, -10,np.log(a))[0]

  def growth_rat(self,a1,a2):
    """Given a set of cosmological parameters returns the ratio of linear growth factors D(a1)/D(a2)
    where D(a) = 5/2 Omega_m * H(a)/H0 int_0^a da' 1/(a'H(a')/H0)^3

    Parameters
    ----------
    a1,a2: float
      scale factor normalised to 1 at present time.
    output: 
      ratio of linear growth factors D(a1)/D(a2)
    """
    return self._growth_integral_a(a1)/self._growth_integral_a(a2)

  def growth(self,a):
    """Given a set of cosmological parameters returns the linear grows factor 
    D(a) = 5/2 Omega_m * H(a)/H0 int_0^a da' 1/(a'H(a')/H0)^3
    Parameters
    ----------
    a: float
      scale factor normalised to 1 at present time.
    output: np.array with a.shape
      linear growth factor D(a)
    """
    return 2.5 * self.om * self._growth_integral_a(a)

  def f(self, a):
    """Given a set of cosmological parameters returns the logarithmic derivative of
    the linear grows factor D(a) = 5/2 Omega_m * H(a)/H0 int_0^a da' 1/(a'H(a')/H0)^3
    with respect to the scale factor 'a'
    Parameters
    ----------
    a: float
      scale factor normalised to 1 at present time.
    output: np.array with a.shape
      d ln(D(a))/d ln(a)
    """
    D = self.growth(a)    #growth factor
    result = D*self._dE_da(a) + 5.*self.om/(a*a*a)   #D dE(a)/da + 5 Omega_m/a^3

    return result*a / (2.*D*self.E_a(a))

