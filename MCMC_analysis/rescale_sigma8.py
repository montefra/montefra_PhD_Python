#!/usr/bin/python
# -*- coding: utf-8 -*-

import cosmology as c
import numpy as np
import optparse as op

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

  p.set_usage("""usage: %prog [options] input_file_name output_file_name
  Read the chain stored in 'input_file_name', rescale sigma_8 from redshift 
  z[1] to z[0] (or a[1] to a[0]) using the cosmological parameters 
  in the chain and save the new chain in 'output_file_name'.
  """)

  p.add_option("-z", "--redshift", action="store", dest="z", type=float, nargs=2, default=[0,0.312782], help="Rescale sigma8 between redshift z[1] and z[0]. [Default: %default]")
  p.add_option("-a", "--scale-factor", action="store", dest="a", type=float, nargs=2, help="Rescale sigma8 between scale factor a[1] and a[0]. If given, the redshifts are not considered.")
 
  p.add_option("-c", "--columns", action="store", dest="cols", type=float, nargs=6, default=[20,6,8,9,24,21], help="Columns containing the parameters: Omega_matter, Omega_k, w_DE, w_a, H_0, sigma_8. Mind the order of the parameters. [Default: %default]")

  return p.parse_args()

if __name__ == "__main__":
  
  (opt, args) = options(op.OptionParser(version="%prog version 1."))

  if( opt.a == None):
    opt.a = 1./(1.+np.array(opt.z))
  else:
    opt.a = np.array(opt.a)

  chain = np.loadtxt(args[0])   #read the input file

  for i,(om,ok,w0,wa,H0,s8) in enumerate(chain[:,opt.cols]):
    chain[i,opt.cols[-1]] = s8* c.Cosmology(om=om,ok=ok,w0=w0,wa=wa,H0=H0).growth_rat(opt.a[0], opt.a[1])

  np.savetxt(args[1], chain, fmt='%9.8e', delimiter='\t')

  exit()
