#!/usr/bin/python
# -*- coding: utf-8 -*-

import contour_func as cf
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
  p.set_usage("""%prog [options] file_root
  Compute the covariance between the different columns of the file.
  """)

  p.add_option("-c", "--chain-file-ext", action="store", dest="cext", default=".0.5_total.txt", help="Extension of the input file containing the chain. [Default: %default]")
  p.add_option("-p", "--parnames-file-ext", action="store", dest="pext", help="Extension of the parameter name file. If not given not used, otherwise the first row of the output file will contain the parameter names in the file.")
  p.add_option("-o", "--outfile-ext", action="store", dest="oext", default=".covmat", help="Extension of the output file. [Default: %default]")

  p.add_option("-n", "--no-weights", action="store_true", dest="nw", default=False, help="If the input file name contains only the parameters. If this option is not given the file is expected do have the first two columns containing the weight and the likelihood or chi square.")

  return p.parse_args()


#############################################
###                 Main                  ###
#############################################

if __name__ == "__main__":   # if is the main

  (opt, args) = options(op.OptionParser(version="%prog version 0.1"))   #create the object optparse
  args = args[0]   #only the first argument to be used

  chain = np.loadtxt(args+opt.cext)
  if(opt.nw == False):  #if the first two columns contain weights and likelihood
    w = chain[:,0]      #save the weights
    chain = chain[:,2:] #save the parameter part of the chain
  else:
    w = None

  fn = open(args+opt.oext, mode="w")   #open the output file

  if(opt.pext != None):   #if the extension of the parameter name file is give, get them
    parnames = cf.get_paramnames([args,], ext=opt.pext)
    for i, p in enumerate(parnames):
      parnames[i] = str(i+1)+") "+ p.strip().split("\t")[0] #extract the parameter name and assign the corresponding number
    fn.write("; ".join(parnames)+"\n")   #if the parameter name file exists, print the parameter names with the numnbers

  #Compute covariance
  cov = np.empty([chain.shape[1], chain.shape[1]+1])   #create a covariance matrix. The first row and column contain simply a counter
  cov[:,0] = np.arange(1, chain.shape[1]+1)   #write the counter in the first column of the matrix

  cov[:,1:] = np.cov(chain, rowvar=0)  

  #save output file
  fn.write("##\t " + "\t\t".join([str(i) for i in range(1,chain.shape[1]+1)]) + "\n")    #write che counter in the first (or second) row of the output file

  fmt = "%02d\t "+"\t".join(["%5.4e",]*chain.shape[1])   #output formatter
  np.savetxt(fn, cov, fmt=fmt)   #write the output file with the nice formatting

  fn.close()    #close the file

  exit()
