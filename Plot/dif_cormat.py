#!/usr/bin/python
# -*- coding: utf-8 -*-

"""This program computes the correlation matrix from a list of files,
extract the variance and make a plot of the variances and of all the 
differences of the correlation matrices.
"""
import numpy as np  #numpy
import matplotlib.pyplot as plt  #plotting stuff
import optparse as op  #import optsparse: allows nice command line option handling
import sys    #mondule sys
import os    #contain OS dependent stuffs: it helps with sistem portability  
import myplotmodule as mpm   # my model used to plot
import itertools

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
  p.set_usage("""%prog [options] fname_cov1 fname_cov2 [fname_cov3 ... fname_covn]
  Given at n file names containg the covariance matrix,
  it computes the correlation matrix, extract the variance, make a plot
  with all the variances and n(n-1)/2 differences of the correlation matrices.
  The name of at least two files must be given.
  The program is expecting a file with the structure: i, j, ki, kj, covariance
  The covariance matrices need to be evaluated at the same wave number
  """)

  p.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Produces more output.")  # verbose option

  p.add_option("-o", "--output-file", action="store", dest="outfile", type="string", help="Output file name root.")   #output file name

  p.add_option("-k", "--kmax", action="store", dest="kmax", type="float", nargs=1, help="Maximum k of the range of covariance matrix: two values must be given. If the option is not used, the range will be automatically determined")   #x range values from the user

  return p.parse_args()

def file2corr(fname, kmax=None):
  """Given a file name read the columns with k and covariance,
  cut the covariance according to the range, if given, and returns 
  the the krange, the variance and correlation matrix

  Parameters:
  fname: string
    file name
  range: two elemnts list (optional)
    range of the covariance matrix

  Output: sequence of 3 arrays
    k, variance, correlation matrix
  """

  k, cov = np.loadtxt(fname, usecols=[3,4]).T   #read the covariance
  dim_cov = np.sqrt(cov.size)  #get the dimension of the covariance assuming that is square

  cov = np.reshape(cov, (dim_cov, dim_cov))   #make the covariace squared
  k = k[:dim_cov]  #get the kvalues
  if(kmax!=None):   #if the maximum value of k is given cut k and cov according to it
    k = k[k<kmax]
    cov = cov[:k.size, :k.size]

  var = np.sqrt(np.diag(cov))
  for i,j in itertools.product(range(k.size), repeat=2):   #loops over the two dimensions of the covariance
    cov[i,j] /= var[i]*var[j]   #compute the correlation matrix Cij/sqrt(Cii*Cjj)

  return k, var, cov  #return

def get_diff(flist, begin="pk", end=".dat"):
  """Given a list of string extract the part between "begin" and "end" and return it

  Paramters:
  ---------
  flist: list
    list of strings
  begin: string (optional)
    string to be found before the different part
  end: string (optional)
    string to be found after the different part

  output: list
    list of strings containing the difference
  """
  output=[]
  for i in flist:
    output.append(i[i.find(begin)+len(begin):i.find(end)])   #extract the part between begin and end
  return output


if __name__ == "__main__":   # if is the main

  (options, args) = options(op.OptionParser(version="%prog version 0.1"))   #create the object optparse

  if(len(args)<2):
    print """At least two file names must be given.
    Type './%s -h' for more informations""" % os.path.basename(sys.argv[0])
    exit()

  n_args = len(args)   #number of input file names

  var=[]    #create list of variances
  cor=[]    #create list of correlation matrices
  vmin=1e+6   #min and max of the variance
  vmax=0.
  for i in args:
    if(options.verbose):
      print "reading file %s, extracting the variance and computing the correlation matrix" %i
    k, tvar, tcor = file2corr(i, kmax=options.kmax)   #give file name and max k and get back k, variance and correlation matrix
    tvx, tvn = np.amax(tvar), np.amin(tvar)   #evaluate the min and max of the variance
    if(tvx>vmax):
      vmax=tvx
    if(tvn<vmin):
      vmin=tvn
    var.append(tvar)
    cor.append(tcor)

  tag = get_diff(args)   #get the tag to be used in the labels

  if(options.verbose):
    print "Plotting the variances together"

  mpm.linestyles.extend(mpm.symbols)   #merge lines and symbols

  #plot the variances
  xs = ys = 9./mpm.inc2cm   # window size
  fig1 = plt.figure(1, figsize=(xs,ys))  # open window figure
  fig1.subplots_adjust(left=0.25, right=0.95, bottom=0.21, top=0.95, wspace=0.2, hspace=0.25) # set window size
  subf1 = fig1.add_subplot(1,1,1, xscale="log", yscale="log")  #add a subplot

  lines=[] #it will contain a list of lines
  for i,v in enumerate(var):
    lines.append(subf1.plot(k, v, mpm.linestyles[i%(mpm.numls+mpm.numsy)], color=mpm.colors[i%mpm.numcol], lw=3))  #create the line

  subf1.set_xlabel(r"$\mathbf{k\, [h/Mpc]}$", fontsize=15)  #axis names
  subf1.set_ylabel(r"$\mathbf{variance\, [(Mpc/h])^3}$", fontsize=15)
  subf1.axis([np.amin(k)*0.9,np.amax(k)*1.1, vmin*0.9,vmax*1.1])  # axix range
  l = subf1.legend(lines, tag,loc= 'upper right', )  #legend
  l.draw_frame(False)

  if(options.outfile != None):   #save the file if required
    fig1.savefig(options.outfile+"var.eps") #output file

  if(options.verbose):
    print "Plotting the covariance differences together"
  #plot the differences between the correlation matices
  n_spl = n_args*(n_args-1)/2  #number of differences to be computed
  n_cols = 3
  if(n_spl < n_cols):     #maximum number of columns = 3
    n_cols = n_spl
  n_rows = np.ceil(n_spl/n_cols)
  dcor=[]  #store the differences between correlation matrices
  cmin=1e+6
  cmax=0.
  for i,j in itertools.combinations(range(n_args),2):   #compute the differences of the
    temp = cor[i] - cor[j]
    tvx, tvn = np.amax(temp), np.amin(temp)   #evaluate the min and max of the difference of the correlation
    if(tvx>cmax):
      cmax=tvx
    if(tvn<cmin):
      cmin=tvn
    dcor.append(temp)

  xs = ys = 9./mpm.inc2cm   # window size
  if(n_cols > 1):
    xs *= 2.
    ys *= n_row*0.66
  fig2 = plt.figure(2, figsize=(xs,ys))  # open window figure
  fig2.subplots_adjust(left=0.25, right=0.95, bottom=0.21, top=0.95, wspace=0.2, hspace=0.25) # set window size
  extent=[k[0],k[-1],k[0],k[-1]]

  for i,c in enumerate(dcor):
    subf2 = fig2.add_subplot(n_cols, n_rows, i+1)  #add a subplot
    subf2.imshow(c.T, extent=extent, interpolation='nearest', origin='lower', vmin=cmin, vmax=cmax, cmap='gray')
    if(i%n_cols == 0):
      subf1.set_ylabel(r"$\mathbf{k\, [h/Mpc]}$", fontsize=15)  #axis names
    if(n_spl-i < n_cols):
      subf1.set_xlabel(r"$\mathbf{k\, [h/Mpc]}$", fontsize=15)  #axis names

  if(options.outfile != None):
    fig1.savefig(options.outfile+"diffcor.eps") #output file

  if(options.outfile == None):
    plt.show()

  exit()

