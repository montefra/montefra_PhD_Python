#!/usr/bin/python

import numpy as np    # import numpy
import matplotlib.pyplot as plt   #import matplotlib.pyplot
import scipy.interpolate as interpolate  # interpolation routines under scipy
import sys    #mondule sys
import myplotmodule as mpm   #my module

if __name__ == "__main__":   # if is the main
  """takes the linear, p1loop ps with Ried 2010 cosmology,
  the measured power spectrum, its covariance matrix and the windows function for SDSS DR7
  and allow to do 'fit by hand'"""
  
  if(len(sys.argv) < 4):
    print "b, kstar and A_MC needed as first 3 command line arguments"
    exit()

  kpl = mpm.readfile("/Users/Ciccio/Data/reid_bf_matterpower.dat")  #reads the linear power spectrum
  kp1 = mpm.readfile("/Users/Ciccio/Data/reid_bf_1loop.dat")  #reads the 1 loop power spectrum
  kpm = mpm.readfile("/Users/Ciccio/Data/LRGdr7/lrgDR7_ccmeasurements_kmax02kmin02_maxLv2.ALL_MAGCOVv3.txt", cols=[0,3])  # measured ps
  kpmvar = mpm.readfile("/Users/Ciccio/Data/LRGdr7/lrgDR7_cov_kmax02kmin02_maxLv2.ALL_MAGCOVv3.txt")  #read the covariance matrix
  kpmvar = np.sqrt(np.diag(kpmvar)) #extract the 
  winfunc = mpm.readfile("/Users/Ciccio/Data/LRGdr7/lrgDR7_windows_kmax02kmin02_maxLv2.ALL_MAGCOVv3.txt")  #read the window function
  kwin = mpm.readfile("/Users/Ciccio/Data/LRGdr7/lrgDR7_kbands_kmax02kmin02_maxLv2.ALL_MAGCOVv3.txt")  #read the k values assocated with the long side of the window function

  # spline the linear and 1 loop ps to kwin values 
  spline = interpolate.splrep(kpl[::,0], kpl[::,1], s=0)  #spline the linear ps
  kpli= interpolate.splev(kwin, spline)
  spline = interpolate.splrep(kp1[::,0], kp1[::,1], s=0)  #spline the 1loop ps
  kp1i= interpolate.splev(kwin, spline)
  
  model = float(sys.argv[1])**2 * ( np.exp(-(kwin/float(sys.argv[2]))**2)*kpli + float(sys.argv[3])*kp1i )  #create the model power spectrum
  model = np.dot(winfunc, model)

  xs=20./mpm.inc2cm   # window size
  ys=xs
  fig=plt.figure(1, figsize=(xs,ys))
  plt.subplots_adjust(left=0.12, right=0.98, bottom=0.12, top=0.98)

  coord=[0.02, 0.21, 2e+3, 4e+4]  #x and y range
  pm = plt.loglog(kpm[:,0], kpm[:,1], 'rD', lw=3)  #measured ps
  pmerp = plt.loglog(kpm[:,0], kpm[:,1]+kpmvar, 'r-.', lw=2)  #and its variance
  pmerm = plt.loglog(kpm[:,0], kpm[:,1]-kpmvar, 'r-.', lw=2)
  pmod = plt.loglog(kpm[:,0], model, 'k', lw=3)
  #plin = plt.loglog(kdl[:,0], kdl[:,1], 'b', lw=3)  # linear ps
  #pdm = plt.loglog(kdDM[:,0], kdDM[:,1], 'rD', lw=3)  # DM ps
  plt.axis(coord)  # axix range
  plt.xlabel(r"$\mathbf{k\, [h/Mpc]}$", fontsize=20)  #axis names
  plt.ylabel(r"$\mathbf{P(k)\, [(Mpc/h)^3]}$", fontsize=20)
  plt.legend([pm, pmerp, pmod],['measured','variance','model'], loc= 'upper right')  #legend

  plt.savefig("/Users/Ciccio/Dottorato/Slides/Student Seminar II/sdss_dr7_lrg.eps")


