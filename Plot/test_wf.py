#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np    # import numpy
import matplotlib.pyplot as plt   #import matplotlib.pyplot
import scipy.interpolate as interpolate  # interpolation routines under scipy
import sys    #mondule sys
import os    #operating system depending things
import myplotmodule as mpm   #my module

if __name__ == "__main__":   # if is the main
  """Test: see how strong is the impact on the model PS 
  selecting different dimensions of the of Reid09"""

  if(os.name == "posix"):
    rootc = "/afs/ipp-garching.mpg.de/home/m/montefra/cambout/"
    rootd = "/afs/ipp-garching.mpg.de/home/m/montefra/Cosmomc/cosmomc.01.2010/data/"
  elif(os.name == "mac"):
    rootc = "/Users/Ciccio/Data/"
    rootd = "/Users/Ciccio/Data/LRGdr7/"
  else:
    print "WARNING: this script could not work on you operating system"
    exit()

  kpl = mpm.readfile(rootc+"reid_halofit_matterpower.dat")  #reads the linear power spectrum
  kp1 = mpm.readfile(rootc+"reid_bf_1loop.dat")  #reads the 1 loop power spectrum
  winfunc = mpm.readfile(rootd+"lrgDR7_windows_kmax02kmin02_maxLv2.ALL_MAGCOVv3.txt")  #read the window function
  kwin = mpm.readfile(rootd+"lrgDR7_kbands_kmax02kmin02_maxLv2.ALL_MAGCOVv3.txt")  #read the k values assocated with the long side of the window function
  km = mpm.readfile(rootd+"lrgDR7_ccmeasurements_kmax02kmin02_maxLv2.ALL_MAGCOVv3.txt", cols=[0])  # k of the measured ps

  # spline the linear and 1 loop ps to kwin values 
  spline = interpolate.splrep(kpl[::,0], kpl[::,1], s=0)  #spline the linear ps
  kpli= interpolate.splev(kwin, spline)
  spline = interpolate.splrep(kp1[::,0], kp1[::,1], s=0)  #spline the 1loop ps
  kp1i= interpolate.splev(kwin, spline)

  b = pow(1.15,2.)  #model parameters
  kstar = 0.2
  amc = 0.6

  #model = b* (np.exp(-(kwin/kstar)**2)*kpli + amc*kp1i )  #create the model power spectrum
  model = kpli   #use only the linear or HALOFIT power spectrum

  xs=20./mpm.inc2cm   # window size
  ys=xs
  fig=plt.figure(1, figsize=(xs,ys))
  plt.subplots_adjust(left=0.12, right=0.98, bottom=0.12, top=0.98)

  coord=[0.02, 0.21, 0.97, 1.02]  #x and y range

  kmax = np.arange(10,40,2)/100.  #max k used to evaluate the window function
  lines=[] #it will contain a list of lines
  j=0  #counter
  mpm.linestyles.extend(mpm.symbols)
  wrefmodel = np.dot(winfunc[:,:], model[:])  #model convolved with the full window function

  for i in kmax:    # loop over different values of kmax for the window function
    num = kwin[kwin < i].size #number of k in the win function smaller than i
    wmodel = np.dot(winfunc[:,:num], model[:num])  #model convolved with the window function
    lines.append(plt.plot(km, wmodel/wrefmodel, mpm.linestyles[j%(mpm.numls+mpm.numsy)], color=mpm.colors[j%mpm.numcol], lw=3))  #create the line with the model
    j=j+1   #increment the counter

  plt.plot(coord[0:2], [0.99, 0.99], 'k--', lw=2)
  plt.plot(coord[0:2], [1.01, 1.01], 'k--', lw=2)
  plt.plot(coord[0:2], [1., 1.], 'k-', lw=2)
  plt.plot([0.15,0.15], coord[2:4], 'k:', lw=2)
  plt.axis(coord)  # axix range
  plt.xlabel(r"$\mathbf{k\, [h/Mpc]}$", fontsize=20)  #axis names
  plt.ylabel(r"$\mathbf{P(k)\, [(Mpc/h)^3]}$", fontsize=20)
  plt.legend(lines, kmax,loc= 'lower left')  #legend

  plt.savefig("test_halofit_wf.eps")


