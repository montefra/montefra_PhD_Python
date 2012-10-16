#!/usr/bin/python

import numpy as np    # import numpy
import matplotlib.pyplot as plt   #import matplotlib.pyplot
import scipy.interpolate as interpolate  # interpolation routines under scipy
import myplotmodule as mpm   #my module

if __name__ == "__main__":   # if is the main
  """plot the dark matter vs linear PS
  and the splined one"""

  kdl = mpm.readfile("/Users/Ciccio/Codes/Data/b6_matterpower_z0_ext.dat")  #reads the linear power spectrum 
  kdnw = mpm.readfile("/Users/Ciccio/Codes/Data/b6_matterpower_z0_nowig_ext.dat")  #reads the linear power spectrum 
  kdDM = mpm.readfile("/Users/Ciccio/Codes/Data/ps_fftw.z0.mpi24_lin_meanTCS_DM.dat")  #read DM power spectrum 

  sp1 = interpolate.splrep(kdnw[::,0], kdnw[::,1], s=0)  #spline the non wiggle ps
  kdnwi = interpolate.splev(kdDM[:,0], sp1)
  stride = 5
  offset = 1

  xs=30./mpm.inc2cm   # window size
  ys=xs/2
  fig=plt.figure(1, figsize=(xs,ys))
  plt.subplots_adjust(left=0.09, right=0.98, bottom=0.14, top=0.98, wspace=0.25)

  plt.subplot(1,2,1) # first plot
  coord=[0.004, 0.5, 1e+2, 5e+4]  #x and y range
  plin = plt.loglog(kdl[:,0], kdl[:,1], 'b', lw=3)  # linear ps
  pdm = plt.loglog(kdDM[:,0], kdDM[:,1], 'rD', lw=3)  # DM ps
  perp = plt.loglog(kdDM[:,0], kdDM[:,1]+kdDM[:,2], 'r-.', lw=2)  #DM variance
  perm = plt.loglog(kdDM[:,0], kdDM[:,1]-kdDM[:,2], 'r-.', lw=2)
  plt.axis(coord)  # axix range
  plt.xlabel(r"$\mathbf{k\, [h/Mpc]}$", fontsize=20)  #axis names
  plt.ylabel(r"$\mathbf{P(k)\, [(Mpc/h)^3]}$", fontsize=20)
  plt.legend([plin, pdm, perp],['linear','DM','variance'], loc= 'lower left')  #legend

  plt.subplot(1,2,2)  # second plot
  coord = [0., 0.3, 0.7, 1.3]   #x and y range
  plt.plot(coord[0:2], [1, 1], 'k--', lw=2)
  plt.plot(kdl[:,0], kdl[:,1]/kdnw[:,1], 'b', lw=3) #lin
  plt.plot(kdDM[:,0], kdDM[:,1]/kdnwi[:],  'rD', lw=3)  #mean DM
  plt.plot(kdDM[:,0], (kdDM[:,1]+kdDM[:,2])/kdnwi[:], 'r-.', lw=2)
  plt.plot(kdDM[:,0], (kdDM[:,1]-kdDM[:,2])/kdnwi[:], 'r-.', lw=2)
  plt.axis(coord)
  plt.xlabel(r"$\mathbf{k\, [h/Mpc]}$", fontsize=20)
  plt.ylabel(r"$\mathbf{P(k)\, }$", fontsize=20)


  #fig=plt.figure()
  #fig1=fig.add_subplot(1,1,1)
  #p=fig1.plot(kd[:,0], kd[:,1])
  #fig1.set_xscale('log')

  plt.savefig("/Users/Ciccio/Dottorato/Slides/Student Seminar II/ps_full_ratio.jpg")

