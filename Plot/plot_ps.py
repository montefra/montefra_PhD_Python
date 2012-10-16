#!/usr/bin/python

import numpy as np    # import numpy
import matplotlib.pyplot as plt   #import matplotlib
import matplotlib.axes as plta   #import matplotlib axes
import myplotmodule as mpm

if __name__ == "__main__":   # if is the main

  kd = mpm.readfile("/Users/Ciccio/Data/b6_matterpower_z0_ext.dat")  #read power spectrum and correlation function files
  cf = mpm.readfile("/Users/Ciccio/Data/test_correlation_z_0_ariel_lin.dat")

  xs=30./mpm.inc2cm   # window size
  ys=xs/2
  plt.figure(1, figsize=(xs,ys))

  plt.subplot(1,2,1)
  plt.loglog(cf[:,0], cf[:,1],'b', lw=3)
  plt.axis([1,150, 1e-4, 5])
  plt.xlabel(r"$\mathbf{r\, [Mpc/h]}$", fontsize=20)
  plt.ylabel(r"$\mathbf{\xi(r)}$", fontsize=20)

  plt.subplot(1,2,2)
  plt.loglog(kd[:,0], kd[:,1],'b', lw=3)
  plt.axis([0.001,5, 1, 5e+4])
  plt.xlabel(r"$\mathbf{k\, [h/Mpc]}$", fontsize=20)
  plt.ylabel(r"$\mathbf{P(k)\, [(Mpc/h)^3]}$", fontsize=20)

  plt.subplots_adjust(left=0.09, right=0.98, bottom=0.14, top=0.98, wspace=0.25)

  #fig=plt.figure()
  #fig1=fig.add_subplot(1,1,1)
  #p=fig1.plot(kd[:,0], kd[:,1])
  #fig1.set_xscale('log')

  plt.savefig("/Users/Ciccio/Dottorato/Slides/Student Seminar II/cf_ps.eps")

