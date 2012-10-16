#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np    # import numpy
import matplotlib.pyplot as plt   #import matplotlib
import scipy.interpolate as interpolate  # interpolation routines under scipy
import optparse as op   #option parse analysis
import os.path as osp   #operating sistem stuff
import myplotmodule as mpm

if __name__ == "__main__":   # if is the main

  fdir = "/home/montefra/data1/CosmoMC/Test_integration"  #directory of operation

  #set the options using the funtion in mymodulplot
  option = op.OptionParser()   #create the object optparse
  option = mpm.options(option)  #send to the function to create the standard options
  #customize the standard options
  option.set_defaults(outfile=fdir+"/cov_test.eps")

  (options, args) = option.parse_args()  #get the option names and the arguments

  n_itr=["800", "1000", "1200"] #, "1200"] #["500", "800", "1000", "1200"] #, "3000"]  #number of steps
  n_ita=["300"] #, "3000"]  #number of steps
  mpm.linestyles.extend(mpm.symbols)  #set line stiles
  n_sym = mpm.numls + mpm.numsy  #number of lines and symbols

  #myp = mpm.readfile(fdir+"/1loop.dat") #read my calculation
  k = mpm.readfile(fdir+"/1loop.3000.3000.dat", cols=[0]) #read the values of k of the interpolation
  myp = mpm.readfile(fdir+"/1loop.3000.3000.dat", cols=[1]) #read the values of k of the interpolation

  #spline = interpolate.splrep(myp[::,0], myp[::,1], s=0)  #spline the linear ps
  #myp = interpolate.splev(k, spline)  #get the intorpolated ps

  xs=20./mpm.inc2cm   # window size
  ys=xs
  fig=plt.figure(1, figsize=(xs,ys))
  subp = fig.add_subplot(1, 1, 1)  #add a subplot
  plt.subplots_adjust(left=0.12, right=0.98, bottom=0.12, top=0.98)
  coord=[0.00, 0.52, 0.99, 1.01]  #x and y range
  lines=[] #it will contain a list of lines
  bin=[]   #it will contain the labels for the legend
  l=0  #counter

  for i in n_itr:
    for j in n_ita:
      pfile = fdir+"/1loop."+i+"."+j+".dat"
      if(osp.isfile(pfile) == False):
	continue   #if the file doesn't exists skip it
      pk = mpm.readfile(pfile, cols=[1], verbose=False) #read the values of k of the interpolation

      lines.append(subp.plot(k, pk/myp, mpm.linestyles[l%n_sym], color=mpm.colors[l%mpm.numcol], lw=2))  #create the line with the model
      bin.append(i+"."+j)
      l=l+1   #increment the counter
  
  #subp.loglog(k,myp)

  subp.plot(coord[0:2], [0.99, 0.99], 'k--', lw=2)
  subp.plot(coord[0:2], [1.01, 1.01], 'k--', lw=2)
  subp.plot(coord[0:2], [1., 1.], 'k-', lw=2)
  subp.axis(coord)  # axix range
  plt.xlabel(r"$\mathbf{k\, [h/Mpc]}$", fontsize=20)  #axis names
  plt.ylabel(r"$\mathbf{P(k)\, [(Mpc/h)^3]}$", fontsize=20)
  subp.legend(lines, bin,loc= 'upper left')  #legend
  plt.draw()
  plt.savefig(options.outfile) #output file
#  plt.show()


  exit()