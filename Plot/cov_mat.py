#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np    # import numpy
import matplotlib.pyplot as plt   #import matplotlib.pyplot
import sys    #mondule sys
import myplotmodule as mpm   #my module

def cor(cov):
  """Given the covariance matrix returns the correlation matrix: 
  C_ij/sqrt(C_ii*C_jj) """
  cor = np.empty_like(cov)  #create the correlation matrix
  for i in range(cov.shape[0]):
    for j in range(cov.shape[1]):
      cor[i,j] = cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])  #fill the correlation matrix
  return cor # return the correlation matrix

if __name__ == "__main__":   # if is the main

  bin = [11, 13, 33]  #bin names
  files = ["/home/montefra/data1/BASICC/MeasuredPS_mpi_24_z0/Rau",]*len(bin)  #create the file names
  for i in range(len(bin)):
    files[i] = files[i] + "/PS_fftw_TCS_" + str(bin[i]) + "a_24_lin/ps_fftw.z0.cov_mat_mpi24_lin_TCS_" + str(bin[i]) + "a.dat"

  xs = 20./mpm.inc2cm   # window size
  ys = xs/3.
  fig = plt.figure(1, figsize=(xs,ys))  # open window figure
  fig.subplots_adjust(left=0.04, right=0.98, bottom=0.12, top=0.98, wspace=0.25) # set window size
  plt.gray()  #grayscale


  n_rows = 1  # number of rows and colums
  n_cols = 3

  for i in range(n_rows*n_cols):
    subp = fig.add_subplot(n_rows, n_cols, i+1)  #add a subplot

    cov = mpm.readfile(files[i], cols=[3,4])   #read the covariance matrix
    dim = np.sqrt(cov.shape[0])   #reshape it and cut at a given k
    k = cov[:dim,0]
    k = k[k<0.5]
    cov = np.reshape(cov[:,1], [dim,dim])
    dim = k.shape[0]
    cov=cov[:dim,:dim]
    extent=[k[0],k[-1],k[0],k[-1]]

    subp.imshow(cor(cov).T, extent=extent, interpolation='nearest',origin='lower') 
    #if(i == 2):
      #plt.colorbar()

  plt.draw()
  plt.savefig("covmat_z0_11_13_33.eps")
  plt.show()

  exit()

