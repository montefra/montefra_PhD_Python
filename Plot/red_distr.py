#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pytplot as plt
import glob

"""Compute the radial density of galaxies in the lasdamas mocks and in the dr7 LRG
"""

def input():
  """Contains the input file name of the mocks and of the SDSS DR7 LRG"""
  fmocks = "~/data1/SDSSDR7/LasDamas/"
  fdr7 = "~/data1/SDSSDR7/DR7/"
  return fdr7, glob.glob(fmocks)

if __name__ == "__main__":   # if is the main

  fdr7, fmocks = input()

  rdzdr7 = np.loadtxt(fdr7, usecols=[])




  exit()
