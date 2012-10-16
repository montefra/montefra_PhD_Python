# -*- coding: utf-8 -*-
"""This module contains functions that create colors or color maps"""

import numpy as np     #import numpy
import wavelen2rgb as w2rgb   #import the module wavelen2grb: given a wavelenght returns a RGB color

#wavelength limits of the colors in the rainbow: lim_rainb[0,:]:min, lim_rainb[1,:]:max
lim_rainbow = np.array([
                       [475,450],   #blue
                       [620,750],    #red
                       [560,500],   #green
                       [449,380],   #violet
                       [590,619],   #orange
                       [570,589],   #yellow
                       [476,494]   #cyan
		       ])

num_rainbow = len(lim_rainbow)    #number of colors in the rainbow

intensity = 255   #maximum intensity required in the function in w2grb

def cols_shades(n_col, n_lev, alpha=1):
  """given the numbers of colors and the number of levels returns a 
  n_col * n_lev * 4 numpy array of GRBA colors

  Parameters
  ----------
  n_col: integer
    number of colors needed: the list of colors is selected from the variable lim_rainbow.
    If n_col > colors in lim_rainbow, the colors will be used recoursively
  n_lev: integer
    number of shades per color
  alpha: list, tuple or array of integers, floats (optional)
    If 0<=alpha<=1 all the colors will have the correponding alpha value
    If alpha > 1, it will create 'floor(alpha)' linearly spaced alpha values. If alpha < n_col the alpha values will be cicled
    If alpha is a list, tuple or array of float 0<=alpha<=1, they will considered as the the alpha values. If alpha < n_col the alpha values will be cicled
    All the n_lev of each colors will have the same value of alpha

  output: (n_col * n_lev * 4) numpy array
    contains the RGBA quadruplet for the n_lev shades of the n_col colors

  -------
  Example:

  >>> cols_shades(2,2)
      array([[[ 0.38,    0., 0.38,  1.],
              [   0., 0 .28,   1.,  1.]],

             [[   0.,  0.28,  1.,  1.],
              [   0.,  0.75,  1.,  1.]]])
  where each horizontal line is a RGBA quadruplet
  """

  alpha = np.array(alpha)   #convert alpha to a numpy array
  if(alpha.size ==1):   #if only one element 
    if(alpha >= 0 and alpha <= 1):    #in [0,1]
      alpha = alpha.repeat(n_col)     #replicate it n_col times
    elif(alpha > 1):   #if alpha in an integer bigger then 1
      alpha = np.linspace(1,0.3,num=np.floor(alpha), endpoint=False) #create n_col linearly spaced alpha values
  else:
    alpha = alpha.flatten()
    if(np.all(alpha>=0)==False or np.all(alpha<=1)==False):  #check that they are [0,1]
      print "The %s must contain only elements in the range [0,1]" % str(alphatype)[6:-1]
      exit()
  n_alpha = alpha.size

  cols = np.zeros([n_col, n_lev, 4])   #create a vector to contain all the colors and their shade

  for i in range(n_col):
    cols[i,:,0] = np.linspace(lim_rainbow[i%num_rainbow, 0], lim_rainbow[i%num_rainbow, 1], num=n_lev, endpoint=False)   #fill each row with the required shade of the color
    for j in range(n_lev):
      cols[i,j,:3] = w2rgb.wavelen2rgb(cols[i,j,0], MaxIntensity=intensity)   #convert the wavelenght in a RGB triplet
    cols[i,:,3] = alpha[i%n_alpha]

  cols[:,:,:3] /= intensity

  return np.around(cols, decimals=2)    #return RGBA color in the range [0-1]

def custom_colors():
  """My custom colors for contourplots"""

  custcol = [
             [[ 0., 0.75, 1., 1.],          #blues
              [ 0., 0.53, 1., 1.]],
             [[ 1., 0.47, 0., 0.85],          #reds
              [ 1., 0., 0., 0.85]],
	     [[ 90./255., 255./255., 0./255., 0.75],       #greens
	      [ 0./255., 210./255., 0./255., 0.75]]
	     ]
	      #[ 0., 0.39, 0., 0.7]]]

  return np.array(custcol)

