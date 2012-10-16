#!/usr/bin/python
# -*- coding: utf-8 -*-
"""This modul contains functions and variables that I need often for plotting"""

import matplotlib.collections as mplcol
import matplotlib.pyplot as plt
import numpy as np    # import numpy
import optparse as op   #option parse analysis

inc2cm=2.54    # 1 inch = 2.54 cm

linestyles = ['-', '--', '-.', ':' ]   #linestyle list
numls = len(linestyles)  #number of linestyles
symbols = ['.',',','o','v','^','<','>','1','2','3','4','s','p','*','h','H','+','x','D','d','|','_']  #symbol list
numsy = len(symbols)  #number of symbols
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k'] #, 'w']  #easy color array
numcol = len(colors)  #number of colors


def options(option):
  """this function set some default parse option that can be used by lot of programs
  """
  option.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Produces more output.")  # verbose option
  option.add_option("-o", "--output-file", action="store", dest="outfile", type="string", help="Output file name [default: %default]")   #output file name

  option.add_option("-x", "--xrange", action="store", dest="xr", type="float", nargs=2, help="Custom x range value: two values must be given. If the option is not used, the x range will be automatically determined")   #x range values from the user

  option.add_option("-y", "--yrange", action="store", dest="yr", type="float", nargs=2, help="Custom y range value: two values must be given. If the option is not used, the y range will be automatically determined")   #y range values from the user

  return option


def linecollection(ax, lines, cmap="prism", lw=1, ls="solid", legend=False):
  """draws a collection of lines in axes 'ax' using a given color map
  Parameters
  ----------
  ax: matplotlib.axes.AxesSubplot or matplotlib.axes.Axes object
    axes reference for the plot
  lines: list
    list of lines to draw
  cmap: string (optional)
    color map name
  lw: float (optional)
    line width
  ls: string (optional)
    line style [ ‘solid’ | ‘dashed’ | ‘dashdot’ | ‘dotted’ ]
  legend: bool (optional)
    return items for the legend

  output: matplotlib.axes.AxesSubplot or matplotlib.axes.Axes object
    axes reference with line collection added
  """
  lc = mplcol.LineCollection(lines, cmap=plt.get_cmap(cmap, lut=len(lines)), linewidths=lw, linestyles=ls)  #create a collection of lines in order to plot them with different colors from a colormap
  lc.set_array(np.arange(len(lines)))   #add an array of length len(lines) in order to color the lines in different ways
  ax.add_collection(lc)   #plot the collection
  
  if(legend == True):
    return ax, lc
  else:
    return ax

def plot_winmat_rows(ax, wm, k, ext=None, spacing=1, cmap="prism", lw=1, ls="solid", legend=False):
  """Given a window matrix 'wm' and the array of 'k' to which the elements in the rows of 'wm' correspond
  make a plot with line collections of the rows of 'wm' on the axes 'ax'
  
  Parameters
  ----------
  ax: matplotlib.axes.AxesSubplot or matplotlib.axes.Axes object
    axes reference for the plot
  wm: NxM array
    window matrix
  k: 1xM array
    values of k of the rows of the window matrix
  ext: list (optional)
    minimum and maximum index of the columns to use
  spacing: integer (optional)
    plot only one every 'spacing' rows
  cmap: string (optional)
    color map name
  lw: float (optional)
    line width
  ls: string (optional)
    line style [ ‘solid’ | ‘dashed’ | ‘dashdot’ | ‘dotted’ ]
  legend: bool (optional)
    return items for the legend

  output: matplotlib.axes.AxesSubplot or matplotlib.axes.Axes object
    axes reference with line collection added
  """
  if(ext==None):
    ext = [0,wm.shape[0]]   #if the extrema not given plot all of them
  lines = []
  for i in np.arange(ext[0], ext[1], spacing):   #create a list of lines
    lines.append(zip(k, wm[i,:]))

  return linecollection(ax, lines, cmap=cmap, lw=lw, ls=ls, legend=legend)

def legend(ax, lines, tags, pos, **kargs):
  """Place the legend in position 'pos' in the 'ax' reference frame

  Parameters
  ----------
  ax: class 'matplotlib.axes.Axes' or class 'matplotlib.axes.AxesSubplot'
    axes or subplot where to draw the legend
  lines: list
    lines of the legend
  tags: list
    tags for the legend
  pos: integer [0,10] or corresponding string
    position of the legend in the plot
  **kargs: dictionary (optional)
    dictionary of options
  """
  l = ax.legend(lines, tags, loc=pos, numpoints=1, borderpad=0.3, columnspacing=0.8, handlelength=1.4, borderaxespad=0.1, labelspacing=0.2)
  l.draw_frame(False)   #delete the frame
  txt = l.get_texts()
  if(kargs.has_key("lsfont")==True and kargs["lsfont"]!=None):
    plt.setp(txt, fontsize=kargs["lsfont"])
  else:
    print "Keep default font size of the legend"
  if(kargs.has_key("textcol")==True and kargs["textcol"] == True):
    for i,l in enumerate(lines):
      txt[i].set_color(list(l)[0].get_color())

  pass
