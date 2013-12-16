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

def MNRAS_fig():
    """
    Set default matplotlibrc parameters for MNRAS figures
    Parameters
    ----------
    figsize: list
        size of the figure. Default for figure with one panel
    """
    from matplotlib import rcParams

    rcParams['axes.labelsize'] = 10
    rcParams['axes.titlesize'] = 10
    rcParams['font.size'] = 10
    rcParams['legend.fontsize'] = 8
    rcParams['lines.linewidth'] = 1.5
    rcParams['lines.markersize'] = 5

def set_rcParams(**kwargs):
    """
    Set the parameters in matplotlib rcParams. 
    If no option is given, the parameters fall back to the default in 'matplotlibrc' file
    ----------
    set_MNRAS: bool
        set the default parameters from 'pythonlib/myplotmodule.py'.
        This option is called first and others change the parameters

    the other parameters have the same name as options in `matplotlib.rcParams`
    with dots substituted by underscores
    """
    from matplotlib import rcParams

    if kwargs.get("set_MNRAS", False):
        mpm.MNRAS_fig()

    # Figure parameters
    if kwargs.get('figure_figsize') is not None:
        rcParams['figure.figsize'] = kwargs['figure_figsize']

    # lines properties
    if kwargs.get("lines_linewidth") is not None:
        rcParams['lines.linewidth'] = kwargs["lines_linewidth"]

    # font properties
    if kwargs.get("axes_labelsize") is not None:
        rcParams['axes.labelsize'] = kwargs["axes_labelsize"]

    if kwargs.get("font_size") is not None:
        rcParams['font.size'] = kwargs["font_size"]

    if kwargs.get('legend_fontsize') is not None:
        rcParams['legend.fontsize'] = kwargs["legend_fontsize"]

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
    print("Keep default font size of the legend")
  if(kargs.has_key("textcol")==True and kargs["textcol"] == True):
    for i,l in enumerate(lines):
      txt[i].set_color(list(l)[0].get_color())

  pass

#########################################
# functions to operate on dictionaries
# related with plotting on dictionaries 
# of axes/subplots
#########################################
class MyWarning(UserWarning):
    pass

def subsitute_labels(new_labels, labels):
    """
    Substitute some of the labels with new ones
    Parameters
    ----------
    new_labels: list of len(2) lists 
        new_labels[i] = [key, new_label]
    labels: dictionary
        key: short name; value: latex string
    output
    ------
    labels: same as input with new labels
    """
    for k, v in new_labels:
        if k not in labels:
            warn("Key '{}' is not in labels dictionary".format(k), MyWarning)
        else:
            labels[k] = v
    return labels

def change_xylim(axs_dic, x_lim=[], y_lim=[]):
    """
    Change the x and/or y limits of desired axes
    Parameters
    ----------
    axs_dic: dictionary
        key: short param names; value: plt.subplot
    x_lim, y_lim: list of len(3) lists 
        x/y_lim[i] = [key, lower limit, upper limit]
    """
    # x limits
    for k, v1, v2 in x_lim:
        if k not in axs_dic:
            warn("Key '{}' is not in axes dictionary".format(k), MyWarning)
        else:
            axs_dic[k].set_xlim([float(v1), float(v2)])
    # y limits
    for k, v1, v2 in y_lim:
        if k not in axs_dic:
            warn("Key '{}' is not in axes dictionary".format(k), MyWarning)
        else:
            axs_dic[k].set_ylim([float(v1), float(v2)])

def draw_legend(fig, axs_dic, labels, whichax, loc):
    """
    Draw the legend in figure 'fig' or one of the axes in 'axs_dic'.
    If `whichax` is a key in `axs_dic`, the corresponding axes is used to draw
    the legend; if not and a `empty` key exists, it must contain a list of
    empty axes and the first is used for the legend and moved to
    `axs_dic['legend']`; else a figure legend is drawn
    Parameters
    ----------
    fig: matplotlib figure
    axs_dic: dic
        dictionary of subplot objects with keys from key_list
    labels: list of string
        legend labels
    whichax: string
        axis where to draw the legend. If not in the axes dictionary, draw
        figure legend
    loc: matplotlib legend locations
    output
    ------
    matplotlib.legend.Legend instance
    """
    # Check if there are empty axes
    # get legend handles and set where to draw the legend
    if whichax in axs_dic:
        ax = axs_dic[whichax]
        handles,_ = ax.get_legend_handles_labels()
    elif 'empty' in axs_dic:
        ax = axs_dic['empty'][0]
        axs_dic['empty'] = axs_dic['empty'][1:]
        axs_dic['legend'] = ax
        handles,_ = ax.get_legend_handles_labels()
        # make x and y axis invisible.
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
    except KeyError:
        ax = list(axs_dic.values())[0]
        handles,_ = ax.get_legend_handles_labels()
        ax = fig

    return ax.legend(handles, labels, loc=loc)

def plot_horiz_vert(axs_dic, horiz=[], vert=[]):
    """
    Plot horizontal and/or vertical lines in desired axes
    Parameters
    ----------
    axs_dic: dict
        key: short param names; value: plt.subplot
    horiz: list of len(2) lists 
        horiz[i] = [key, y]
    vert: list of len(2) lists 
        vert[i] = [key, x]
    """
    for k, v in horiz:
        if k not in axs_dic:
            warn("Key '{}' is not in axes dictionary".format(k), MyWarning)
        else:
            axs_dic[k].axhline(y=float(v), color='k', ls=':')
    for k, v in vert:
        if k not in axs_dic:
            warn("Key '{}' is not in axes dictionary".format(k), MyWarning)
        else:
            axs_dic[k].axvline(x=float(v), color='k', ls=':')

