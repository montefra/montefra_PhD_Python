#!/usr/bin/python
# -*- coding: utf-8 -*-

import contour_func as cf
import itertools as it
import matplotlib.pyplot as plt  #plotting stuff
import matplotlib.ticker as tic  #axis formatter
import my_statistic as ms
import myplotmodule as mpm   # my model used to plot
import numpy as np  #numpy
import optparse as op  #import optsparse: allows nice command line option handling
import os    #contain OS dependent stuffs: it helps with sistem portability  
import sys    #mondule sys

def options(p):
  """
  This function accept the option parser istance and fill it with options.

  The version must be set when creating the optparse istance

  Parameters
  ----------
  p: option parser istance

  output: tuple with the options and the arguments
  ---------
  """

  p.set_usage("""%prog [options] min max file_root1 [file_root2 .. [file_rootn]] [legtag1 [legtag2 .. [legtagn]]]
  This program plot the values of the given parameters as function of the maximum scale
  used to perform the fit. 
  File names list is created as 'file_rooti+numpy.arange(min,max+1,stride)/100+"."+extension'
  'pw' and 'extension' of the parameter and chain files can be set with the options.
  If more than a file root given, the files must have the same structure, and the constraints for the
  chosen parameters will be overplotted.
  If the legend option is selected, a number of strings equal to the number of file_root are espected 
  after the file_roots and will be used as legend tags.
  """) #set usage

  p.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="Produces more output.")  # verbose option
  p.add_option("-o", "--output-file", action="store", dest="outfile", type="string", help="Output file name [default: %default]")   #output file name

  p.add_option("-c", "--columns", action="append", dest="cols", type=int,
      help="""Allows to choose the columns used for the plot. WARNING: the input
      chain has the following structrure: chain[:,0]=counter,
      chain[:,1]=likelihood, chain[:,2:]=parameters. If this option is not
      called all the parameters will be plotted.""")   #number of parameters to plot
  p.add_option("--stride", action="store", dest="stride", type=int,
      default="1", help="""Set the stride of the values of k_max used in the
      computation of the chain. [Default: %default]""")

  p.add_option("-n", "--nsigma", action="store", dest="ns", type=int,
      default="1", help="Plot the ns-sigma error bars. [default: %default]")

  p.add_option("-m", "--marker-size", action="store", type=float,
      dest="smarker", default="4", help="Marker size. [Default: %default]")
  p.add_option("--font-size", action="store", type=int, dest="fsize",
      default="15", help="Axis font size. [Default: %default]")
  p.add_option("-w", "--line-width", action="store", type=float, dest="lwidth",
      default=1, help="Line width. [Default: %default]")

  p.add_option("--shift", action="store", type=float, dest="shift", default=0,
      help="""If errorbars drawn, for it allows to displace the x axis values by
      the given factor if more than one fileroot is given. [Default: %default]""")
  p.add_option("-f", "--fill-between", action="store", dest="fill", type=float,
      nargs=2, help="""If two floats between [0,1] given, instead of teh
      errorbars, the standard deviation is plotted with fill_between with alpha
      in the given range. No alpha for eps.""")

  p.add_option("-s", "--skip", action="store_true", dest="skip", default=False,
      help="Skip non existing files")

  p.add_option("--paramlist", action="store_true", dest="plist", default=False,
      help="Read the parameter list, print it out and close the program")

  p.add_option("--ext-chain", action="store", type="string", dest="ec",
      default="0.5_total.dat", help="""Extension of finale part of the input file
      containing the chain. [Default: %default]""")
  p.add_option("--ext-parmn", action="store", type="string", dest="ep",
      default="paramnames", help="""Extension of finale part of the input file
      containing the parameters' names. [Default: %default]""")

  p.add_option("--horizontal", action="append", type=float, nargs=2,
      dest="horiz", help="Plot an horizontal line at y=horiz[1] in the horiz[0] plot.")

  p.add_option("-l", "--legend", action="store", dest="legend", type=int,
      help="""Turn on the legend on the 'legend' plot, using as tags the second
      half of the command line arguments.""")
  p.add_option("--figlegend", action="store", dest="fleg", help="""Turn on the
      figure legend in position 'fleg'. It overides the standard legend""")

  p.add_option("-b", "--bounds", action="store", dest="bounds", type=float,
      nargs=4, default=[0.20, 0.20, 0.97, 0.97], help="""Left, bottom, right and
      top limits of the plotting area within the figure window. [Default: %default]""")

  p.add_option("-r", "--rescale", action="append", dest="rescale", type=float,
      nargs=2, help="""Rescale the plot rescale[0] by multiplying the values by
      rescale[1]. The same rescaling factor is shown in the y label.""")

  p.add_option("--y-label", action="append", dest="ylab", nargs=2, help="""Set
      custom y axis label 'ylab[1]' in subplot number 'ylab[0]'""")

  p.add_option("--y-range", action="append", dest="yrange", nargs=3,
      help="""Set the yrange to ('yrange[1]', 'yrange[2]') in plot 'yrange[0]'""")

  return p.parse_args()

def skip(flist, slist=False):
  """
  Skip non existing files on a given list

  Parameters
  ----------
  flist: list
    list of file names to be checked
  slist: bool (optional)
    if given, returns the list with the index of discarded files
    
  ----------
  output: list [, list]
    list with non existing files excluded [, optional: list of indeces of the non existing files]
  """
  discard=[]  #inizialize the list of discarded files
  for i,fn in enumerate(ifiles):   #check if the file exists
    if(os.path.isfile(fn) == False):    #remove the non existing ones
      discard.append(i)
  if len(discard) > 0:
    for disc in discard:
      flist.pop(disc)
    
  if(slist == True):
    return flist, discard
  else:
    return flist


if __name__ == "__main__":   # if is the main

  (options, args) = options(op.OptionParser(version="%prog version 1.0"))   #create the object optparse

  k_max = np.arange(int(args[0]),int(args[1])+1,options.stride)/100.  #create array with the k_max

  if(options.legend != None or options.fleg != None):
    n_args = len(args)-2   #number of arguments excuded the firt two
    if(n_args%2 == 1):   #if there is an odd number of args there is some tag missing or exceding
      print """Something doesn't fit: if you want the legend you better check
      that the number of file roots is equal to the number of legend tags"""
      sys.exit(2)
    leg = args[n_args/2+2:]  #save the legend tags
    args = args[:n_args/2+2] #redefine args with only min, max and the file_root

  files = [args[2],]*k_max.size  # parameter name files
  for i,km in enumerate(k_max):   #create the file names
    files[i] += "%3.2f.%s" %(km,options.ep)
  if(options.skip == True):  #check if to skip files
    files = skip(files)
  paramnames = cf.get_paramnames(files, ext="", verbose=options.verbose)   # and get parameter names

  if(options.plist == True):    #print out parameters name and number
    for i,pn in enumerate(paramnames):
      print i+2, pn.strip().split("\t")[0]
    sys.exit(1)

  if(options.cols == None):   #if no columns given all the columns read
    options.cols = list(np.arange(len(paramnames))+2)
  else:
    discard = []
    for i in options.cols:   #check that the columns are correct
      if(i<2 or i>len(paramnames)+1):
        print "The column %d is not acceptable or do not exists. I'll skip it" %i
	discard.append(i)
    for i in discard:   #discard the non existing columns
      options.cols.remove(i)
    if(len(options.cols) == 0):
      print "All the column number are too big or to small"
      sys.exit(2)
  options.cols.insert(0,0)  #read also the first column
  n_plots = len(options.cols)-1   #number of plots

  horiz = [[None, 0.] for i in range(n_plots)] #initialise the list of lists for the horizontal lines
  if(options.horiz != None):   #if a tuple of values to draw the horizontal line is given order it
    for h in options.horiz:
      horiz[int(h[0])] = list(h)
  rescale = [1 for i in range(n_plots)] #initialise the list of lists for the rescaling of the plot y axes
  if(options.rescale != None):   #if a tuple of values to draw the horizontal line is given order it
    for h in options.rescale:
      rescale[int(h[0])] = h[1]
  ylab = [[None, None] for i in range(n_plots)] #initialise the list of lists for the custom y labels
  if(options.ylab != None):   #if a tuple of values to draw the horizontal line is given order it
    for y in options.ylab:
      ylab[int(y[0])] = list(y)
  yrange = [[None, None, None] for i in range(n_plots) ] #initialise the list of lists for the custom y limits
  if(options.yrange != None):
    for y in options.yrange:
      yrange[int(y[0])] = list(y)


  #create the numpy arrays that will contain the mean and the stddev of the columns
  mean = np.empty([len(files), n_plots])
  stddev = np.empty([len(files), n_plots])

  #initialise the figure
  xs=10./mpm.inc2cm   # window size
  ys= 2.5*(n_plots+1)/mpm.inc2cm
  fig=plt.figure(1, figsize=(xs,ys))  #open the figure
  rheight = (options.bounds[3]-options.bounds[1])/n_plots     #height of each subplot

  if(options.fill != None):
    alphas = np.linspace(options.fill[0], options.fill[1], len(args[2:]))    #alpha valued to be used in the fill_between
  else:
    alphas = [1 for i in args[2:]]
  tk_max = np.copy(k_max)  #make a copy of the k_max to make the plots
  lines=[] #initialise the line list
  for i,fr,col,ls,m,a in it.izip(it.count(), args[2:], it.cycle(mpm.colors),
      it.cycle(mpm.linestyles), it.cycle(mpm.symbols[2:-2]), alphas):   #loop through the file roots

    files = [fr,]*k_max.size  # chain file names
    for j,km in enumerate(k_max):   #create the file names
      files[j] += "%3.2f.%s" %(km,options.ec)
    if(options.skip == True):  #check if to skip files
      files, discard = skip(files, slist=True)
      tk_max = np.delete(k_max, discard)   #discard values of kmax if some of the file is not there

    #read each chain and obtain mean and std
    for j, fn in enumerate(files):
      if(options.verbose == True):
	print "Reading file %s" % fn
      chain = np.loadtxt(fn, usecols=options.cols)  #read
      if(options.verbose == True):
	print "Computing the mean of the parameters in file %s" % fn
      mean[j,:] = np.average(chain[:,1:], axis=0, weights=chain[:,0])  #mean of the columns containing parameters
      if(options.verbose == True):
	print "Computing the standard deviation of the parameters in file %s\n" % fn
      stddev[j,:] = ms.stddev(chain[:,1:], weights=chain[:,0], axis=0)  #stddev of all the columns containing parameters

    for j,h,r,yl,yr,c in it.izip(it.count(), horiz, rescale, ylab, yrange, options.cols[1:]):   #loop through the various columns to plot
      box = [options.bounds[0], options.bounds[1]+(n_plots-j-1)*rheight,
	  options.bounds[2]-options.bounds[0], rheight]   #create the plot box
      subpl = fig.add_axes(box)  #create a subplots

      if(np.absolute(r-1) > 1e-3):   #if a rescale required
        mean[:,j] *= r
        stddev[:,j] *= r
      if(options.fill == None):
	if(j==0):
	  lines.append(subpl.errorbar(k_max+i*options.shift, mean[:,j],
	    yerr=stddev[:,j], c=col, ls=ls, lw=options.lwidth, marker=m,
	    ms=options.smarker)[0])
	else:
	  subpl.errorbar(k_max+i*options.shift, mean[:,j], yerr=stddev[:,j],
	      c=col, ls=ls, lw=options.lwidth, marker=m, ms=options.smarker)
      else:
	if(j==0):
	  lines.append(subpl.plot(k_max+i*options.shift, mean[:,j], c=col,
	    ls=ls, lw=options.lwidth, marker=m, ms=options.smarker))
	else:
	  subpl.plot(k_max+i*options.shift, mean[:,j], c=col, ls=ls,
	      lw=options.lwidth, marker=m, ms=options.smarker)
	if(options.outfile != None and os.path.splitext(options.outfile)[1] == '.eps'):
	  subpl.fill_between(k_max+i*options.shift, mean[:,j]+stddev[:,j],
	      mean[:,j]-stddev[:,j], color=[140./255., 130./255., 255./255.,
		1.], zorder=1)   #this is when eps is generated
	else:
	  subpl.fill_between(k_max+i*options.shift, mean[:,j]+stddev[:,j],
	      mean[:,j]-stddev[:,j], color=col, alpha=a, edgecolor="w")
      if(h[0] != None):
	subpl.plot([k_max[0]*0.95,(k_max[-1]+i*options.shift)*1.02], [h[1],h[1]], 'k--')
  
      if(yl[0] == None):
	label = "$"+paramnames[c-2].strip().split()[1]+"$"
	if(np.absolute(r-1) > 1e-3):   #if a rescale required
	  if(r.is_integer()==True):
	    strr = str(int(r))
	  else:
	    strr = str(r)
	  label = "$"+strr+"$"+label
      else:
	label = yl[1]
      subpl.set_ylabel(label, fontsize=options.fsize)

      if( yr[0] != None ):  #set the custom yrange
	subpl.set_ylim( bottom=float(yr[1]), top=float(yr[2]) )
      #exclude the first and the last visible tick labels
      ymin,ymax = subpl.get_ylim()
      ytl = subpl.get_ymajorticklabels()
      tl = subpl.yaxis.get_majorticklocs()
      ytl[(tl<ymin).sum()].set_visible(False)
      ytl[-(tl>ymax).sum()-1].set_visible(False)
      if(len(ytl)-2>5):
        for iy in range(2,len(ytl)-1,2):
	  ytl[iy].set_visible(False)
      for iy in range(len(ytl)):   #set tick label font size
	ytl[iy].set_fontsize(options.fsize)
      if(j < n_plots-1):
	subpl.set_xticklabels([])
      else:
	subpl.set_xlabel("$k_{\mathrm{max}}\,[h/Mpc]$", fontsize=options.fsize)
	xtl = subpl.get_xmajorticklabels()  #get xlabels and change font size
	for ix in range(len(xtl)):   #set tick label font size
	  xtl[ix].set_fontsize(options.fsize)
      for label in subpl.get_xticklabels():
	label.set_fontsize(options.fsize)
      for label in subpl.get_yticklabels():
	label.set_fontsize(options.fsize)

      subpl.set_xlim([k_max[0]*0.95, (k_max[-1]+i*options.shift)*1.02])  #set the axis size

  if(options.legend != None and options.fleg == None):   #draw the legend in the desired axes
    box = [options.bounds[0],
	options.bounds[1]+(n_plots-options.legend-1)*rheight,
	options.bounds[2]-options.bounds[0], rheight]   #create the plot box
    subpl = fig.add_axes(box)  #create a subplots
    l = subpl.legend(lines, leg, loc="best", numpoints=1, borderpad=0.3,
	columnspacing=0.8, handlelength=1.4, borderaxespad=0.1,
	labelspacing=0.2)
    l.draw_frame(False)
  if( options.fleg != None):  #draw figure legend
    l = fig.legend(lines, leg, loc=options.fleg, numpoints=1, borderpad=0.3,
	columnspacing=0.8, handlelength=1.4, borderaxespad=0.1,
	labelspacing=0.2, frameon=False)

  if(options.outfile == None):
    plt.show()
  else:
    plt.savefig(options.outfile)

  sys.exit()
