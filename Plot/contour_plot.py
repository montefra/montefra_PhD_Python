#!/usr/bin/python
# -*- coding: utf-8 -*-

import colmaps as cm
import contour_func as cf
import itertools as it  #optimized iteration tools
import matplotlib.pyplot as plt  #plotting stuff
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

  p.set_usage("""%prog [options] file_root1 [file_root2 ... file_rootn] [legtag1 [legtag2 .. [legtagn]]]
  This program reads n MCMC chains and plot 2D contour plots for two of the parameter
  WARNING: the file name and structrure is based on COSMOMC
  Each line of the input chain must be 'counter, likelihood, parameters'
  The arguments of the program are intended as the root of the different chains to plot. 
  At least one must be provided
  The choise of the two parameters to plot can be given with the column numbers 
  Otherwise the program list the parameter names as read from 'file_root1.params'
  and read from stdin the 2 numbers given by the user.
  If the legend option is selected, a number of strings equal to the number of file_root are espected 
  after the file_roots and will be used as legend tags.

  WARNING: the number of colums in the input chains and in the various 'file_rooti.params'
  must be consistent""")  # program explanation

  p= mpm.options(p)  #send to the function to create the standard options

  #add new options
  p.add_option("-c", "--columns", action="store", dest="cols", type="int", nargs=2, help="Allows to choose the two column used for the plot. WARNING: the input chain has the following structrure: chain[:,0]=counter, chain[:,1]=likelihood, chain[:,2:]=parameters. If this option is not called the number of the columns to plot will be asked.")   #x range values from the user

  p.add_option("-b", "--num_bins", action="store", dest="n_bins", type="int", default="80", help="Number of bins for the 2D histogram. Assumed square [default: %default]")    #number of bins to create the 2D histogram for the contourplots

  p.add_option("-l", "--level", action="append", dest="level", type="float", help="Confidence level for the contour plot. 'LEVEL' must be [0,100]. To give more that one level call this keyword more than once.")    #number of bins to create the 2D histogram for the contourplots

  p.add_option("-s", "--smooth", action="store_true", default="False", dest="smooth", help="Turn on the smothing of the 2D histogram [default: %default]")

  p.add_option("-f", "--fraction", action="store", dest="fraction", default="0.5", help="The first 'FRACTION' of the chain that has been discarted: it enters in the total chain file name. 'FRACTION' must be [0,1]. [default: %default]")
  
  p.add_option("-e", "--extension", action="store", type="string", dest="ext", default="txt", help="Extension of the input chains. [Default: %default]")

  p.add_option("-a", "--alpha", action="store_true", dest="alpha", default=False, help="If this option is selected, transparency turned on.")
  
  p.add_option("-w", "--line-width", action="store", dest="lw", type=float, default=2, help="Set the line width for the contour plots. [Default: %default]")
  p.add_option("--font-size", action="store", type=float, dest="sfont", default=20, help="Font size of the axis. [Default: %default]")
  p.add_option("--legend-fsize", action="store", type=float, dest="lsfont", help="Font size of the legend tags. If not given set as the axis font size")

  p.add_option("--legend", action="store", dest="legend", type=int, help="Turn on the legend on the given plot, using as tags the second half of the command line arguments. The option requires an integer value from 0 to 10 giving the position of the legend. 0='best' (http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.legend).")
  p.add_option("--color-text", action="store_true", dest="textcol", default=False, help="If selected, writes the text of the legend with the same color of the lines of the contour plots.")
  p.add_option("--invert-legend", action="store_true", dest="invleg", default=False, help="Invert the legend entries")
  p.add_option("--legend-frame", action="store", dest="legfr", help="Set the legend frame to the given matplotlib color. Default no frame")

  p.add_option("--figure-size", action="store", type=float, nargs=2, dest="fsize", default=[10.,10.], help="x and y size of the plot in cm. [Default: %default]")
  p.add_option("--bounds", action="store", dest="bounds", type=float, nargs=4, default=[0.20, 0.20, 0.78, 0.78], help="axes rect [left, bottom, width, height]. [Default: %default]")

  p.add_option("--horizontal", action="store", type=float, dest="hline", help="Draw a horizontal line at y='hline'")
  p.add_option("--vertical", action="store", type=float, dest="vline", help="Draw a vertical line at x='vline'")
  p.add_option("--diagonal", action="store", type=float, nargs=2, dest="dline", help="Draw a straight line with y='dline[0]'*x + 'dline[1]'")
  p.add_option("--line", action="store", type=float, nargs=3, dest="line", help="Draw a line defined as y='line[0]'*x^'line[1]' + 'line[2]'")

  p.add_option("-t", "--text", action="append", nargs=3, dest="text", help="Writes the text 'text[2]' at coordinates x='text[0]', y='text[1]'. To have more than one text element call the option more than once")
  p.add_option("--text-fsize", action="store", type=float, dest="tsfont", help="Font size of the texts. If not given set as the axis font size")
  p.add_option("--text-bkgr", action="store", dest="textfr", help="Set the text background to any matplotlib color. Default no background")

  p.add_option("--xlabel", action="store", dest="xlabel", help="Allows custom x label")
  p.add_option("--ylabel", action="store", dest="ylabel", help="Allows custom y label")

  return p.parse_args()


#############################################
###                 Main                  ###
#############################################

if __name__ == "__main__":   # if is the main

  (options, args) = options(op.OptionParser(version="%prog version 1.1"))   #create the object optparse

  #check if the input options and parameters are ok   
  if(len(args) < 1):   #file_root required
    print """At least one argument is required by the program
    Type './%s -h' for more informations""" % os.path.basename(sys.argv[0])
    exit()
  else:
    if(options.legend != None):   #if the legend required, check if the number of args is even, and then save the second half of the args
      if(len(args)%2 == 1):
        print "If the legend is required, I expect the second half of the argument list to be the legend tags, but I have an odd number of args."
	sys.exit(2)
      else:
        leg = args[len(args)/2:]  #save the legend tags
	args = args[:len(args)/2] #save the file roots

  if(options.level != None):  #check if it's fine
    for i in options.level: 
      if(i < 0 or i > 100):
	print "The confidence level must be in the interval [0,100]"
	exit ()
    options.level = np.sort(np.array(options.level))/100.  #get the confidence levels as %
    n_levs = options.level.size   #number of levels

  t0,t1 = mpm.linestyles[:2]   #invert the first and second elements of the line style
  mpm.linestyles[:2] = t1,t0
  #t0,t1 = mpm.linestyles[2:4]   #invert the third and fourth elements of the line style
  #mpm.linestyles[2:4] = t1,t0

  paramnames = cf.get_paramnames(args, verbose=options.verbose)   # get parameter names

  if(options.cols == None):    #if the option is not used ask which columns one wants
    for i,pn in enumerate(paramnames):
      print "%i) %s" %(i+2, pn.strip().split("\t")[0])
    temp = raw_input("Please give me the number corresponding to the two parameters you want to plot separated by one or more space (if more the two are given the others will be discarded):\n")
    options.cols = map(lambda x: int(x), temp.split()[0:2])   #convert the first 2 elementes of raw input into integers incremented by two
  else:    #otherwise get the option values
    options.cols = list(options.cols)
  for i in options.cols:
    if(i<2 or i > len(paramnames)+2):
      print "Please give a correct column number"
      exit()
  options.cols.insert(0,0)

  for i in range(len(args)):
    if options.fraction != "None":
      args[i] += "."+str(options.fraction)
  totchain = cf.get_chains(args, options.cols, verbose=options.verbose, ext="_total."+options.ext)   #get the chains and the index where they are singolarly stored

  if(options.verbose == True):
    print "Computing x and y range"
  xr, yr = [np.infty,-np.infty], [np.infty,-np.infty]
  for ch in totchain:
    xr = [min(xr[0], np.amin(ch[:,1])), max(xr[1], np.amax(ch[:,1]))]
    yr = [min(yr[0], np.amin(ch[:,2])), max(yr[1], np.amax(ch[:,2]))]

  if(options.verbose == True):
    print "Converting the MCMC chains in 2d histograms"

  h2D, xe, ye = cf.hist2D(totchain, xr=xr, yr=yr, bins=options.n_bins, smooth=options.smooth) #convert the chains to 2D histograms

  if(options.level != None):  #if the confidence levels are given
    amplitude = cf.h2D2conflev(h2D, options.level) # return an array containing the ampliudes relative to the maxima corresponding to the confidence levels
    n_shades = 1   #if alpha shades wanted give the number of them
    if(options.alpha == True):
      n_shades = len(args)
    colors = cm.cols_shades(len(args), n_levs, alpha=n_shades)    #create authomatic colors
    #colors = cm.custom_colors()
    
  options.fsize = tuple(np.array(options.fsize)/mpm.inc2cm)
  fig=plt.figure(1, figsize=options.fsize)  #open the figure
  subpl = fig.add_axes(options.bounds)
  #fig=plt.figure(1, figsize=(10./mpm.inc2cm, 10./mpm.inc2cm))  #open the figure
  #fig.subplots_adjust(left=0.19, right=0.93, bottom=0.15, top=0.95) # set subplot 
  #subpl = fig.add_subplot(1,1,1)  #create a subplots
  for i,h,x,y in it.izip(it.count(), h2D, xe, ye):   #plot the filled contours plots
    if(options.level != None):  #if the confidence levels are given
      col = tuple([tuple(colors[i,j,:]) for j in range(n_levs)])
      alphas = colors[i,0,3]
      subpl.contourf(x, y, h, colors=col, alpha=alphas, levels=amplitude[i], origin='lower')  #plot the contour plot
    else:
      subpl.contourf(x, y, h, origin='lower')  #contour plot

  lines = []   #inizialise the lines and texts for the legend
  for i,h,x,y in it.izip(it.count(), h2D, xe, ye):   #plot the filled contours plots
    if(options.level != None):  #if the confidence levels are given
      lines.append(subpl.contour(x, y, h, colors=(tuple(colors[i,-1,:3]),), levels=amplitude[i], linestyles=mpm.linestyles[i%mpm.numls], linewidths=options.lw, origin='lower').collections[0])  #plot the contours and get the outermost line
    else:
      lines.append(subpl.contour(x, y, h, linestyles=mpm.linestyles[i%mpm.numls], linewidths=options.lw, origin='lower').collections[0])  #plot the contours and get the outermost line

  if(options.hline != None):
    subpl.plot(x, options.hline*np.ones_like(x), ":k", lw=1.5)
  if(options.vline != None):
    subpl.plot(options.vline*np.ones_like(y), y, ":k", lw=1.5)
  if(options.dline != None):
    subpl.plot(x, options.dline[0]*x+options.dline[1], ":k", lw=1.5)
  if(options.line != None):
    subpl.plot(x, options.line[0]*np.power(x,options.line[1])+options.line[2], ":k", lw=1.5)

  if(options.text != None):
    if(options.tsfont == None):
      options.tsfont = options.sfont
    for x,y,t in options.text:
      txt = subpl.text(float(x),float(y),t, fontsize=options.tsfont)
    if(options.textfr != None):
      txt.set_backgroundcolor(options.textfr)

  if(options.xr != None):
    subpl.set_xlim(options.xr)  #set the axis size
  if(options.yr != None):
    subpl.set_ylim(options.yr)

  if(options.xlabel == None):
    options.xlabel = "$"+paramnames[options.cols[1]-2].strip().split("\t")[1]+"$"
  if(options.ylabel == None):
    options.ylabel = "$"+paramnames[options.cols[2]-2].strip().split("\t")[1]+"$"
  plt.xlabel(options.xlabel, fontsize=options.sfont)  #axis names
  plt.ylabel(options.ylabel, fontsize=options.sfont)
  for label in subpl.get_xticklabels():
    label.set_fontsize(options.sfont)
  for label in subpl.get_yticklabels():
    label.set_fontsize(options.sfont)

  if(options.legend != None):  #write the legend if required
    if(options.invleg):
      lines.reverse()
      leg.reverse()
    l = subpl.legend(lines, leg, options.legend, numpoints=1, borderpad=0.3, columnspacing=0.8, handlelength=1.4, borderaxespad=0.1, labelspacing=0.2)
    if(options.legfr == None):
      l.draw_frame(False)   #delete the frame
    else:
      l.get_frame().set_edgecolor(options.legfr)
      l.get_frame().set_facecolor(options.legfr)

    txt = l.get_texts()
    if(options.lsfont == None):
      options.lsfont = options.sfont
    plt.setp(txt, fontsize=options.lsfont)
    if(options.textcol == True):
      for i,line in enumerate(lines):
        txt[i].set_color(list(line.get_color()[0]))

  if(options.outfile == None):
    plt.show()
  else:
    plt.savefig(options.outfile)

  exit()
