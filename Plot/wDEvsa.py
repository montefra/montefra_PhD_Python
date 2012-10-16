#!/usr/bin/python
# -*- coding: utf-8 -*-

import colmaps as cm
import itertools as it  #optimized iteration tools
import matplotlib.pyplot as plt  #plotting stuff
import myplotmodule as mpm   # my model used to plot
import numpy as np  #numpy
import optparse as op  #import optsparse: allows nice command line option handling

default_wwa=(-1.,0.)   #default values for wde and wa
default_swwa=(0.15,0.5)   #default values for wde and wa
default_cov=[0.,]
default_sigma=[1.,]

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

  p.set_usage("""%prog [options]
  Plot the mean and error areas as from the cosmomc for w_{DE} as function of redshift or scale factor a.
  In order to overplot multiple lines call '-b', '-s' and '-c' multiple times (their length are checked)
  """)

  p.add_option('-b', '--best-fit', dest="bf", action="append", nargs=2, type=float, help="Best fit values for w_0 and w_a. [Default: "+", ".join([str(x) for x in default_wwa])+"]")
  p.add_option('-s', '--stddev', dest="stddev", action="append", nargs=2, type=float, help="Standard deviation for w_0 and w_a. [Default: "+", ".join([str(x) for x in default_swwa])+"]")
  p.add_option('-c', '--covariance', dest='cov', action='append', nargs=1, type=float, help="Covariance between w_0 and w_a. [Default: "+str(default_cov[0])+"]")

  p.add_option('--sigma', dest='sigma', action='append', nargs=1, type=float, help="Values of sigma to plot. Call this option more than ones to have multiple sigma areas. [Default: "+str(default_sigma[0])+"]")

  p.add_option('-a', "--scale-factor", dest='use_a', action="store_true", default=False, help="Uses the scale factor 'a' instead of the redshift 'z' as x axes")

  p.add_option('-x', '--xrange', dest='xr', action='store', nargs=2, type=float, default=[0.,2.], help="Set x range. w_DE(z) or W_DE(a) will be computed on this range.")
  p.add_option('-y', '--yrange', dest='yr', action='store', nargs=2, type=float, help="Set custom y range")

  p.add_option("--legend", action="store", dest="legend", type=int, help="Turn on the legend in the plot. If given a number of arguments equal to the number of lines to plot must be given. The option requires an integer value from 0 to 10 giving the position of the legend. 0='best' (http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.legend).")

  p.add_option("-f", "--font-size", action="store", type=float, dest="fsize", default=20, help="Font size of the axis. [Default: %default]")
  p.add_option("--legend-fsize", action="store", type=float, dest="lsfont", default=13, help="Font size of the legend tags. [Default: %default]")

  p.add_option("-w", "--line-width", action="store", dest="lw", type=float, default=2, help="Set the line width for the contour plots. [Default: %default]")

  p.add_option("--bounds", action="store", dest="bounds", type=float, nargs=4, default=[0.20, 0.15, 0.96, 0.94], help="[left, bottom, right, top] in units of the figure canvas. [Default: %default]")

  p.add_option("--horizontal", action="append", type=float, dest="hline", help="Draw a horizontal line at y='hline'. If more than one needed call this option more times")

  p.add_option("-o", "--output-file", action="store", dest="outfile", type="string", help="Output file name [default: %default]")   #output file name

  return p.parse_args()

def check_opt_args(opt, args):
  """Set the default values of some of the options

  Parameters
  ----------
  opt: dict
    options

  output: dict
    options with some of the None substituted by default values
  """

  if(opt.bf == None):    #set default values for w_de and w_a
    opt.bf = [default_wwa]
  if(opt.stddev == None):    #set default standard deviations for w_de and w_a
    opt.stddev = [default_swwa]
  if(opt.cov == None):    #set default correlation between w_de and w_a
    opt.cov = [tuple(default_cov)]
  if(len(opt.bf)!=len(opt.stddev) or len(opt.bf)!=len(opt.cov)):
    raise SystemExit("The number of elements in the best fit values, standard deviations and covariances is not the same")

  if(opt.sigma == None):
    opt.sigma = default_sigma
  opt.sigma.sort()
  opt.sigma.reverse()

  if(opt.legend != None and len(args) != len(opt.bf)):
    raise SystemExit("If the legend desired, the number of tags, given through command line arguments, must be the same of the number of w_DE to plot")

  return opt, args

if __name__ == "__main__":   # if is the main

  (opt, args) = options(op.OptionParser(version="%prog version 0.1"))   #create the object optparse
  opt, args = check_opt_args(opt, args)   #set the defaults

  t0,t1 = mpm.linestyles[:2]   #invert the first and second elements of the line style
  mpm.linestyles[:2] = t1,t0

  z = np.linspace(opt.xr[0], opt.xr[1], num=100)    #redshifts
  a = 1./(1.+z)    #scale factors

  mean = lambda x,w0,wa: w0+(1.-x)*wa   #local function with w_DE(a)
  stddev = lambda x,s0,sa,c: np.sqrt(s0*s0 + 2.*(1.-x)*c + (1+x*x-2.*x)*sa*sa)   #local function with the associated standard deviation
  wde, sdwde = [],[]   #initialise the lists containing the values of the best fit w_DE and the corresponding standard deviation as function of redshift
  for i, bf, sd, c in it.izip(it.count(), opt.bf, opt.stddev, opt.cov):
    wde.append( mean(a, bf[0], bf[1]) )  #w_DE(a)
    sdwde.append( stddev(a, sd[0], sd[1], c) )   #variance of w_DE as function of 'a'
    ap = 1.+c/(sd[1]*sd[1])   #pivot scale factor
    print r"{0:d}) $a_{{\mathrm{{p}}}} = {1:5.4f}$, $z_{{\mathrm{{p}}}} = {2:5.4f}$ and $w_{{\mathrm{{DE}}}}(a_{{mathrm{{p}}}}) = {3:5.4f}\pm{4:5.4f}$".format(i+1, ap, 1./ap-1., mean(ap, bf[0], bf[1]), stddev(ap, sd[0], sd[1], c))

  if(opt.use_a == True):   #use a as x axis if required
    z = a
    xlabel="a"   #label for the x axes
  else:
    xlabel="z"

  #get colors according to my color map
  #colors = cm.cols_shades(len(wde), len(opt.sigma)+1, alpha=len(wde))
  colors = cm.custom_colors()

  fig = plt.figure(1, figsize=(10./mpm.inc2cm, 10./mpm.inc2cm))
  fig.subplots_adjust(left=opt.bounds[0], bottom=opt.bounds[1], right=opt.bounds[2], top=opt.bounds[3])
  spl = fig.add_subplot(111)

  lines = []   #lines for the legend
  for i,w,sw in it.izip(it.count(), wde, sdwde):
    sp = opt.sigma[-1]
    spl.fill_between(z, w+sp*sw, y2=w-sp*sw, facecolor=colors[i,-2,:3], alpha=colors[i,-2,3], edgecolor='none')
    spl.plot(z, w+sp*sw, linewidth=1, linestyle=mpm.linestyles[i%mpm.numls], color=colors[i,-1,:3])
    spl.plot(z, w-sp*sw, linewidth=1, linestyle=mpm.linestyles[i%mpm.numls], color=colors[i,-1,:3])
    if(len(opt.sigma)>1):
      for j,sp,sm in it.izip(it.count(), opt.sigma[:-1], opt.sigma[1:]):
	spl.fill_between(z, w+sp*sw, y2=w+sm*sw, facecolor=colors[i,j,:3], alpha=colors[i,j,3], edgecolor='none')
	spl.fill_between(z, w-sp*sw, y2=w-sm*sw, facecolor=colors[i,j,:3], alpha=colors[i,j,3], edgecolor='none')
	spl.plot(z, w+sp*sw, linewidth=1, linestyle=mpm.linestyles[i%mpm.numls], color=colors[i,-1,:3])
	spl.plot(z, w-sp*sw, linewidth=1, linestyle=mpm.linestyles[i%mpm.numls], color=colors[i,-1,:3])
  for i,w,sw in it.izip(it.count(), wde, sdwde):
    lines.append(spl.plot(z, w, color=colors[i,-1,:3], linewidth=opt.lw, linestyle=mpm.linestyles[i%mpm.numls]))

  if(opt.hline != None):   #if require draw horizontal lines
    for hl in opt.hline:
      spl.plot(z, hl*np.ones_like(z), ":k", lw=1.5)


  spl.set_xlabel(xlabel, fontsize=opt.fsize)    #set x and y label
  spl.set_ylabel(r"w$_{\mathrm{DE}}$("+xlabel+")", fontsize=opt.fsize)
  for label in spl.get_xticklabels():
    label.set_fontsize(opt.fsize)
  for label in spl.get_yticklabels():
    label.set_fontsize(opt.fsize)
  if(opt.yr != None):
    spl.set_ylim(opt.yr)
  
  if(opt.legend != None):
    mpm.legend(spl, lines, args, opt.legend, **opt.__dict__)

  if(opt.outfile == None):
    plt.show()
  else:
    plt.savefig(opt.outfile)

  exit()
