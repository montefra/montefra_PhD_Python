#!/usr/bin/python
# -*- coding: utf-8 -*-

import contour_func as cf   #functions used in contour plot
import itertools as it  #optimized iteration tools
import matplotlib.pyplot as plt  #plotting stuff
import myplotmodule as mpm   # my model used to plot
import numpy as np  #numpy
import optparse as op  #import optsparse: allows nice command line option handling

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
  This program reads n MCMC chains and plot 1D histograms for given set of the parameter
  WARNING: the file name and structrure is based on COSMOMC
  Each line of the input chain must be 'counter, likelihood, parameters'
  The arguments of the program are intended as the root of the different chains to plot. 
  At least one must be provided
  If the legend option is selected, a number of strings equal to the number of file_root are espected 
  after the file_roots and will be used as legend tags.

  WARNING: the number of colums in the input chains and in the various 'file_rooti.params'
  must be consistent""")  # program explanation

  p.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Produces more output.")  # verbose option
  p.add_option("-o", "--output-file", action="store", dest="outfile", type="string", help="Output file name [default: %default]")   #output file name

  #add new options
  p.add_option("-c", "--columns", action="append", dest="cols", type="int", help="Column numbers used for the plot. Call this options more than one to read more columns. The input chain has the following structrure: chain[:,0]=counter, chain[:,1]=likelihood, chain[:,2:]=parameters.")   #x range values from the user

  p.add_option("-b", "--num_bins", action="store", dest="n_bins", type="int", default="40", help="Number of bins for the histogram. [default: %default]")    #number of bins to create the 2D histogram for the contourplots

  p.add_option("--chain-ext", action="store", dest="cext", default=".0.5_total.txt", help="'file_rooti+cext': file containing the chains. [Default: %default]")
  p.add_option("--parname-ext", action="store", dest="pext", default=".paramnames", help="'file_rooti+cext': file containing the name of the parameters. [Default: %default]")

  p.add_option("-w", "--line-width", action="store", dest="lw", type=float, default=2, help="Set the line width for the contour plots. [Default: %default]")
  p.add_option("--font-size", action="store", type=float, dest="sfont", default=20, help="Font size of the axis. [Default: %default]")
  p.add_option("--legend-fsize", action="store", type=float, dest="lsfont", help="Font size of the legend tags. If not given set as the axis font size")
  p.add_option("--ticklabel-fsize", action="store", type=float, dest="tsfont", help="Font size of the tick labels. If not given set as the axis font size")

  p.add_option("--legend", action="store", dest="legend", type=int, nargs=2, help="Turn on the legend on the legend[0] plot, using as tags the second half of the command line arguments. legend[1] in the range [0,10] gives the position of the legend. 0='best' (http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.legend).")
  p.add_option("--color-text", action="store_true", dest="textcol", default=False, help="If selected, writes the text of the legend with the same color of the lines of the contour plots.")

  p.add_option("--figure-size", action="store", type=float, nargs=2, dest="fsize", default=[10.,10.], help="x and y size of the plot in cm. [Default: %default]")
  p.add_option("--bounds", action="store", dest="bounds", type=float, nargs=6, default=[0.20, 0.20, 0.98, 0.98, 0, 0.2], help="[left, bottom, right, top, horizonal space, vertical space] in units of the figure canvas. [Default: %default]")
  p.add_option("-r", "--rows-number", action="store", dest="nrows", type=int, default=1, help="Number of rows to be plot. The number of columns is computed according to the number of columns to read. [Default: %default]")

  p.add_option("--vertical", action="append", type=float, nargs=2, dest="vline", help="Draw a vertical line at x='vline[1]' in plot # vline[0]")

  return p.parse_args()

def check_optargs(opt, args):
  """
  Check if the input options and parameters are ok.
  
  Parameters
  ----------
  options: dictionary
    dictionary of options returned by optparse
  args: list
    list of arguments returned by optparse

  output:
    opt: dictionary of the options
    file_roots: list of file roots
    tags: list of tags
    paramnames: list of parameter names
  """

  #check and split the args
  if(len(args)<1):
    raise SystemExit("At least one argument required")
  if(opt.legend != None):
    if(len(args)%2 != 0):
      raise SystemExit("The legend is required: the second half of the argument list should contain the legend tags, but I have an odd number of args.")
    else:
      file_roots = args[:len(args)/2] #save the file roots
      tags = args[len(args)/2:]  #save the legend tags
  else:
    file_roots = args[:] #save the file roots
    tags = []  #save the legend tags

  #read parameter file names
  paramnames = cf.get_paramnames(file_roots, ext=opt.pext, verbose=opt.verbose)
  #if no column number is given, get it from std in
  if(opt.cols == None):
    for i,pn in enumerate(paramnames):
      print "%i) %s" %(i+2, pn.strip().split("\t")[0])
    temp = raw_input("Please give me the number corresponding to the parameters you want to plot separated by one or more space:\n")
    opt.cols = [int(i) for i in temp.split()] #save the chosen columns
  for i in opt.cols:
    if(i<2 or i > len(paramnames)+1):
      raise IndexError("Index out of the range [2,%d]" %(len(paramnames)+1))
  opt.cols.insert(0,0)
  paramnames = [paramnames[i] for i in np.array(opt.cols[1:])-2] #rearrange the parameter names and discard the unused ones

  opt.ncols = int(np.ceil(len(opt.cols[1:])/float(opt.nrows)))   #compute the number of columns

  if(opt.nrows == 1):   #adjust the boundaried 
    opt.bounds[5] = 0.
  if(opt.ncols == 1):
    opt.bounds[4] = 0.

  if(opt.lsfont == None):    #set legend font size if not given
    opt.lsfont = opt.sfont
  if(opt.lsfont == None):    #set tick labes font size if not given
    opt.tsfont = opt.sfont*0.7

  opt.fsize = tuple(np.array(opt.fsize)/mpm.inc2cm)   #convert the figure size in inches

  vertical = [[None, 0.] for i in opt.cols[1:]] #initialise the list of lists for the vertical lines
  if(opt.vline != None):   #if a tuple of values to draw the vertical line is given order it
    for v in opt.vline:
      vertical[int(v[0])-1] = list(v)
  opt.vline = vertical   #replace the option with the full list

  return opt, file_roots, tags, paramnames

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
    plt.setp(txt, fontsize=opt.lsfont)
  else:
    print "Keep default font size of the legend"
  if(kargs.has_key("textcol")==True and kargs["textcol"] == True):
    for i,l in enumerate(lines):
      txt[i].set_color(list(l)[0].get_color())

  pass

#############################################
###                 Main                  ###
#############################################

if __name__ == "__main__":   # if is the main

  (opt, args) = options(op.OptionParser(version="%prog version 0.1"))   #create the object optparse

  opt, froots, tags, paramnames = check_optargs(opt, args)   #check the input
  
  c_s, r_s = (opt.bounds[2]-opt.bounds[0])/float(opt.ncols), (opt.bounds[3]-opt.bounds[1])/float(opt.nrows)  #dimentions of the columns dnd rows
  lines=[]  #initialise list of lines for the legend, if required

  fig=plt.figure(1, figsize=opt.fsize)  #open the figure
  fig.subplots_adjust(left=opt.bounds[0], bottom=opt.bounds[1], right=opt.bounds[2], top=opt.bounds[3], wspace=opt.bounds[4], hspace=opt.bounds[5])

  for i,f in enumerate(froots):  #goes to the files
    chain = list(np.loadtxt(f+opt.cext, usecols=opt.cols).T)   #read the file containing the chain
    weights = chain[0]    #separate the weights from the values of the chains
    chain = chain[1:]
    for j,c,p,v in it.izip(it.count(),chain,paramnames,opt.vline):   #loop through the parameters to plot
      h, xe = np.histogram(c, bins=opt.n_bins, weights=weights)   #compute the histogram
      h /= np.amax(h)   #rescale the istogram in order to have peak =1
      xe = (xe[1:]+xe[:-1])/2.   #compute the mean of each bin

      #x_st, y_st = j%opt.ncols, opt.nrows-1-np.floor(j/opt.ncols)   #sublplot index
      #spl = fig.add_axes([opt.bounds[0]+c_s*x_st+opt.bounds[4], opt.bounds[1]+r_s*y_st+opt.bounds[5], c_s-opt.bounds[4], r_s-opt.bounds[5]])   #subplot
      spl = fig.add_subplot(opt.nrows, opt.ncols, j+1)

      l = spl.plot(xe, h, ls=mpm.linestyles[i%mpm.numls], c=mpm.colors[i%mpm.numcol], lw=opt.lw)
      if(j==0):  #save the line from the first parameters from each input file
        lines.append(l)
      if(i==0):  #for the first chain, write the x and y labels
	spl.set_xlabel("$"+p.strip().split("\t")[1]+"$", fontsize=opt.sfont)  #x axis names
	if(j%opt.ncols == 0):
	  spl.set_ylabel(r"$\mathcal{L}/\mathcal{L}_{\mathrm{max}}$", fontsize=opt.sfont)  #axis names
	else:
	  spl.set_yticklabels([])
	yr= spl.get_ylim()  #get y limits
	spl.set_ylim(yr[0],yr[1]*1.05)  #set ylimits
	#do not show the outmost x tick labes
	xmin,xmax = spl.get_xlim()
	xtl = spl.get_xmajorticklabels()
	tl = spl.xaxis.get_majorticklocs()
	xtl[(tl<xmin).sum()].set_visible(False)
	xtl[-(tl>xmax).sum()-1].set_visible(False)
	for xt in xtl:
	  xt.set_fontsize(opt.tsfont)
	for yt in spl.get_ymajorticklabels():
	  yt.set_fontsize(opt.tsfont)
	if(v[0] != None):
	  spl.plot([v[1],v[1]], spl.get_ylim(), "k--")

  if(opt.legend != None):  #write the legend if required
    #x_st, y_st = opt.legeng[0]%opt.ncols, opt.nrows-1-np.floor(opt.legeng[0]/opt.ncols)   #sublplot index
    #legend(fig.add_axes([opt.bounds[0]+c_s*x_st+opt.bounds[4], opt.bounds[1]+r_s*y_st+opt.bounds[5], c_s-opt.bounds[4], r_s-opt.bounds[5]]), lines, tags, opt.legend[1], opt)   #draw the legend
    legend(fig.add_subplot(opt.nrows, opt.ncols, opt.legend[0]), lines, tags, opt.legend[1], **opt.__dict__)   #draw the legend


  if(opt.outfile == None):
    plt.show()
  else:
    plt.savefig(opt.outfile)

  exit()
