#!/usr/bin/python
# -*- coding: utf-8 -*-

"""This program computes the correlation matrix from a list of files,
extract the variance and make a plot of the variances and of all the 
differences of the correlation matrices.
"""
import glob
import itertools as it
import matplotlib as mpl  #plotting stuff
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
  p.set_usage("""%prog [options] fname_cov1 [fname_cov2 ... fname_covn]
  Given n file names also with wildcards containg the covariance matrix, 
  it computes the correlation matrix and extract the variance. According 
  to the required make a plot of the variance, correlation and/or the
  variances and n(n-1)/2 differences of the correlation matrices.
  The default is the covariance. If the difference required, he name
  of at least two files must be given. The program is expecting a
  file with the structure: i, j, ki, kj, covariance. The covariance
  matrices need to be evaluated at the same wave number.""")

  p.add_option("--verbose", action="store_true", dest="verbose", help="Produces more output.")  # verbose option

  p.add_option("-k", "--krange", action="store", dest="krange", type=float, nargs=2, default=[None,None], help="Range in k to be considered. If the option is not used, the range will be automatically determined")   #x range values from the user

  p.add_option("-m", "--matrix", action="store", dest="matrix", help="If a valid file name is given, covariances are assumed to be saved a square matrices, whose k values are taken from the first column of 'matrix'")

  p.add_option("-l", "--legend", action="store", dest="legend", nargs=4, help="If give, from the file names is estracted a string between legend[0] and legend[1] and put between legend[2] and legend[3]. This list of strings is used for legend.")

  p.add_option("--vfigure-size", action="store", type=float, nargs=2, dest="vfsize", default=[10.,10.], help="x and y size of the plot of the variance in cm. [Default: %default]")
  p.add_option("--vbounds", action="store", dest="vbounds", type=float, nargs=4, default=[0.20, 0.20, 0.98, 0.98], help="Boudaries of the variance plot [left, bottom, right, top]. [Default: %default]")

  p.add_option("--cfigure-size", action="store", type=float, nargs=2, dest="cfsize", help="x and y size of the plot of the correlation matrices in cm.")
  p.add_option("--cbounds", action="store", dest="cbounds", type=float, nargs=6, help="Boudaries of the correlation plot [left, bottom, right, top, space, bar_width]. The 'right' keyword refers to the right limit for the covariances; the color bar will be outside.")

  p.add_option("--dfigure-size", action="store", type=float, nargs=2, dest="dfsize", default=[10.,10.], help="x and y size of the plot of the differnce of the correlation matrices in cm.")
  p.add_option("--dbounds", action="store", dest="dbounds", type=float, nargs=6, help="Boudaries of the difference of the covariances plot [left, bottom, right, top, space, bar_width]. The 'right' keyword refers to the right limit for the covariances; the color bar will be outside.")
  
  p.add_option("--cmap", action="store", type="string", dest="cmap", default="gist_rainbow_r", help="Matplotlib colormap name. See e.g. 'http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps' for a list of colormaps. [Default: %default]")

  p.add_option("-f", "--font-size", action="store", type=int, dest="fsize", default="15", help="Axis font size. [Default: %default]")
  p.add_option("-w", "--line-width", action="store", type=float, dest="lwidth", default=1, help="Line width. [Default: %default]")

  p.add_option("-o", "--output-file", action="store", dest="outfile", type="string", help="Output file name root.")   #output file name
  p.add_option("-e", "--extension", action="store", type="string", dest="ext", default="eps", help="Extension of the output files. [Default: %default]")

  p.add_option("-v", "--variance", action="store_true", dest="var", default=False, help="Plot the variances.")
  p.add_option("-c", "--covariance", action="store_true", dest="cov", default=False, help="Plot the covariances. If none of the other is given, the covariance is plotted.")
  p.add_option("-d", "--diffcov", action="store_true", dest="diffcov", default=False, help="Plot the difference of the covariances.")

  p.add_option("--n-cols", action="store", dest="n_cols", default=3, type=int, help="Number of columns in the plots of the correlation matrices. [Default: %default]")
  
  p.add_option("--ref", action="store", dest="ref", type=int, help="If the difference of covariances is required and this option give, the 'ref' elements is subtracted from all the other covariances.")

  return p.parse_args()

def files2corr(fnames, kmin=None, kmax=None, verbose=False, matrix=None):
  """Given a file name read the columns with k and covariance, 
  or the file with the k and the covariance matrix,
  cut the covariance according to the range, if given, and returns
  the the krange, the variance, correlation matrix and the range of the latter two

  Parameters:
  fname: list of strings
    file names
  kmin, kmax: two elemnts list (optional)
    range of the covariance matrix

  Output: 
    (k, kmin, kmax), (variance, varmin, varmax), (correlation matrix cormin, cormax)
  """

  k, var, cor = [],[],[]   #initialise the lists of k values, variances, correlation matrices
  kmm, varmm, cormm = [np.inf, -np.inf], [np.inf, -np.inf], [np.inf, -np.inf]  #extreme
  if(matrix!=None):   #if the the input files are already matricies read the file containing the values of k
    kk = np.loadtxt(matrix, usecols=[0])   #read the values of k
      
  for fn in fnames:
    if(options.verbose):
      print "reading file '%s', extracting the variance and computing the correlation matrix" %fn
    if(matrix==None):   #if the covariance is in a column of the file
      kt, cov = np.loadtxt(fn, usecols=[3,4]).T   #read the covariance
      dim_cov = np.sqrt(cov.size)  #get the dimension of the covariance assuming that is square
      cov = np.reshape(cov, (dim_cov, dim_cov))   #make the covariace squared
      kt = kt[:dim_cov]  #get the kvalues
    else:   #if the covariance is already a matrix
      cov = np.loadtxt(fn)
      kt = kk
    if(kmax!=None):   #if the maximum value of k is given cut k and cov according to it
      kt = kt[kt<kmax]
      cov = cov[:kt.size, :kt.size]
    if(kmin!=None):   #if the maximum value of k is given cut k and cov according to it
      temp = kt[kt<kmin]
      cov = cov[temp.size:, temp.size:]
      kt = kt[temp.size:] 
    tvar = np.sqrt(np.diag(cov))
    for i,j in it.product(range(kt.size), repeat=2):   #loops over the two dimensions of the covariance
      cov[i,j] /= tvar[i]*tvar[j]   #compute the correlation matrix Cij/sqrt(Cii*Cjj)

    k.append(kt)
    var.append(tvar)  #append the variance and the correlation matrices to the output lists
    cor.append(cov)

    kmm = [min(kmm[0], np.amin(kt)), max(kmm[1], np.amax(kt))]   #find the mimimum un maximum of the variance and covariance
    varmm = [min(varmm[0], np.amin(tvar)), max(varmm[1], np.amax(tvar))]   #find the mimimum un maximum of the variance and covariance
    cormm = [min(cormm[0], np.amin(cov)), max(cormm[1], np.amax(cov))]

  return (k, kmm), (var, varmm), (cor, cormm)  #return

def get_tags(flist, begin="pk", end=".fkp", pre="$p_{\mathrm{w}}=", post="$"):
  """Given a list of string extract the part between "begin" and "end" and return it

  Paramters:
  ---------
  flist: list
    list of strings
  begin: string (optional)
    string to be found before the different part
  end: string (optional)
    string to be found after the different part
  pre: string (optional)
    string preceding to the extracted tag
  post: string (optional)
    string following to the extracted tag

  output: list
    list of strings containing the difference
  """
  output=[]
  for i in flist:
    output.append(pre+i[i.find(begin)+len(begin):i.find(end)]+post)   #extract the part between begin and end
  return output

def cor2diff(cor, ref=None, tag=None):
  """Given a list of arrays containing the correlation matrices (at least 2)
  return the n(n-1)/2 differences as a list of arrays

  Paramaters
  ----------
  cor: list of arrays
    list of correlation matrices, all with the same shape
  ref: integer (optional)
    if given cor[ref] subtracted from all the others
  tag: None|list (optional)
    if given the tags for the differences created and returned

  output: list of arrays, list of two elements, list or None
    n(n-1)/2 with (n=len(cor)) containing the differences of the correlation matrices and 
    the absolut maximum and minimum of the differences. The last element is returned as the tags of the difference if
    tag is given and none otherwise
  """
  if(ref!=None and ref>len(cor)-1):
    print "The index %d is out of range" %ref
    exit()
  if(tag!=None and len(tag)!=len(cor)):
    print "The numbers of tags and of matrices do not correspond"
    exit()
  dcor=[]  #store the differences between correlation matrices
  dtag=[]  #tags for the differences
  dcormm=[1e+6,0.]
  if(ref==None):
    iterable = it.combinations(range(len(cor)),2)
  else:
    iterable = it.product(range(len(cor)),[ref])
  for i,j in iterable:   #compute the differences of the
    if(ref!=None and i==j):
      continue
    temp = cor[i] - cor[j]
    if(tag!=None):
      dtag.append(tag[i]+" - "+tag[j])
    dcormm= [min(dcormm[0], np.amin(temp)), max(dcormm[1], np.amax(temp))]
    dcor.append(temp)

  if(tag==None):
    return dcor, dcormm, None
  else:
    return dcor, dcormm, dtag


def ncr(n_el, m_c=3, m_r=None):
  """given the number of elements returns the number of columns and
  rows needed if m_c or m_r is the maximum number of columns or
  rows required. If both a given and n_el > m_c*m_r the extra elements
  not considered; otherwise adjust optimally the size prefering the
  columns. If the number of elements is smaller than the maximum, a single
  column/row returned. If m_c and m_r == None, a single column returned

  Parameters
  ----------
  n_el: integer
    number of elements
  m_c: integer (optional)
    maximum number of columns
  m_r: integer (optional)
    maximum number of rows

  output: list
    [n_columns, n_rows]
  """

  if(m_r==None and m_r==None):
    n_col, n_row = 1, n_el
  if(m_r==None):
    n_col = min(n_el,m_c)
    n_row = np.ceil(n_el/float(n_col))
  elif(m_c==None):
    n_row = min(n_el,m_r)
    n_col = np.ceil(n_el/float(n_row))
  else:
    if(n_el > m_c*m_r):
      n_col, n_row = m_c, m_r
    else:
      for i in range(1,min(m_c,m_r)):
        if(n_el > (m_c-i)*(m_r-i)):
          n_col, n_row = m_c-i+1, np.ceil(n_el/float(m_c-i+1))
          break

  return [int(n_col), int(n_row)]


def plotvar(fig, x, y, extent=None, outfile=None, leg=None, fsize=12, lw=1, bounds=None):
  """Given a figure object, an array of x values and a list of arrays of y values,
  it plot them on a single figure.

  Parameters
  ----------
  fig: figure object
    figure for the plot
  x: array 
    contains the x axis
  y: list of arrays
    contains the y axis: each element of y must be an array with the same dimension of x
  extent: list (optional)
    extrema of the plot
  outfile: string (optional)
    file name where to save the plot. If not given the figure is not saved and can be shown
  leg: list (optional)
    if not 'None' plot the legend with leg as tags. len(leg)==len(y) is required
  fsize: int (optional)
    size of the fonts used in the plot
  lw: float (optional)
    line width
  bounds: list
    bounds of the subplot

  ------
  output: nothing
  """
  if(leg!=None and len(y)!=len(leg)):
    print "The number of tags and of lines do not correspond"
    exit()

  if(bounds == None):
    fig.subplots_adjust(left=0.25, right=0.95, bottom=0.21, top=0.95) # set window size
  else:
    fig.subplots_adjust(left=bounds[0], right=bounds[2], bottom=bounds[1], top=bounds[3]) # set window size

  subf = fig.add_subplot(1,1,1, xscale="log", yscale="log")  #add a subplot

  lines=[] #it will contain a list of lines
  for i,k,v in it.izip(it.count(),x,y):
    lines.extend(subf.plot(k, v, mpm.linestyles[i%(mpm.numls+mpm.numsy)], color=mpm.colors[i%mpm.numcol], lw=lw))  #create the line

  subf.set_xlabel(r"$k\, [h/Mpc]$", fontsize=fsize)  #axis names
  subf.set_ylabel(r"$\sigma(k)\, [(Mpc/h])^3$", fontsize=fsize)
  subf.axis(extent)  # axix range
  if(leg!=None):
    l = subf.legend(lines, leg, loc= 'best')  #legend
    l.draw_frame(False)

  if(options.outfile != None):   #save the file if required
    fig.savefig(outfile, dpi=150) #output file
  pass

def plotmat(fig, mat, ncr, extent=None, vmin=None, vmax=None, outfile=None, leg=None, cmap='gray', fsize=12, bounds=None):
  """Given a figure object, an array of x values and a list of arrays of matrices 'mat',
  it plot them on a single figure on different subplots.

  Parameters
  ----------
  fig: figure object
    figure for the plot
  mat: list of matrices
    contains the y axis: each element of y must be a square matrix with the same dimensions of x
  ncr: list
    number of columns and rows to be plotted
  extent: list (optional)
    extent for matplotlib.pyplot.imshow
  vmin/vmax: scalar (optional)
    Used to scale a luminance image to 0-1. If either is None, the min and max of the luminance values will be used
  outfile: string (optional)
    file name where to save the plot. If not given the figure is not saved and can be shown
  leg: list (optional)
    if given and of the same lenght of mat, the string contained in leg are annotated on top of each image
  cmap: string (string)
    color map used for the plot
  fsize: int (optional)
    size of the fonts used in the plot
  bounds: list
    bounds[0:3]=bounds of the matrixes plots, bounds[4]=distance between them, bounds[5]=width of the bar

  ------
  output: nothing
  """

  if(bounds == None):
    if(n_cr[0]==1):
      llxp, urxp = 0.2, 0.7
    else:
      llxp, urxp = 0.13, 0.8
    llyp,uryp = 0.09, 0.98  #positions of the plots
    space = 0.01
    wbar = 0.05
  else:
    (llxp, llyp, urxp, uryp, space, wbar) = tuple(bounds)

  c_s, r_s = (urxp-llxp)/float(ncr[0]), (uryp-llyp)/float(ncr[1]) #dimentions of the columns dnd rows

  for i,c in enumerate(mat):
    x_st, y_st = i%ncr[0], ncr[1]-1-np.floor(i/ncr[0])
    box = [llxp+c_s*x_st+space, llyp+r_s*y_st+space, c_s-space, r_s-space]
    sfig = fig.add_axes(box)
    sfig.imshow(c.T, extent=extent, interpolation='nearest', origin='lower', vmin=vmin, vmax=vmax, cmap=cmap)
    if(leg!=None and len(leg)==len(mat)):
      sfig.text(np.mean(extent[:2]), extent[1]*1.03, leg[i], horizontalalignment='center', size=10)
    if(i%ncr[0] == 0):
      sfig.set_ylabel("k [h/Mpc]", fontsize=fsize)  #axis names
    else:
      sfig.set_yticklabels([])
    if(len(mat)-i <= ncr[0]):
      sfig.set_xlabel("k [h/Mpc]", fontsize=fsize)  #axis names
    else:
      sfig.set_xticklabels([])
    for label in sfig.get_xticklabels():
      label.set_fontsize(fsize)
    for label in sfig.get_yticklabels():
      label.set_fontsize(fsize)

  #colorbar
  if(ncr[1]==1):
    rect_bar = [urxp+3*space, llyp, wbar, r_s-3*space]  #positions of the bar
  else:
    rect_bar = [urxp+3*space, 0.5-r_s, wbar, 2*r_s-3*space]  #positions of the bar

  cbar=np.outer(np.arange(vmin,vmax,0.01),np.ones(10))   #crate the bar
  axBar = fig.add_axes(rect_bar)  #open the box for the plot
  axBar.imshow(cbar, aspect='auto', cmap=cmap, origin="lower", extent=(10,20,vmin,vmax))

  axBar.yaxis.set_ticks_position("right") # ticks and ticklabels
  axBar.yaxis.set_label_position("right") # axis label
  axBar.set_ylabel("Correlation", fontsize=fsize)
  for label in sfig.get_yticklabels():
    label.set_fontsize(fsize)

  #xticklabels = axBar.get_xticklabels()
  axBar.set_xticks([])

  if(options.outfile != None):
    fig.savefig(outfile, dpi=150) #output file

  pass

if __name__ == "__main__":   # if is the main

  (options, args) = options(op.OptionParser(version="%prog version 0.1"))   #create the object optparse

  if(options.var==False and options.cov==False and  options.diffcov==False):
    options.cov = True     #set the covariance as default if no other is required

  flist=[]   #convert args list in file name list
  for a in args:
    for f in glob.iglob(a):
      flist.append(f)

  if(len(flist)<1):
    print """At least one file name must be given.
    Type './%s -h' for more information""" % os.path.basename(sys.argv[0])
    exit()
  if(options.diffcov==True and len(flist)<2):
    print """To compute the difference of convariancees at least two file names must be given.
    Type './%s -h' for more information""" % os.path.basename(sys.argv[0])
    exit()

  mpm.linestyles.extend(mpm.symbols)   #merge lines and symbols

  n_args = len(args)   #number of input file names

  (k, kmm), (var, varmm), (cor, cormm) = files2corr(flist, kmin=options.krange[0], kmax=options.krange[1], verbose=options.verbose, matrix=options.matrix)   #give file name and max k and get back k, variance and correlation matrix

  if(options.legend==None):   #if no legend requred
    tag = None
  else:
    tag = get_tags(flist, begin=options.legend[0], end=options.legend[1], pre=options.legend[2], post=options.legend[3])

  if(options.var==True):
    if(options.verbose):
      print "Plotting the variances together"
    #plot the variances
    outfile = options.outfile  #output file name
    if(outfile!=None):
      outfile = options.outfile+"var." + options.ext
    plotvar(plt.figure(1, figsize=tuple(np.array(options.vfsize)/mpm.inc2cm)), k, var, extent=[kmm[0]*0.95,kmm[1]*1.05, varmm[0]*0.95,varmm[1]*1.05], outfile=outfile, leg=tag, fsize=options.fsize, lw=options.lwidth, bounds=options.vbounds)
    # open window figure

  mpl.rcParams['xtick.direction'] = 'out'
  mpl.rcParams['ytick.direction'] = 'out' 

  if(options.cov==True):
    if(options.verbose):
      print "Plotting the covariances"

    n_cr = ncr(len(cor), m_c=options.n_cols)
    if(options.cfsize == None):
      xs = ys = 9./mpm.inc2cm   # window size
      if(n_cr[0]==1):
	xs, ys = xs+4./mpm.inc2cm, ys*n_cr[1]
      if(n_cr[0]>1):
	xs, ys = xs*2., ys*n_cr[1]*2./n_cr[0]
      options.cfsize = (xs, ys)
    else:
      options.cfsize = tuple(np.array(options.cfsize)/mpm.inc2cm)
    outfile = options.outfile  #output file name
    if(outfile!=None):
      outfile = options.outfile+"cor." + options.ext
    plotmat(plt.figure(2, figsize=options.cfsize), cor, n_cr, extent=[kmm[0],kmm[1],kmm[0],kmm[1]], vmin=cormm[0], vmax=cormm[1], outfile=outfile, leg=tag, cmap=options.cmap, fsize=options.fsize, bounds=options.cbounds)


  if(options.diffcov==True):
    if(options.verbose):
      print "Plotting the covariance differences together"

    dif, difmm, dtag = cor2diff(cor, ref=options.ref, tag=tag)   #get the differences between covariances

    n_cr = ncr(len(dif), m_c=options.n_cols)
    if(options.dfsize == None):
      xs = ys = 9./mpm.inc2cm   # window size
      if(n_cr[0]==1):
	xs, ys = xs*1.2, ys*n_cr[1]
      elif(n_cr[0]>1):
	xs, ys = xs*2., ys*n_cr[1]*2./n_cr[0]
      options.dfsize = (xs, ys)
    else:
      options.dfsize = tuple(np.array(options.dfsize)/mpm.inc2cm)
    outfile = options.outfile  #output file name
    if(outfile!=None):
      outfile = options.outfile+"difcor." + options.ext
    plotmat(plt.figure(3, figsize=options.dfsize), dif, n_cr, extent=[kmm[0],kmm[1],kmm[0],kmm[1]], vmin=difmm[0], vmax=difmm[1], outfile=outfile, leg=dtag, cmap=options.cmap, fsize=options.fsize, bounds=options.dbounds)

  if(options.outfile == None):
    plt.show()

  exit()

