#!/usr/bin/python3
# -*- coding: utf-8 -*-

#import itertools as it
#import matplotlib.pyplot as plt  #plotting stuff
#import matplotlib.ticker as tic  #axis formatter
#import my_statistic as ms
#import myplotmodule as mpm   # my model used to plot
#import os    #contain OS dependent stuffs: it helps with sistem portability  

import contour_func as cf
import numpy as np
import sys 


def parse(argv):
    """
    This function accept a list of strings, create and fill a parser istance 
    and return a populated namespace

    Parameters
    ----------
    argv: list of strings
        list to be parsed

    output: namespace
    ---------
    """

    import argparse as ap 
    import argparse_custom as apc

    description = """    This program plot the values of the given parameters as
    function of the maximum scale used to perform the fit. 
    To build the file names first the maximum scale is computed:
    k_max=numpy.arange(range[0],range[1]+1,stride)/100.
    Then the full file name is build if the following way:
      1) if 'froot' contains only one '{}', this is substituted with
      k_max and the extension is appended
      2) otherwise the file name is build as froot+kmax+extension
    If more than a file root given and all the columns are plotted,
    the files must have the same structure. If only selected columns 
    are desired, then they must be present in all the files.
    """

    p = ap.ArgumentParser(description=description,
            formatter_class=apc.RawTextArgDefHelpFormatter)

    p.add_argument('range', nargs=2, type=float, action='store',
            help="Maximum and minimum wavenumber multiplied by 100")
    p.add_argument('froot', nargs='+', action='store',
            help="Root of the file names")

    p = apc.version_verbose(p, '2')

    p.add_argument("--paramlist", action="store_true", default=False,
            help="Read the parameter list from the fist file, print it out and close the program")

    # file related options
    pf = p.add_argument_group(description='Input-output file options')

    pf.add_argument("--stride", action="store", type=int, default="1",
            help="""Stride in 'range'""")

    pf.add_argument("--ext-chain", action="store", default="0.5_total.txt",
            help="""Extension of the input file containing the chain.""")
    pf.add_argument("--ext-parmn", action="store", default="paramnames",
            help="""Extension of the input file containing the parameters' names.""")

    pf.add_argument("-s", "--skip", action="store_true", default=False,
            help="Skip non existing files")

    pf.add_argument("-o", "--ofile", action="store", 
            help="Output file name. If not given, the plot is shown on screen.")

    # plot options
    pp = p.add_argument_group(description="Plot related options")

    pp.add_argument("-c", "--columns", action="store", nargs='+', 
            help="""Name of the columns to plot, as appears in the first column of
each parameter file name. If not given, all columns are read""")

    pp.add_argument("-n", "--nsigma", action="store", type=int, default="1",
            help="Plot the ns-sigma error bars.")

    pp.add_argument("-m", "--marker-size", action="store", type=float,
            default="4", help="Marker size.")
    pp.add_argument("--font-size", action="store", type=int, 
            default="15", help="Axis font size.")
    pp.add_argument("-w", "--line-width", action="store", type=float, default=1,
            help="Line width.")

    pp.add_argument("--shift", action="store", type=float, default=0,
            help="""Displace the errorbars by '%(dest)s' to avoid overlap.
Ignored if 'fill-between' is used""")

    pp.add_argument("-f", "--fill-between", action=apc.required_range(0, 1), 
            type=float, nargs=2, 
            help="""Substitute errorbars with matplotlib 'fill_between' 
with transparencies in the range [0,1]. 
No transparency for eps.""")

    pp.add_argument("--horizontal", type=float, nargs='+',
            action=apc.multiple_of(2, reshape=True), 
            help="""Plot an horizontal line in subplot '%(dest)s[0]' at
y='%(dest)s[1]'. Multiple lines can be drawn providing couple of
subplot-line.""")

    pp.add_argument("-b", "--bounds", action="store", type=float, nargs=4,
            default=[0.20, 0.20, 0.97, 0.97], 
            help="""Left, bottom, right and top limits of the plotting area
within the figure window.""")

    pp.add_argument("-r", "--rescale", type=float, nargs='+',
            action=apc.multiple_of(2, reshape=True), 
            help="""Rescale the abscissa in subplot '%(dest)s[0]' by
'%(dest)s[1]'. The same rescaling factor is shown in the y label.
Multiple rescaling can be drawn providing couple of
subplot-rescaling.""")
            
    pp.add_argument("--y-label", nargs='+', action=apc.multiple_of(2, reshape=True), 
            help="""Set y axis label in subplot '%(dest)s[0]' to '%(dest)s[1]'.
Multiple labels can be drawn providing couple of 
subplot-label.""")

    pp.add_argument("--y-range", nargs='+', action=apc.multiple_of(2, reshape=True), 
            help="""Set y axis range in subplot '%(dest)s[0]' to '%(dest)s[1]'.
Multiple ranges can be set providing couple of 
subplot-range.""")

    #legend options
    pl = p.add_argument_group(description='Legend options')

    pl.add_argument("-l", "--legend", nargs='+', 
            help="""Legend tags. If given, the number of elements must be the
same as the number of input root and in the same order.""")

    pl.add_argument("--legend-plot", action="store", type=int, default=0,
            help="""Put the legend in plot '%(dest)s'. If negative, figure
legend drawn.""")

    pl.add_argument("--loc", type=apc.int_or_str, default=0,
            help='Legend location (see matplotlib legend help)')

    pl.add_argument("--legend-fsize", type=float, 
            help="""Legend tags font size. Defaults to axis font size""")


    #correction of the variance due to number of mocks and bins: Percival ... 2013
    description = """Correct the estimated variance by sqrt(m1) or sqrt(m2), as
    in Percival et al. 2013"""
    pc = p.add_argument_group(description=description)
    with open('par_vs_k.ini', 'r') as f:
        lines = f.readlines()
    lines = [l.replace('%', '%%') for l in lines]
    lines.insert(0, "Inifile with the following structure\n\n")
    lines = "".join(lines)
    pc.add_argument('--ini', type=ap.FileType('r'), help=lines)

    return p.parse_args(args=argv)
#end def parse(argv)

class IniNumber(Exception):
    def __init__(self, key, n):
        """gets the keyword name and the number of elements it should have"""
        string = "Keyword '{}' can have only 1 or {} elements"
        self.string = string.format(key, n)
    def __str__(self):
        return repr(self.string)

def read_ini(ini, n_froot):
    """read inifile
    Parameters
    ----------
    ini: file object
        inifile to parse
    n_froot: int
        number of file roots
    output
    ------
    n_params: int
        number of model parameters
    n_mocks: list of int of len n_froot
        number of mocks used to create che covariance matrix for each froot
    k_from_files: list of len n_froot of ndarrays 
        firs columns of the files given in keyword pk_files
    ind_kmin: list of int of len n_froot
        indeces of the first elements in k_from_files larger than value(s) in
        the keyword kmin 
    """
    from configparser import ConfigParse
    # initialise and get the default section (no other used)
    config = ConfigParser()
    config.read_file(ini)
    config = config['default']  
    # get the parameters
    # n_params
    n_params = config.getint('n_params')
    # n_mocks
    n_mocks = [int(nm) for nm in config.get('n_mocks').split()]
    if len(n_mocks) == 1:
        n_mocks = n_mocks*n_froot
    if len(n_mocks) != n_froot:
        raise IniNumber("n_mocks", n_froot)
    # pk_files: read the first column
    pk_files = config.get('pk_files').split()
    if len(pk_files) == n_froot:
        k_from_files = [np.loadtxt(fn, usecols=[0,]) for fn in pk_files]
    elif len(pk_files) == 1:
        k_from_files = [np.loadtxt(pk_files[0], usecols=[0,]),]*n_froot
    else:
        raise IniNumber("pk_files", n_froot)
    # kmin
    kmin = [float(km) for km in config.get('kmin').split()]
    if len(kmin) == 1:
        kmin = kmin*n_froot
    if len(kmin) != n_froot:
        raise IniNumber("kmin", n_froot)
    else:
        ind_kmin = [(k<km).sum() for km, k in zip(kmin, k_from_files)]

    return n_params, n_mocks, k_from_files, ind_kmin

def kmax_correction(krange, stride, n_froot, ini=None):
    """
    Create the values of kmax and if ini is not none uses the information in
    the configuration file to estimate the correction factor needed to unbias
    the estimated variance as in Percival et al. 2013
    Parameters
    ----------
    krange: list of two floats
        minimum and maximum k value multiplied by 100 
        (in my applications these are in reality integers
    stride: int
        stride between min and max
    n_froot: int
        number of file root
    ini: file object
        configuration file that is passed to configparser
    output
    ------
    kmax: nd array
        maximum value of k used
    correction: nd array of shape k_max.shape
        correction to be applied according to the above paper. All 1 if ini==None
    """
    kmax = np.arange(int(args[0]),int(args[1])+1,options.stride)/100.  #create array with the k_max
    #array containing the corrections. All 1, unless filled if ini is not none
    correction = np.ones([kmax.size, n_froot]) 
    #read the ini file and compute the corrective factor
    if ini is not None:
        n_params, n_mocks, k_from_files, ind_kmin = read_ini(ini, n_froot)
        #Iterate trough corr and do the required operations
        # create the iterator with the indeces
        it = np.nditer(correction, flags=['multi_index'], op_flags=['writeonly'])
        while not it.finished: #loop untill the end of the iterator
            ind0, ind1 = it.multi_index # the the two indeces
            #number of bins in given file with give k_max
            nbins = np.sum(k_from_files[ind1]<k_max[ind0])-ind_kmin[ind1]
            # compute m1 and save into corr
            denominator = (n_mocks[ind1]-nbins-1.) * (n_mocks[ind1]-nbins-4.)
            A = 2./denominator
            B = (n_mocks[ind1]-nbins-2.)/denominator
            it[0] = (1.+B*(nbins-n_params)) / (1.+A+B*(n_params+1.))
            it.iternext()

    return np.tile(kmax[:,None], n_froot), correction

class FrootException(Exception):
    """exception related with file roots"""
    def __init__(self, string):
        self.string=string
    def __str__(self):
        return repr(self.string)
def make_full_root(froots, kmax):
    """
    Create a nd array of full file roots using the value of kmax
    Parameters
    ----------
    froots: list
        list of file roots
    kmax: nd array
        values of kmax
    output
    ------
    full_roots: 2d array
        full file roots with kmax inserted
    """
    #if a '{}' present substitute with '{0:3.2f}'
    #otherwise append it to froots
    to_subst = '{0:3.2f}'
    temp_roots = []
    for fr in froots:
        n_brackets = fr.count('{}') #count the number of brackets
        if n_brackets == 0:
            temp_roots.append(fr+to_subst)
        elif n_brackets == 1:
            temp_roots.append(fr.replace('{}', to_subst))
        else:
            raise FrootException("The input file roots must have zero or one '{}'")
    #make an empty string 2d numpy array 
    full_roots = np.empty_like(kmax, dtype=unicode)
    it = np.nditer([full_roots, kmax], flags=['multi_index'], op_flags=['readwrite'])
    while not it.finished: #loop untill the end of the iterator
        ind0, ind1 = it.multi_index # the the two indeces
        it[0] = temp_roots[j].format(it[1]) #substitute {} with kmax
        it.iternext()

    return full_roots

def print_paramlist(fname):
    """
    Print the parameter name file and exit
    Parameters
    ----------
    fname: string
        file name
    """
    with open(fname, 'r') as f:
        print(f.read())
    exit()

def skip(fchain, fparamnames, correction):
    """
    Check if the pairs of chain/paramnames files exist. If not set to -999 the
    both file names and the correction
    Parameters
    ----------
    fchain: nd array of strings
        chain files
    fparamnames: nd array of strings
        paramname files
    correction: nd array
        correction to apply to each covariance
    """
    it = np.nditer([fchain, fparamnames, correction], flags=['multi_index'],
            op_flags=['readwrite'])
    while not it.finished: #loop untill the end of the iterator
        try:  #if the files exists no problem opening them
            with open(it[0], 'r'), open(it[1], 'r'):
                pass
        except IOError: #if any of them does note exist set -999
            it[0] = '-999'
            it[1] = '-999'
            it[2] = -999
        it.iternext()
    return fchain, fparamnames, correction

def files_to_mean_std(paramname, chains, params=None, verbose=True, skip=False):
    """Read the paramname and chain files and compute mean and std for all the
    required columns in the chains
    Parameters
    ----------
    paramname: nd array
        parameter file name files
    chains: nd array
        chains array
    params: list (optional)
        name of the parameters of interest as in the first column of the
        parameter name files. If None all parameters saved
    verbose: bool (optional)
        more output
    skip: bool (optional)
        skip non existing files
    output
    ------
    """
    # these dictionaries will have as key the latex name of the variable 
    # and each value is a nd array of float of paramnames.shape 
    mean, std = {}, {}
    value_prototipe = np.empty_like(paramnames, dtype=float)

    it = np.nditer([paramname, chains], flags=['multi_index'],
            op_flags=['readonly'])
    while not it.finished: #loop untill the end of the iterator
        ind0, ind1 = it.multi_index # the the two indeces
        try:
            paramindex = cf.get_paramnames(it[0], params=params, ext="",
                verbose=verbose, skip=args.skip)
            for pi in paramindex:

        except cf.ContourError:


        it.iternext()
    

if __name__ == "__main__":   # if is the main

    import sys
    args = parse(sys.argv[1:])

    # command line parameter checks
    if args.legend is not None and len(args.legend)!=len(args.froot):
        print("The number of legend tags must be the same as the number of file roots")
        exit()

    if args.verbose:
        print("Computing kmax and correction")
    k_max, corr = kmax_correction(args.range, args.stride, len(args.froot),
            ini=args.ini)

    if args.verbose:
        print("Create the the file names")
    args.froot = make_full_root(args.froot, k_max)

    #check that they have the same shape
    assert args.froot.shape == corr.shape, 
            "The array of file roots and corrections must have the same shape"

    #create the parameter and chain file names
    #substitute with np.add when will become available
    paramnames = np.core.defchararray.add(args.froot, '.'+args.ext_parmn)
    chains = np.core.defchararray.add(args.froot, '.'+args.ext_parmn)

    if args.paramlist:
        print_paramlist(paramnames[0])

    if args.skip:
        if args.verbose:
            print("Skipping non existing files")
        chains, paramnames, corr = skip(chains, paramnames, corr)

    #read the parameter name files

    #check that all the elements in "paramnames" are the same
    if args.columns is None:



    exit()
    

#  if(options.cols == None):   #if no columns given all the columns read
#    options.cols = list(np.arange(len(paramnames))+2)
#  else:
#    discard = []
#    for i in options.cols:   #check that the columns are correct
#      if(i<2 or i>len(paramnames)+1):
#        print "The column %d is not acceptable or do not exists. I'll skip it" %i
#	discard.append(i)
#    for i in discard:   #discard the non existing columns
#      options.cols.remove(i)
#    if(len(options.cols) == 0):
#      print "All the column number are too big or to small"
#      sys.exit(2)
#  options.cols.insert(0,0)  #read also the first column
#  n_plots = len(options.cols)-1   #number of plots
#
#  horiz = [[None, 0.] for i in range(n_plots)] #initialise the list of lists for the horizontal lines
#  if(options.horiz != None):   #if a tuple of values to draw the horizontal line is given order it
#    for h in options.horiz:
#      horiz[int(h[0])] = list(h)
#  rescale = [1 for i in range(n_plots)] #initialise the list of lists for the rescaling of the plot y axes
#  if(options.rescale != None):   #if a tuple of values to draw the horizontal line is given order it
#    for h in options.rescale:
#      rescale[int(h[0])] = h[1]
#  ylab = [[None, None] for i in range(n_plots)] #initialise the list of lists for the custom y labels
#  if(options.ylab != None):   #if a tuple of values to draw the horizontal line is given order it
#    for y in options.ylab:
#      ylab[int(y[0])] = list(y)
#  yrange = [[None, None, None] for i in range(n_plots) ] #initialise the list of lists for the custom y limits
#  if(options.yrange != None):
#    for y in options.yrange:
#      yrange[int(y[0])] = list(y)
#
#
#  #create the numpy arrays that will contain the mean and the stddev of the columns
#  mean = np.empty([len(files), n_plots])
#  stddev = np.empty([len(files), n_plots])
#
#  #initialise the figure
#  xs=10./mpm.inc2cm   # window size
#  ys= 2.5*(n_plots+1)/mpm.inc2cm
#  fig=plt.figure(1, figsize=(xs,ys))  #open the figure
#  rheight = (options.bounds[3]-options.bounds[1])/n_plots     #height of each subplot
#
#  if(options.fill != None):
#    alphas = np.linspace(options.fill[0], options.fill[1], len(args[2:]))    #alpha valued to be used in the fill_between
#  else:
#    alphas = [1 for i in args[2:]]
#  tk_max = np.copy(k_max)  #make a copy of the k_max to make the plots
#  lines=[] #initialise the line list
#  for i,fr,col,ls,m,a in it.izip(it.count(), args[2:], it.cycle(mpm.colors),
#      it.cycle(mpm.linestyles), it.cycle(mpm.symbols[2:-2]), alphas):   #loop through the file roots
#
#    files = [fr,]*k_max.size  # chain file names
#    for j,km in enumerate(k_max):   #create the file names
#      files[j] += "%3.2f.%s" %(km,options.ec)
#    if(options.skip == True):  #check if to skip files
#      files, discard = skip(files, slist=True)
#      tk_max = np.delete(k_max, discard)   #discard values of kmax if some of the file is not there
#
#    #read each chain and obtain mean and std
#    for j, fn in enumerate(files):
#      if(options.verbose == True):
#	print "Reading file %s" % fn
#      chain = np.loadtxt(fn, usecols=options.cols)  #read
#      if(options.verbose == True):
#	print "Computing the mean of the parameters in file %s" % fn
#      mean[j,:] = np.average(chain[:,1:], axis=0, weights=chain[:,0])  #mean of the columns containing parameters
#      if(options.verbose == True):
#	print "Computing the standard deviation of the parameters in file %s\n" % fn
#      stddev[j,:] = ms.stddev(chain[:,1:], weights=chain[:,0], axis=0)  #stddev of all the columns containing parameters
#
#    for j,h,r,yl,yr,c in it.izip(it.count(), horiz, rescale, ylab, yrange, options.cols[1:]):   #loop through the various columns to plot
#      box = [options.bounds[0], options.bounds[1]+(n_plots-j-1)*rheight,
#	  options.bounds[2]-options.bounds[0], rheight]   #create the plot box
#      subpl = fig.add_axes(box)  #create a subplots
#
#      if(np.absolute(r-1) > 1e-3):   #if a rescale required
#        mean[:,j] *= r
#        stddev[:,j] *= r
#      if(options.fill == None):
#	if(j==0):
#	  lines.append(subpl.errorbar(k_max+i*options.shift, mean[:,j],
#	    yerr=stddev[:,j], c=col, ls=ls, lw=options.lwidth, marker=m,
#	    ms=options.smarker)[0])
#	else:
#	  subpl.errorbar(k_max+i*options.shift, mean[:,j], yerr=stddev[:,j],
#	      c=col, ls=ls, lw=options.lwidth, marker=m, ms=options.smarker)
#      else:
#	if(j==0):
#	  lines.append(subpl.plot(k_max+i*options.shift, mean[:,j], c=col,
#	    ls=ls, lw=options.lwidth, marker=m, ms=options.smarker))
#	else:
#	  subpl.plot(k_max+i*options.shift, mean[:,j], c=col, ls=ls,
#	      lw=options.lwidth, marker=m, ms=options.smarker)
#	if(options.outfile != None and os.path.splitext(options.outfile)[1] == '.eps'):
#	  subpl.fill_between(k_max+i*options.shift, mean[:,j]+stddev[:,j],
#	      mean[:,j]-stddev[:,j], color=[140./255., 130./255., 255./255.,
#		1.], zorder=1)   #this is when eps is generated
#	else:
#	  subpl.fill_between(k_max+i*options.shift, mean[:,j]+stddev[:,j],
#	      mean[:,j]-stddev[:,j], color=col, alpha=a, edgecolor="w")
#      if(h[0] != None):
#	subpl.plot([k_max[0]*0.95,(k_max[-1]+i*options.shift)*1.02], [h[1],h[1]], 'k--')
#  
#      if(yl[0] == None):
#	label = "$"+paramnames[c-2].strip().split()[1]+"$"
#	if(np.absolute(r-1) > 1e-3):   #if a rescale required
#	  if(r.is_integer()==True):
#	    strr = str(int(r))
#	  else:
#	    strr = str(r)
#	  label = "$"+strr+"$"+label
#      else:
#	label = yl[1]
#      subpl.set_ylabel(label, fontsize=options.fsize)
#
#      if( yr[0] != None ):  #set the custom yrange
#	subpl.set_ylim( bottom=float(yr[1]), top=float(yr[2]) )
#      #exclude the first and the last visible tick labels
#      ymin,ymax = subpl.get_ylim()
#      ytl = subpl.get_ymajorticklabels()
#      tl = subpl.yaxis.get_majorticklocs()
#      ytl[(tl<ymin).sum()].set_visible(False)
#      ytl[-(tl>ymax).sum()-1].set_visible(False)
#      if(len(ytl)-2>5):
#        for iy in range(2,len(ytl)-1,2):
#	  ytl[iy].set_visible(False)
#      for iy in range(len(ytl)):   #set tick label font size
#	ytl[iy].set_fontsize(options.fsize)
#      if(j < n_plots-1):
#	subpl.set_xticklabels([])
#      else:
#	subpl.set_xlabel("$k_{\mathrm{max}}\,[h/Mpc]$", fontsize=options.fsize)
#	xtl = subpl.get_xmajorticklabels()  #get xlabels and change font size
#	for ix in range(len(xtl)):   #set tick label font size
#	  xtl[ix].set_fontsize(options.fsize)
#      for label in subpl.get_xticklabels():
#	label.set_fontsize(options.fsize)
#      for label in subpl.get_yticklabels():
#	label.set_fontsize(options.fsize)
#
#      subpl.set_xlim([k_max[0]*0.95, (k_max[-1]+i*options.shift)*1.02])  #set the axis size
#
#  if(options.legend != None and options.fleg == None):   #draw the legend in the desired axes
#    box = [options.bounds[0],
#	options.bounds[1]+(n_plots-options.legend-1)*rheight,
#	options.bounds[2]-options.bounds[0], rheight]   #create the plot box
#    subpl = fig.add_axes(box)  #create a subplots
#    l = subpl.legend(lines, leg, loc="best", numpoints=1, borderpad=0.3,
#	columnspacing=0.8, handlelength=1.4, borderaxespad=0.1,
#	labelspacing=0.2)
#    l.draw_frame(False)
#  if( options.fleg != None):  #draw figure legend
#    l = fig.legend(lines, leg, loc=options.fleg, numpoints=1, borderpad=0.3,
#	columnspacing=0.8, handlelength=1.4, borderaxespad=0.1,
#	labelspacing=0.2, frameon=False)
#
#  if(options.outfile == None):
#    plt.show()
#  else:
#    plt.savefig(options.outfile)
#
#  sys.exit()
