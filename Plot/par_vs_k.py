#!/usr/bin/python3
# -*- coding: utf-8 -*-

#import matplotlib.ticker as tic  #axis formatter
#import os    #contain OS dependent stuffs: it helps with sistem portability  

import contour_func as cf
import itertools as it
import matplotlib.pyplot as plt
import myplotmodule as mpm
import my_statistic as ms
import numpy as np
import sys 
from warnings import warn

def get_ini():
    """
    Sample ini file
    """
    inifile = """
    [default]
    # test config parser for Python/Plots/par_vs_k.py
    # number of parameters
    n_params = 4
    # number of mocks used: can be either 1 element (the same for all froot) or
    # as many as the number of froot. If negative the correction not computed
    #n_mocks = 600
    n_mocks = -1 300
    # if the covariance is computed as the mean of two covariances with
    # correlation r_mocks correct for it. Can be either 1 element (the same for
    # all froot) or as many as the number of froot. If not given, or negative
    # or correction ignored
    #r_mocks = 0.33
    r_mocks = -1 0.33

    # measurement file name used to compute the number of bins used (assumed to
    # be the first column) 
    # Can be either one or as many as the number of froot
    pk_dir = /data01/montefra/PThalos/V5.2/PS_win
    pk_files = %(pk_dir)s/ps_pthalos-v5.2-irmean.south.4000.1200.500.pw20000.masTSC.corr2.wmissboss.sectcomp.even.dat_rancorr.dat 
               %(pk_dir)s/ps_pthalos-v5.2-irmean.south.4000.1200.100.pw20000.masTSC.corr2.wmissboss.sectcomp.even.dat_rancorr.dat 
    # minimum value of k used. Can be either one or as many as the number of froot
    kmin = 0.02
    #kmin = 0.02 0.05
    """
    return inifile

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
    import textwrap as tw

    description = tw.dedent("""\
    This program plot the values of the given parameters as
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
    """)

    p = ap.ArgumentParser(description=description,
            formatter_class=apc.RawTextArgDefHelpFormatter)

    p.add_argument('range', nargs=2, type=float, action='store',
            help="Maximum and minimum wavenumber multiplied by 100")
    p.add_argument('froot', nargs='+', action='store',
            help="Root of the file names")

    p = apc.version_verbose(p, '2')

    p.add_argument("--paramlist", action="store_true", default=False,
            help=tw.dedent("""\
                    Read the parameter list from the fist file, print it out and
                    close the program"""))

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

    pp.add_argument("--variance-only", action="store_true", 
            help="""Plot only the variance.""")

    pp.add_argument("-c", "--columns", action="store", nargs='+', 
            help=tw.dedent("""\
                    Name of the columns to plot, as appears in the first column
                    of each parameter file name. If not given, all columns are read"""))

    pp.add_argument("-n", "--nsigma", action="store", type=int, default="1",
            help="Plot the ns-sigma error bars.")

    pp.add_argument("--set-MNRAS", action='store_true', 
            help='Use MNRAS presets for plotting')
    pp.add_argument("-m", "--marker-size", action="store", type=float,
            default="4", help="Marker size.")
    pp.add_argument("--font-size", action="store", type=int, 
            default="15", help="Axis font size.")
    pp.add_argument("-w", "--line-width", action="store", type=float, default=1,
            help="Line width.")

    pp.add_argument("--shift", action="store", type=float, default=0,
            help=tw.dedent("""\
                    Displace the errorbars by '%(dest)s' to avoid overlap.
                    Ignored if 'fill-between' is used"""))

    pp.add_argument("-f", "--fill", action=apc.required_range(0, 1), 
            type=float, nargs=2, help=tw.dedent("""\
                    Substitute errorbars with matplotlib 'fill_between' with
                    transparencies in the range [0,1].  No transparency for
                    eps."""))

    pp.add_argument("--horizontal", nargs='+',
            action=apc.multiple_of(2, reshape=True), help=tw.dedent("""\
                    Plot an horizontal line in subplot '%(dest)s[0]' at
                    y='%(dest)s[1]'. The sublot is identified with the short
                    name from the parameter file. Multiple lines can be drawn
                    providing couple of subplot-line."""))

    pp.add_argument("--figsize", type=float, nargs=2,
            action=apc.Cm2Inch, help="Figure size in cm")
    pp.add_argument("-b", "--bounds", action="store", type=float, nargs=4,
            default=[0, 0, 1, 1], help=tw.dedent("""\
                    (left, bottom, right, top) in the normalized figure
                    coordinate passed to 'plt.tight_layout'"""))

    pp.add_argument("-r", "--rescale", nargs='+',
            action=apc.multiple_of(2, reshape=True), help=tw.dedent("""\
                    Rescale the abscissa in subplot '%(dest)s[0]' by
                    '%(dest)s[1]'. The same rescaling factor is shown in the y
                    label. The sublot is identified with the short name from
                    the parameter file. Multiple rescaling can be drawn
                    providing couple of subplot-rescaling."""))
            
    pp.add_argument("--x-label", default="$k_{\mathrm{max}}\,[h/Mpc]$", help='x axis label')
    pp.add_argument("--y-label", nargs='+', action=apc.multiple_of(2, reshape=True), 
            help=tw.dedent("""\
                    Set y axis label in subplot '%(dest)s[0]' to '%(dest)s[1]'.
                    The sublot is identified with the short name from the
                    parameter file. Multiple labels can be drawn providing
                    couple of subplot-label."""))

    pp.add_argument("--y-range", nargs='+', action=apc.multiple_of(3, reshape=True), 
            help=tw.dedent("""\
                    Set y axis range in subplot '%(dest)s[0]' to '%(dest)s[1]'.
                    The sublot is identified with the short name from the
                    parameter file. Multiple ranges can be set providing couple
                    of subplot-range."""))

    #legend options
    pl = p.add_argument_group(description='Legend options')

    pl.add_argument("-l", "--legend", nargs='+', help=tw.dedent("""\
            Legend tags. If given, the number of elements must be the same as
            the number of input root and in the same order."""))

    pl.add_argument("--legend-plot", action="store", help=tw.dedent("""\
            Put the legend in plot '%(dest)s'. The sublot is identified with
            the short name from the parameter file. If '%(dest)s' is not in the
            list of parameters, the figure legend is drawn. If not given, the
            legend is drawn in the first plot"""))

    pl.add_argument("--loc", type=apc.int_or_str, default=0,
            help='Legend location (see matplotlib legend help)')

    pl.add_argument("--legend-fsize", type=float, 
            help="""Legend tags font size. Defaults to axis font size""")


    #correction of the variance due to number of mocks and bins: Percival ... 2013
    description = tw.dedent("""\
            Correct the estimated variance by sqrt(m1) or sqrt(m2), as in
            Percival et al. 2013""")
    pc = p.add_argument_group(description=description)
    lines = get_ini().replace('%', '%%')
    lines = ["Inifile with the following structure:", lines]
    lines = "".join(lines)
    pc.add_argument('--ini', type=ap.FileType('r'), help=lines)


    class Save_Ini(ap.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            lines = get_ini()
            lines = lines.split('\n')
            lines = [l[4:] for l in lines[1:]]
            lines = '\n'.join(lines)
            with open("sample.ini", 'w') as f:
                f.write(lines)
            print("'sample.ini' written in the current directory")
            sys.exit(0)
    pc.add_argument('--print-ini', nargs=0, action=Save_Ini, 
            help="Save the inifile to file 'sample.ini' in the current directory")

    return p.parse_args(args=argv)
#end def parse(argv)


# custom errors and warnings
class IniNumber(Exception):
    def __init__(self, key, n):
        """gets the keyword name and the number of elements it should have"""
        string = "Keyword '{}' can have only 1 or {} elements"
        self.string = string.format(key, n)
    def __str__(self):
        return repr(self.string)
class IniNotFound(Exception):
    def __init__(self, key):
        """Key not found"""
        string = "Keyword '{}' has not been found"
        self.string = string.format(key)
    def __str__(self):
        return repr(self.string)
class MyWarning(UserWarning):
    pass

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
    n_mocks: list of tuples of size 2 of len n_froot
        number of mocks and correlation, if a mean computed, used to create che
        covariance matrix for each froot
    k_from_files: list of len n_froot of ndarrays 
        firs columns of the files given in keyword pk_files
    ind_kmin: list of int of len n_froot
        indeces of the first elements in k_from_files larger than value(s) in
        the keyword kmin 
    """
    from configparser import ConfigParser

    # initialise and get the default section (no other used)
    config = ConfigParser()
    config.read_file(ini)
    config = config['default']  

    # get the parameters
    # n_params
    if 'n_params' in config:
        n_params = config.getint('n_params')
    else:
        raise IniNotFound('n_params')

    # n_mocks
    if 'n_mocks' in config:
        n_mocks = [int(nm) for nm in config.get('n_mocks').split()]
        if len(n_mocks) == 1:
            n_mocks = [n_mocks[0] for n in range(n_froot)]
        if len(n_mocks) != n_froot:
            raise IniNumber("n_mocks", n_froot)
    else:
        raise IniNotFound('n_mocks')

    # correlation between mocks when computing the covariance as the average
    # of two. Not mandatory
    r_mocks = config.get('r_mocks')
    if r_mocks is None:
        r_mocks = [-1 for n in n_froot]
    else:
        r_mocks = [float(r) for r in r_mocks.split()]
        if len(r_mocks) == 1:
            r_mocks = [r_mocks[0] for n in range(n_froot)]
        if len(r_mocks) != n_froot:
            raise IniNumber("r_mocks", n_froot)

    # zip n_mocks and r_mocks and save to n_mocks
    n_mocks = list(zip(n_mocks, r_mocks))

    # pk_files: read the first column
    if 'pk_files' in config:
        pk_files = config.get('pk_files').split()
        if len(pk_files) == n_froot:
            k_from_files = [np.loadtxt(fn, usecols=[0,]) for fn in pk_files]
        elif len(pk_files) == 1:
            k_from_files = np.loadtxt(pk_files[0], usecols=[0,])
            k_from_files = [k_from_files.copy() for n in range(n_froot)]
        else:
            raise IniNumber("pk_files", n_froot)
    else:
        raise IniNotFound('pk_files')

    # kmin
    if 'kmin' in config:
        kmin = [float(km) for km in config.get('kmin').split()]
        if len(kmin) == 1:
            kmin = [kmin for n in range(n_froot)]
        if len(kmin) != n_froot:
            raise IniNumber("kmin", n_froot)
        else:
            ind_kmin = [(k<km).sum() for km, k in zip(kmin, k_from_files)]
    else:
        raise IniNotFound('kmin')

    return n_params, n_mocks, k_from_files, ind_kmin
#end def read_ini(ini, n_froot):

def kmax_correction(krange, stride, n_froot, shift=None, ini=None):
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
    shift: float
        shift in kmax of each file root with respect to the previous one to
        avoid superposition of errorbars
    ini: file object
        configuration file that is passed to configparser
    output
    ------
    kmax: nd array
        maximum value of k used
    correction: nd array of shape kmax.shape
        correction to be applied according to the above paper. All 1 if ini==None
    """
    kmax = np.arange(int(krange[0]),int(krange[1])+1,stride)/100.  #create array with the k_max
    kmax = np.tile(kmax[:,None], n_froot).T
    #array containing the corrections. All 1, unless filled if ini is not none
    correction = np.ones_like(kmax) 
    #read the ini file and compute the corrective factor
    if ini is not None:
        n_params, n_mocks, k_from_files, ind_kmin = read_ini(ini, n_froot)
        # Iterate trough correction and kmax and do the required operations
        # create the iterator with the indeces
        ndit = np.nditer([correction, kmax], flags=['multi_index'],
                op_flags=['readwrite'])
        while not ndit.finished: #loop untill the end of the iterator
            ind0, ind1 = ndit.multi_index # the two indeces
            # number of bins in given file with give k_max
            nbins = np.sum(k_from_files[ind0]<ndit[1])-ind_kmin[ind0]
            # compute A and B
            n_m, r_m = n_mocks[ind0]
            if n_m <= 0:
                ndit[0] = 1.
            else:
                denominator = (n_m-nbins-1.) * (n_m-nbins-4.)
                A = 2./denominator
                B = (n_m-nbins-2.)/denominator
                if r_m >= 0: # do the correction if there is any correlation
                    r_c = (1.+r_m**2)/2.
                    A *= r_c
                    B *= r_c
                # m1 and save into correction
                ndit[0] = np.sqrt((1.+B*(nbins-n_params)) / (1.+A+B*(n_params+1.)))
            ndit.iternext()

    # shift the kmax
    if shift is not None:
        shifta = np.arange(n_froot)*shift
        kmax += shifta[:,None]

    return kmax, correction
#end def kmax_correction(krange, stride, n_froot, ini=None):

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
    full_roots = np.tile(np.array(temp_roots)[:,None], kmax.shape[1])
    ndit = np.nditer([full_roots, kmax], flags=['multi_index'], op_flags=['readwrite'])
    while not ndit.finished: #loop untill the end of the nditerator
        ind0, ind1 = ndit.multi_index # the the two indeces
        ndit[0] = str(ndit[0]).format(float(ndit[1])) #substnditute {} with kmax
        ndit.iternext()

    return full_roots
#end def make_full_root(froots, kmax):

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
    ndit = np.nditer([fchain, fparamnames, correction], flags=['multi_index'],
            op_flags=['readwrite'])
    while not ndit.finished: #loop untill the end of the nditerator
        try:  #if the files exists no problem opening them
            with open(ndit[0], 'r'), open(ndit[1], 'r'):
                pass
        except IOError: #if any of them does note exist set -999
            ndit[0] = '-999'
            ndit[1] = '-999'
            ndit[2] = -999
        ndit.iternext()
    return fchain, fparamnames, correction
#end def skip(fchain, fparamnames, correction):

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
    key_list: list
        list of keys of the following dictionaries: short names from the
        parameter files
    labels: dict
        dictionary with key: short name; value: latex string
    mean, stddev: dict
        dictionaries with key: short name; value: ndarray of paramnames.shape
        containing the mean and standard deviations
    """
    key_list=[] #to keep it ordered the dict keys ordered
    # these dictionaries will have as key the parameter short names
    labels = {} #list of y labels. key: short name; value latex code within $$
    # each value is a nd array of float of paramnames.shape 
    mean, stddev = {}, {}
    value_prototipe = np.empty_like(paramname, dtype=float)
    value_prototipe.fill(np.nan)

    # iterate over the parameter and the chain file names
    ndit = np.nditer([paramname, chains], flags=['multi_index'],
            op_flags=['readonly'])
    while not ndit.finished: #loop untill the end of the iterator
        ind0, ind1 = ndit.multi_index # the two indeces
        try:  #try to read files and save the mean and stddev
            paramindex = cf.get_paramnames(str(ndit[0]), params=params, ext="",
                verbose=verbose, skip=skip)
            indeces, short_name, long_name = [], [], []
            for pi in paramindex:
                indeces.append(pi[0]+2)
                short_name.append(pi[1])
                long_name.append(pi[2])
            #insert the index of the first columns with the weights
            indeces.insert(0, 0)
            #read the chains
            chain = np.loadtxt(str(ndit[1]), usecols=indeces)

            #compute mean and average of the chains
            chain_mean = np.average(chain[:,1:], axis=0, weights=chain[:,0])
            chain_stddev = ms.stddev(chain[:,1:], weights=chain[:,0], axis=0)

            #loop over the long_names. Save the long_names in key_list and the
            #mean and std into the corresponding dictionary
            for tmean, tstd, sn, ln in zip(chain_mean, chain_stddev, short_name, long_name):
                # if the key is not already in the dictionary add the label and
                # initialise the value in mean and stddev to value_prototipe
                if sn not in key_list: 
                    key_list.append(sn)
                    labels[sn] = '$'+ln+'$'
                    mean[sn] = np.copy(value_prototipe)
                    stddev[sn] = np.copy(value_prototipe)
                #save the mean and stddev
                mean[sn][ind0, ind1] = tmean
                stddev[sn][ind0, ind1] = tstd
        #if it fails in reading the parameter file names, do just go to the next loop
        except cf.ContourError: 
            pass
        ndit.iternext()

    # if nothing found in the files
    if len(key_list) == 0:
        print("No parameter found")
        sys.exit(10)

    return key_list, labels, mean, stddev
#end def files_to_mean_std(paramname, chains, params=None, verbose=True, skip=False):
    
def multiply_dic(d, corr):
    """
    Multiply all the elements in the dictionary by the factor 'corr'
    Parameters
    ----------
    d: dict
        dictionary with value: ndarray 
    corr: number or ndarray of shape d[key].shape
        multiplicative correction
    output
    ------
    d: dict
        input dictionary corrected
    """
    for k, v in d.items():
        d[k] = v*corr
    return d

def rescale_y(rescale, mean, std, labels):
    """
    Rescale the mean and standard deviation multiplying the desired parameters
    by the desired factor. Modify accordingly the labels 
    Parameters
    ----------
    rescale: list of lists
        list containing two element lists of short paramname names and rescaling
    mean, std: dict
        dictionaries with key: short name; value: ndarray of paramnames.shape
        containing the mean and standard deviations
    labels: dict
        dictionary with key: short name; value: latex string
    output
    ------
    modified mean, std and labels
    """
    for (k, r) in rescale:
        if k not in mean:
            warn("Key '{}' is not in mean and std dictionaries".format(k), MyWarning)
        else:
            r = float(r)
            mean[k] = mean[k]*r
            std[k] = std[k]*r
            if(r.is_integer()==True):
                strr = str(int(r))
            else:
                strr = str(r)
            labels[k] = "$"+strr+"$"+labels[k] 
    return mean, std, labels

def set_mpl_defaults(fsize, leg_fsize, lw, ms):
    """
    Set the defaults rc parameters for matplotlib
    Parameters
    ----------
    fsize: size string or font size in point
        font size
    leg_fsize: size string or font size in point
        legend font size
    lw: number
        lines width
    ms: number
        marker size
    """
    import matplotlib as mpl

    mpl.rcParams['font.size'] = fsize
    mpl.rcParams['axes.labelsize'] = fsize
    if leg_fsize is None:
        mpl.rcParams['legend.fontsize'] = fsize
    else:
        mpl.rcParams['legend.fontsize'] = leg_fsize
    mpl.rcParams['lines.linewidth'] = lw
    mpl.rcParams['lines.markersize' ] = ms

def make_figure(key_list, xlabel, ylabels, figsize=None):
    """
    Create the figure with len(key_list) subplots. Save the subplots into a
    dictionary and make the x and y axis labels
    Parameters
    ----------
    key_list: list
        list of keys of the following dictionary: short names from the
        parameter files
    labels: dict
        dictionary with key: short name; value: y label
    figsize: two floats
        custom figure size in inches
    output
    ------
    fig: matplotlib figure
    axs_dic: dic
        dictionary of subplot objects with keys from key_list
    """
    n_subplots = len(key_list)
    # window size
    if figsize is None:
        xs= 10./mpm.inc2cm
        ys= 2.6*(n_subplots+1)/mpm.inc2cm
    else:
        xs, ys = figsize

    fig, axs = plt.subplots(nrows=n_subplots, ncols=1, sharex=True,
            figsize=(xs, ys))
    if n_subplots==1:
        axs = np.array([axs])
    # move the axs to a dictionary with the elements of key_list as key. Also
    # set the ylabels and remove the x labels except in the last plot
    axs_dic={}
    for k, ax in zip(key_list, axs):
        ax.label_outer()
        ax.set_ylabel(ylabels[k])
        axs_dic[k] = ax

    #set the x label for the last plot
    axs_dic[key_list[-1]].set_xlabel(xlabel)

    return fig, axs_dic

def plot(axs_dic, kmax, mean, stddev, fill=None):
    """
    Populate the axes with plots 
    Parameters
    ----------
    axs_dic: dict
        key: short param names; value: plt.subplot
    kmax: nd_array
        values of kmax to plot
    mean: dict
        key: short param names; value: nd_array of kmax.shape with the means
    stddev: dict
        key: short param names; value: nd_array of kmax.shape with the standard deviations
    fill_between: 2 element list
    """
    # set up the alpha and assosiate the correct function to plot the errorbars 
    if fill is not None:
        alphas = np.linspace(*fill, num=kmax.shape[0])  #create the range of alpha
    else:
        alphas = np.ones(kmax.shape[0])

    for k, ax in axs_dic.items():  #loop over the axes
        #reset colors, line style and markers for every axes
        colors = it.cycle(mpm.colors)
        linestyles = it.cycle(mpm.linestyles)
        markers = it.cycle(mpm.symbols[2:-2])

        for k, m, s, c, ls, ma, alpha in zip(kmax, mean[k], stddev[k], colors,
                linestyles, markers, alphas):
            toplot = np.isfinite(m)  # if any of the elements does not exists, skip it
            if fill is None:
                ax.errorbar(k[toplot], m[toplot], yerr=s[toplot], c=c, ls=ls,
                        marker=ma, alpha=alpha, label='temp')
            else:
                from errorfill import errorfill
                errorfill(k[toplot], m[toplot], yerr=s[toplot], c=c, ls=ls,
                        alpha_fill=alpha, ax=ax, label='temp')

    #set the x limit a bit larger than the plotted ones
    ax.set_xlim(kmax.min()*0.95, kmax.max()*1.04)
#end def plot(axs_dic, kmax, mean, stddev, fill=None):
def plot_variance(axs_dic, kmax, stddev):
    """
    Populate the axes with plots with only the variance
    Parameters
    ----------
    axs_dic: dict
        key: short param names; value: plt.subplot
    kmax: nd_array
        values of kmax to plot
    stddev: dict
        key: short param names; value: nd_array of kmax.shape with the standard deviations
    """
    for k, ax in axs_dic.items():  #loop over the axes
        #reset colors, line style and markers for every axes
        colors = it.cycle(mpm.colors)
        linestyles = it.cycle(mpm.linestyles)
        markers = it.cycle(mpm.symbols[2:-2])

        for k, s, c, ls, ma in zip(kmax, stddev[k], colors,
                linestyles, markers):
            toplot = np.isfinite(s)  # if any of the elements does not exists, skip it
            ax.plot(k[toplot], s[toplot], c=c, ls=ls, marker=ma, label='temp')

    #set the x limit a bit larger than the plotted ones
    ax.set_xlim(kmax.min()*0.95, kmax.max()*1.04)
#end def plot_variance(axs_dic, kmax, stddev, fill=None):

def shave_y_tick_labels(axs_dic):
    """
    automatically removes the first and last tick labels to avoid overlap and
    if ther are more than 5 labels, print only one out of two
    """
    return
    for ax in axs_dic.values():
        ymin,ymax = ax.get_ylim()
        # get y tick labels and locks
        ytl = ax.get_ymajorticklabels()
        tl = ax.yaxis.get_majorticklocs()
        # make invisible the labels outside the plotted limits
        ytl[(tl<ymin).sum()].set_visible(False)
        ytl[-(tl>ymax).sum()-1].set_visible(False)
        if(len(ytl)-2>5):  # if necessary make invisible one every second
            for iy in range(2,len(ytl)-1,2):
                ytl[iy].set_visible(False)

def main(argv):
    """
    """
    args = parse(argv)

    # command line parameter checks
    if args.legend is not None and len(args.legend)!=len(args.froot):
        print("The number of legend tags must be the same as the number of file roots")
        exit()

    if args.verbose:
        print("Computing kmax and correction")
    if args.fill is not None:  #ignore shift 
        args.shift = None
    k_max, corr = kmax_correction(args.range, args.stride, len(args.froot),
            shift=args.shift, ini=args.ini)

    if args.verbose:
        print("Create the the file names")
    args.froot = make_full_root(args.froot, k_max)

    #create the parameter and chain file names
    #substitute with np.add when will become available
    paramnames = np.core.defchararray.add(args.froot, '.'+args.ext_parmn)
    chains = np.core.defchararray.add(args.froot, '.'+args.ext_chain)

    if args.paramlist:
        print_paramlist(paramnames[0,0])

    if args.skip:
        if args.verbose:
            print("Skipping non existing files")
        chains, paramnames, corr = skip(chains, paramnames, corr)

    if args.verbose:
        print("Read the parameter name and chain files and compute the mean and stddev")
    short_names, labels_dic, mean_dic, std_dic = files_to_mean_std(paramnames,
            chains, params=args.columns, verbose=args.verbose, skip=args.skip)

    if (corr!=1).any(): # if any of the correction is different from 1 apply it
        if args.verbose:
            print("Correcting the errors")
        multiply_dic(std_dic, corr)

    if args.nsigma!=1:
        if args.verbose:
            print("{}sigma wanted in the plot".format(args.nsigma))
        multiply_dic(std_dic, args.nsigma)

    if args.rescale is not None:
        if args.verbose:
            print("rescaling the mean and stddev in the desired subplots")
        mean_dic, std_dic, labels_dic = rescale_y(args.rescale, mean_dic,
                std_dic, labels_dic)

    if args.y_label is not None:
        if args.verbose:
            print("Substituting the y labels in the desired subplots")
        labels_dic = subsitute_labels(args.y_label, labels_dic)

    if args.verbose:
        print("Set up the plot area")
    set_mpl_defaults(args.font_size, args.legend_fsize, args.line_width,
            args.marker_size)
    if args.set_MNRAS:  # use MNRAS presets
        mpm.MNRAS_fig()

    fig, daxs = make_figure(short_names, args.x_label, labels_dic, figsize=args.figsize)

    if args.verbose:
        print("Plot mean and stddev")
    if args.variance_only:
        plot_variance(daxs, k_max, std_dic)
    else:
        plot(daxs, k_max, mean_dic, std_dic, fill=args.fill)

    if args.horizontal is not None:
        mpm.plot_horiz_vert(daxs, horiz=args.horizontal)

    if args.y_range is not None:
        mpm.change_xylim(daxs, y_lim=args.y_range)

    shave_y_tick_labels(daxs)

    # legend
    if args.legend is not None:
        if args.legend_plot is None:
            args.legend_plot = short_names[0]
        legend_ = mpm.draw_legend(fig, daxs, args.legend, args.legend_plot, args.loc)

    plt.tight_layout(h_pad=0, pad=0.2, rect=args.bounds)
    if args.ofile is None:
        plt.show()
    else:
        fig.savefig(args.ofile)

if __name__ == "__main__":   # if is the main

    main(sys.argv[1:])
    exit()
    



