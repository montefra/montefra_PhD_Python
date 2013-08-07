# -*- coding: utf-8 -*-
"""It contains functions to prepare the contour plots from the 
MCMC chains produced by COSMOMC"""

import numpy as np  #numpy
import smooth as sm  #import the gaussian kernel density estimate function 
import sys

class ContourError(Exception):
    """exception in this module and Plot/contour_plot.py"""
    def __init__(self, string):
        self.string=string
    def __str__(self):
        return repr(self.string)

def _check_type(string_iter):
    """
    Check if froot is a string or an iterable with a 'append' method (in this order)
    Returns boolean if it is a string or like a list
    Parameters
    ----------
    froots: list, string
        list of file root of the chain to be plotted
    output
    ------
    is_string: bool
    is_listlike: bool
    """
    is_string = False
    is_listlike = False

    import six
    if isinstance(string_iter, six.string_types):
        is_string=True
    elif hasattr(string_iter, 'append'):
        is_listlike = True
    else:
        raise TypeError("Only string and list-like object accepted")

    return is_string, is_listlike

def _get_paramnames(fname, params=None, verbose=False, skip=False):
    """ 
    Read the parameter file name 'fname'. 
    If fname does not exist and skip==True return None
    Parameters
    ----------
    fname: string
        file name
    params: list (optional)
        name of the parameters of interest as in the first column of the
        parameter name files. If None all parameters saved
    verbose: bool (optional)
        more output
    skip: bool (optional)
        skip non existing files
    output
    ------
    parindex: list of lists
        parindex[i] = [index_in_file, first_column_string, second_column_string]
        with len(parindex)==len(params) or the whole file
    """
    try:
        with open(fname, 'r') as f: # open each file
            if(verbose == True):
                print("Reading file {}".format(f.name))
            temp = f.readlines()   # and read it all
    except IOError as e:
        if skip:
            return None
        else:
            raise IOError(e)

    #save the paramnames and theid index into a list of lists
    parindex = []
    if params is None:
        for i, mn in enumerate(temp):
            tmn = [m.strip() for m in mn.split('\t')]
            parindex.append([i, tmn[0], tmn[1]]) 
    else:
        for p in params:  #loop the parameters to be used
            for i, mn in enumerate(temp): #loop the lines of the parameter file
                # split each line to separate the small name and the latex name
                tmn = [m.strip() for m in mn.split('\t')]
                if(p==tmn[0].strip('*')):  # shave space and '*'
                    # save the index and the latex name. 
                    parindex.append([i, tmn[0], tmn[1]]) 
                    break
        if len(parindex) != len(params):
            ContourError("""File '{}' doesn't have at least one of the
                    parameters '{}'""".format(fname, params))
    return parindex

def get_paramnames(file_roots, params=None, ext=".paramnames", verbose=False, skip=False):
    """For each file 'file_roots+ext' finds the indeces of the elements in params
    that matcht what is stored in the first columns of the parametere file names

    Parameters
    ----------
    file_roots: string, list-like (must have append method)
        list of file root of the chain to be plotted
    params: list (optional)
        name of the parameters of interest as in the first column of the
        parameter name files. If None all parameters saved
    ext: string (optional)
        extention of the parameters file
    verbose: bool (optional)
        more output
    skip: bool (optional)
        skip non existing files
    output
    ------
    parindex: 
        if file_roots is string: len(parindex) == len(params) or the whole file
        if file_roots is list-like: parindex is a a list len(file_roots) 
        In the last case len(parindex[i])==len(params) or the whole file.
        In any case:
        parindex[i][j] = [index_in_file, first_column_string, second_column_string]
    """
    #check the type of the input
    is_string, is_listlike = _check_type(file_roots)

    if is_string: 
        parindex = _get_paramnames(file_roots+ext, params=params,
                verbose=verbose, skip=skip)
        if parindex is None:
            ContourError("The file '{}' does not exists".format(file_roots+ext))
        else:
            return parindex
    elif is_listlike:
        parindex = [_get_paramnames(fr+ext, params=params, verbose=verbose,
            skip=skip) for fr in file_roots]
        return(parindex)


def get_chains(file_roots, cols, ext=".txt", verbose=False):
    """Reads the paramnames files of the cosmomc

    Parameters
    ----------
    file_roots: list
    list of file root of the chain to be plotted
    cols: list of lists with len(cols)==len(file_roots)
    list containing the list of columns to read from the files
    ext: string (optional)
    extention of the chain files
    verbose: bool (optional)
    if True print more output
    output
    ------
    chains: list of array
    list of array containing all the input chains
    """

    chains = []   #save the arrays of the chains
    for fr, c in zip(file_roots, cols):     #loop through the total chains
        c.insert(0,0) #insert the index of the first column (with the weights)
        with open(fr+ext, 'r') as f:
            if verbose:
                print("Reading file '{}'".format(f.name))
            chains.append(np.loadtxt(f, usecols=c))   #load the chain file
    return chains


def hist2D(a, xr=None, yr=None, bins=30, smooth=False):
    """
    Given a list of arrays a containing len(a) indipendent chains returns a
    list of array containing all the 2D histograms 
    Smoothing is applied if required

    Parameters
    ----------
    a: list of arrays
        each element must have 2 or 3 columns. If has 3 columns the first one
        will be considered as containing the weights.
    xr: list (optional)
        x_range of the histogram: if not given will be determined from 'a'
        itself
    yr: list (optional)
        y_range of the histogram
    bins: int (optional)
        number of bins for the 2D histogram. This numer will be doubled when
        smothing
    smooth: bool (optional)
        if False numpy.histogram2D used, otherwise fast_kde (gaussian kernel
        density estimate)

    output: list of arrays, list of arrays, list of arrays
        hist2D: list of arrays containing all the 2D histograms
        xe: list of arrays containing the bin means along the first dimension
        ye: list of arrays containing the bin means along the second dimension
    """

    h2D, xe, ye = [], [], []  #lists of outputfiles
    for c in a:  #loop over the input list of arrays
        if(c.shape[1] == 2):  #if a has 2 columns add 1 in the first one
            c = insert(c, 0, 1., axis=1)
        elif(c.shape[1] == 3):  #if a has 3 columns do nothing
            pass
        else:
            print("Too many or too few columns")
            sys.exit(20)

        txr, tyr = xr, yr
        if(xr == None):  #if not given get the x and y ranges
            txr = [np.amin(c[:,1]), np.amax(c[:,1])]
        if(yr == None):
            tyr = [np.amin(c[:,2]), np.amax(c[:,2])]

        if smooth:
            h2D.append(np.flipud(sm.fast_kde(c[:,1], c[:,2],
                gridsize=(bins,bins), extents=(txr[0],txr[1],tyr[0],tyr[1]),
                weights=c[:,0])))  #create a smooth 2D histogram
            txe, xst = np.linspace(txr[0], txr[1], num=bins, endpoint=False,
                    retstep=True)
            tye, yst = np.linspace(tyr[0], tyr[1], num=bins, endpoint=False,
                    retstep=True)
            txe, tye = txe+0.5*xst, tye+0.5*yst
            #arrays containing the mean of x and y axis of the 2d histogram
            xe.append(txe)
            ye.append(tye)
        else:
            th2D, txe, tye = np.histogram2d(c[:,1], c[:,2], bins=bins,
                    range=[txr, tyr], weights=c[:,0])  #create the 2D histogram
            h2D.append(th2D.T)  #append the histrogram
            txe, tye = (txe[:-1] + txe[1:])/2., (tye[:-1] + tye[1:])/2.
            xe.append(txe)
            ye.append(tye)
    return h2D, xe, ye


def h2D2conflev(a, levels):
    """
    This function gets a list of 2D histograms and of confidence levels
    and returns an a list of array containing the amplitudes (relative to the maxima)
    of the histogram that correspond to the given levels

    a: list if numpy arrays
        list of array containing 2D histograms.
    levels: list
        required confidence levels

    output: list of 1D arrays
        amplitude, relative to the maximum, of the 2D histograms corresponding
        to the confidence levels 'levels'
    """

    levels.sort()   #sort so there's no problem
    n_levels = len(levels)   #number of levels

    output = []  #create a list containing all the 
    for c in a:
        [biny, binx] = c.shape
        toutput = np.zeros(n_levels+1)   #allocate the temporary output
        sc = np.sum(c)  #sum of all the elements into the 2d histogram  
        toutput[0] = np.amax(c)
        tlev = sc*levels
        #flatten the histogram, sort its arguments and invert the order
        sortc = np.argsort(c.flatten(), kind='mergesort')[::-1]

        tot = 0.  #initialize the sum
        l = 0  #integer to go through the levels
        for ind in sortc:  #go through the indeces
            tc=c[ind/binx,ind%binx]
            #the sum of the amplitude starting from mimimum
            tot += tc    
            if(tot == tlev[l]):
                toutput[l+1] = tc
                l += 1
            elif(tot > tlev[l]):
                toutput[l+1] = ((tc-tco)*tlev[l] + tco*tot - tc*(tot-tc)) / tc
                l += 1
            if(l >= n_levels):
                break
            #save the value of amplitude of the histogram for the next loop if
            #interpolation needed
            tco=tc

        output.append(toutput[::-1])   #append the levels of the 2D histogram
    return output

