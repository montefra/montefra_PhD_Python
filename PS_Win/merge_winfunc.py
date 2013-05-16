#!/usr/bin/python
# -*- coding: utf-8 -*-
#merge window functions in a unique file 

import numpy as np
import smooth

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

    description = """ 
    Merge the input window matrices or power spectra files.
    The first two columns must be k and W(k) or P(k). The files are ordered
    according to k[0] and all k>k[-1]*fkN are ignored.  When two or more window
    functions/power spectra overlap the one with smaller k[0] will be
    considered.  
    """

    p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument("ofname", action='store', help="Output file name")

    p.add_argument("ifnames", action=apc.file_exists(), nargs='+', 
            help="Input file name(s)")

    p = apc.version_verbose(p, '0.1')
    
    p.add_argument("-o", "--overwrite", action="store_true", 
            help="Overwrite the output file name")

    p.add_argument("--fraction-kN", action="store", type=float, 
            default=0.65, dest='fkN',
            help="All the modes larger than fkN*k[-1] are discarded.")

    smooth = p.add_argument_group(title='Smoothing', 
            description="""Smooth (or fit) the merged files""")
    excsmooth = smooth.add_mutually_exclusive_group()
    
    excsmooth.add_argument('-m', '--mean', action='store_true', 
            help="""Subsitute the window function with its mean""")

    excsmooth.add_argument('-s', '--smooth', action='store_true',
            help="""Smooth the window function""")

    excsmooth.add_argument('-f', '--fit', action='store_true',
            help="""Perfom a power law fit""")

    smooth.add_argument('-k', '-k-sub', action='store', type=float, nargs=2,
            default=[-1,-1], help="""Minimum and maximum wavenumbers to
            substitute. If t""")

    smooth.add_argument('--k-fit', action='store', type=float, nargs=2,
            help="""Minimum and maximum wavenumbers to use to compute the mean
            or do the fit. If not given, use 'k_sub'. Ignored in 'smooth'.""")

    wchoises = ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
    smooth.add_argument('-w', '-window', action='store', choices=wchoises,
            default='hanning', help="Only for 'smooth'. Window type")
    smooth.add_argument('-l', '--window-len', action='store', default=11, type=int, 
            help="Only for 'smooth'. Dimension of the smoothing window; should be an odd integer")

    return p.parse_args(args=argv)

def read_order(ifnames):
    """
    Read and order the file in increasing order of k[0]
    Parameters
    ----------
    ifnames: list of strings
        file names
    output
    ------
    win: list of 2D numpy arrays
        ordered window functions or power spectra
    """
    win = [np.loadtxt(fn) for fn in ifnames] # read the files
    win.sort(key=lambda w: w[0,0])  # sort in-place according to the value of k[0]=w[0,0]
    return win

def unify_win(wins, fkN=0.65):
    """
    Unify the window functions. The windows must be ordered from with
    increasing k = wins[:][0,0] For each window function all rows where
    wins[:][0,:] > fkN*wins[i][0,-1] will be discarder. When two window
    functions overlap the one on larger scale will be considered.

    Parameters
    ----------
    wins: list of 2D arrays
        list of window functions. Each element with at least 2 columns: k, W(k)
    fkN: float (optional)
        fraction of the maximum wave number to discard

    output
    ------
    merged: 2D array
        array with the same number of columns of the input ones, containing the merged window function
    """
    k_min = 0.
    win = []   #output total window
    for w in wins:
        w[:,1] -= noise   #subtract the shot noise
        k = w[:,0]    #get the k
        imin, imax = np.sum(k<k_min), np.sum(k<(fkN*k[-1]))   #get the first and the last bins
        win.extend(w[imin:imax,:])    #add the part of interest of the window function 
        k_min = k[imax]  #substitute the k_min for the next loop with the largest k of this loop

    return np.array(win)  #return the full window

if __name__ == "__main__":   # if is the main

    import sys
    args = parse(sys.argv[1:])

    import io_custom as ioc
    # check output file
    if args.overwrite and ioc.file_exists(args.ofname):
        print("The output file already exists")
        exit()

    win = read_order(args.ifnames) # read the files


    win = unify_win(wins, fkN=opt.fkN)   #unify window functions

    if(opt.kav != None):  #substitute every value for k>kav with the mean of the window for those k
        kindex = (win[:,0]>opt.kav).sum()
        if(opt.smooth == None):
            win[-kindex:,1] = np.mean(win[-kindex:,1])
        else:
            print("understand how to implement it properly")
            win[-kindex:,1] = smooth.smooth1D(win[-kindex:,1], window=opt.smooth[0], window_len=int(opt.smooth[1]))[:kindex]

    np.savetxt(args[0], win, fmt='%7.6e', delimiter='\t')

    exit()

