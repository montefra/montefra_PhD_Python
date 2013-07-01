#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Rebin the power spectrum in the input files
"""

import my_functions as mf
import numpy as np

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

    description = """    Rebin the power spectra for the input files.
    The structure of the input files must be like the following:
        # 	 sum(w) 	 sum(w^2n(z)) 	 sum(w^2)
        #data	 ###             ###             ###   
        #random	 ###             ###             ###   
        #	k	P(k)    	P(k)+noise      n_modes
        [...]
    The output files will have the same structure
    If the number of bins in the input file are not multiples of 'n_bins',
    the last bins will be discarded
    """
    p = ap.ArgumentParser(description=description, formatter_class=apc.RawDescrArgDefHelpFormatter)

    p.add_argument('merge_bins', action='store', type=int, 
            help="Number of bins to merge")
    p.add_argument("ifname", action=apc.file_exists(), nargs='+', 
            help="""Input file name(s).""")

    p = apc.version_verbose(p, '1')

    p, group = apc.insert_or_replace(p)
    p, group = apc.overwrite_or_skip(p)

    p.add_argument('-f', '--from-bin', type=int, action='store',
            help='Ignore the first %(dest)s bins')

    p.add_argument("--fmt", default="%7.6e", action=apc.StoreFmt, nargs='+',
            help="Format of the output files")

    return p.parse_args(args=argv)
#end def parse(argv):

def read_file(fname):
    """
    read the file 'fname' saving the full header and the power spectrum table
    Parameters
    ----------
    fname: string
        file name
    output
    ------
    full_header: list
        full header from the file
    pk: ndarray
        power spectrum table
    """
    #read the file name and the header
    with open(fname, 'r') as f:
        header_lines = mf.n_lines_comments(f) #get the size of the header
        full_header = [f.readline() for i in xrange(header_lines)] #save the header as it is
        pk = np.loadtxt(f)  #read the power spectrum
    return full_header, pk

def rebin(fname, mbins, **kwargs):
    """
    Rebin the power spectrum in file 'fn' and save it into the output file
    Parameters
    ----------
    fname: string
        file name
    mbins: int
        number of bins to merge to create the input

    accepted kwargs that affects the function:
    +from_bin: skip the first 'from' bins
    +replace: replace string *replace[0]* with *replace[1]* in f.name
    +insert: insert string *insert[0]* before *insert[1]* in f.name
    +skip: existing file names skipped [True|False]
    +overwrite: existing file names overwritten [True|False]
    +fmt: format of the output file
    """

    ofile = mf.create_ofile_name(fname, **kwargs) # create the output file name
    if ofile == None:
        return None

    header_lines, ps = read_file(fname) #read the file

    # if required to skip some line
    fbin = kwargs.get('from_bin')
    if fbin is not None: 
        ps = ps[fbin:,:]
    # check that the number of remaining bins is multiple of mbins and if not
    # cut the last bins
    pssize = ps.shape[0]
    psremaining = pssize%mbins
    if psremaining!=0:
        ps = ps[:-psremaining,:]
    # separate the power spectrum into the various components
    k, pk, pksn, nmodes = ps.T  
    # reshape to simplify rebinning
    finalbins = ps.shape[0]/mbins
    kr = k.reshape([finalbins, mbins])
    pkr = pk.reshape([finalbins, mbins])
    pksnr = pksn.reshape([finalbins, mbins])
    nmodesr = nmodes.reshape([finalbins, mbins])

    # get the values for the rebinned power spectrum
    kr = np.average(kr, axis=1)
    pkr = np.average(pkr, axis=1, weights=nmodesr)
    pksnr = np.average(pksnr, axis=1, weights=nmodesr)
    nmodesr = nmodesr.sum(axis=1)

    # save the rebinned power spectrum
    with open(ofile, 'w') as of:
        of.writelines(header_lines)
        np.savetxt(of, np.vstack([kr, pkr, pksnr, nmodesr]).T, delimiter='\t',
                fmt=kwargs['fmt'])
#end def rebin(fname, mbins, **kwargs):

if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])

    for fn in args.ifname:
        if args.verbose:
            print("Process file '{}'".format(fn))
        rebin(fn, args.merge_bins, **vars(args))



    exit()
