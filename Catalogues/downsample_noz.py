#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Given a file with the fraction for downsample as function of redshift and one or more
catalogues, downsample them according to the given fractions"""

import my_functions as mf
import numpy as np
import pandas as pd
import scipy.interpolate as spip


interp_type = ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'] 
def select_interp(string):
    """Check if string is an integer or one of these strings
    'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic' and return it with the correct type.
    Otherwise throw an error.
    
    Parameters
    ----------
    string: string
        input string to check
        
    output
    ------
    value: int or string
    """

    try:
        value = int(string)  #check if it is an int
    except ValueError:   #if not check if it's a string among the possible choises
        if( string in interp_type):
            value = string
        else:
            msg = "'{}' is neither an integer nor a string among those: {}".format(string, stchoise)
            raise ap.ArgumentTypeError(msg)
    return(value)

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
    
    description = """Given a file with the fraction for downsample as function of redshift and 
            one or more catalogues, downsample them according to the given fractions"""
    p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument("downsample", action=apc.file_exists(),
            help="""File with the fraction to downsample as function of redshift.""")
    p.add_argument("ifname", action=apc.file_exists(), nargs='+', 
            help="Input file name(s), containing z in one of the columns")

    p = apc.version_verbose(p, '0.1')

    p.add_argument("-z", "--z-column", action="store", type=int, default=2,
            help="Column containing the redshift")

    p.add_argument("-k", "--kind", action="store", type=select_interp, default="linear",
            help="""Specifies the kind of interpolation of the downsample file
            as a string {} or as an integer specifying the order of the spline
            interpolator to use. Default is ‘linear’.""".format(interp_type))

    p, group = apc.insert_or_replace(p)
    p, group = apc.overwrite_or_skip(p)

    p, pandas = apc.pandas_group(p)

    p.add_argument("--fmt", default="%7.6e", action=apc.StoreFmt, nargs='+', help="Format of the output files")

    return p.parse_args(args=argv)

def downsample_noz(f, interpf, zcol, **kwargs):
    """
    read the input files, downsample it and save it
    Parameters
    ----------
    f: file object or string
        file containing the catalogue
    interpf: callable
        function that returns the fraction used for the downsample
    zcol: int
        column in the input file to use for the downsampling
    output
    ------

    accepted kwargs that affects the function
    +substitute: substitute the content of the columns involved with operations
        with this value, if not None
    +verbose: verbose mode [True|False] 
    +replace: replace string *replace[0]* with *replace[1]* in f.name
    +insert: insert string *insert[0]* before *insert[1]* in f.name
    +skip: existing file names skipped [True|False]
    +overwrite: existing file names overwritten [True|False]
    +pandas: use pandas for the input
    +chunks: chunksize in pandas.read_table
    +fmt: format of the output file
    """
    ofile = mf.create_ofile_name(f, **kwargs) # create the output file name
    if kwargs['verbose']:
        print("Processing file '{}'".format(f))

    # read the input catalogue
    if kwargs['pandas']:
        if kwargs['chunks'] is None:
            cat = pd.read_table(f, header=None,
                    skiprows=mf.n_lines_comments(f), sep='\s')
            to_keep = cat[zcol]  # save the column to use for the selection 
            to_keep = np.random.rand(to_keep.size) <= interpf(to_keep)  #random select objects to keep
            np.savetxt(ofile, cat[to_keep], fmt=kwargs['fmt'], delimiter='\t')
        else:
            chunks = pd.read_table(f, header=None, sep='\s',
                    skiprows=mf.n_lines_comments(f), chunksize=kwargs['chunks'])
            with open(ofile, 'w') as of:
                for cat in chunks:
                    to_keep = cat[zcol]  # save the column to use for the selection 
                    to_keep = np.random.rand(to_keep.size) <= interpf(to_keep)  #random select objects to keep
                    np.savetxt(of, cat[to_keep], fmt=kwargs['fmt'], delimiter='\t')
    else:
        cat = np.loadtxt(f)
        to_keep = cat[:, zcol]  # save the column to use for the selection 
        to_keep = np.random.rand(to_keep.size) <= interpf(to_keep)  #random select objects to keep
        np.savetxt(ofile, cat[to_keep, :], fmt=kwargs['fmt'], delimiter='\t')
# end def downsample_noz(f, interpf, zcol, **kwargs):

def interp_func(fname, kind):
    """
    Read the file containing the fractions to use to downsample and 
    return a 1d interpolation functions
    Parameters
    ----------
    fname: string
        file name
    kind: int or string
        interpolation kind
        
    output
    ------
    func: callable
        function that performs the interpolation
    """
    #read the fraction of according to which the downsample has to be done
    z, fdown = np.loadtxt(args.downsample).T  
    #function with the interpolation. 0 is returned outside x and y range
    func = spip.interp1d(z, fdown, kind=kind, fill_value=0., bounds_error=False)
    return func

if __name__ == "__main__":   #if it's the main

    import sys  #system stuff
    args = parse(sys.argv[1:])

    if args.verbose:
        print("Reading the file for the downsampling")
    interpf = interp_func(args.downsample, args.kind)

    for fn in args.ifname:  
        downsample_noz(fn, interpf, args.z_column, **vars(args))

    exit()

        


