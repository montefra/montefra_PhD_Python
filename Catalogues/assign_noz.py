#!/usr/bin/python
# -*- coding: utf-8 -*-
#Assign n(z) read from a file

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

    description = """Given a table of 'z,n(z)', subsitute a column in the
    catalogue(s) with the value of n(z), using the values of z from an other
    column of the same file """

    p = ap.ArgumentParser(description=description,
            formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument("noz", action="store", type=ap.FileType('r'), help="File
            containing a table of z, n(z)")

    p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'),
            help="Input file name(s), containing z in one of the columns")

    p = apc.version_verbose( p, '1' )

    p, group = apc.insert_or_replace1(p)
    p, group = apc.overwrite_or_skip(p)

    p.add_argument("-z", "--z-column", action="store", type=int, default="-1",
            help="Column in the catalogue with the redshift.")
    p.add_argument("-n", "--nz-column", action="store", type=int, default="5",
            help="Column in the catalogue in which n(z) is saved.")
    p.add_argument("--zn-file-col", action="store", type=int, nargs=2,
            default=[0,1], help="""Columns in 'noz' containing z and n(z).""")

    p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+',
            help="Format of the output files")

    description = """Parameters related to the parallel computation"""
    p, parallel = apc.parallel_group( p, description=description )

    return p.parse_args(args=argv)  
#end def parse(argv)

def assign_noz( f, nofz, **kwargs):
    """read file 'f', substitute a columns with noz(z), with z in 'f' itself, and
    save in a file.
    Parameters
    ----------
    f: file object or string
        file containing the catalogue
    nofz: function
        returns n(z) for all the z in the input file
    output
    ------
    none

    accepted kwargs that affects the function
    +verbose: verbose mode [True|False] 
    +replace: replace string *replace[0]* with *replace[1]* in f.name
    +insert: insert string *insert[0]* before *insert[1]* in f.name
    +skip: existing file names skipped [True|False]
    +overwrite: existing file names overwritten [True|False]
    +fmt: format of the output file
    """
    if( type(f) == file ):  #if f is a file object
        fname = f.name  #get the file name
    else:  #it's alread the file name
        fname = f

    if(kwargs['verbose'] == True):
        print("Process catalogue '{0}'.".format(fname))

    #create the output file name and check it
    if(kwargs['replace'] == None):
        ofile, skip = mf.insert_src(fname, kwargs['insert'],
            overwrite=kwargs['overwrite'], skip=kwargs['skip'])
    else:
        ofile, skip = mf.replace_src(fname, kwargs['replace'],
            overwrite=kwargs['overwrite'], skip=kwargs['skip'])
    if(skip == True):
        print("Skipping file '{0}'".format(fname))
    return None

    cat = np.loadtxt( f )  #read the input catalogue

    cat[:,kwargs['nz_column']] = nofz( cat[:,kwargs['z_column']] )

    np.savetxt(ofile, cat, fmt=kwargs['fmt'], delimiter='\t')
    return None
#end assign_noz( f, noz, **kwargs)

if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])

    #read n(z)
    z, noz = np.loadtxt(args.noz, usecols=args.zn_file_col ).T
    #create an interpolation function
    import scipy.interpolate as spi
    #interp_noz = lambda x: np.interp( x, z, noz, left=0, right=0)
    interp_noz = spi.InterpolatedUnivariateSpline( z, noz ) #, fill_value=0, bounds_error=False ) 

    #if parallel computation required, check that Ipython.parallel.Client 
    #is in installed and that the ipycluster has been started
    if args.parallel :
        import ipython_parallel as IPp
        #command to run on all the engines
        imports = [ 'import numpy as np', 'import my_functions as mf', ]
        args.parallel, lview = IPp.start_load_balanced_view( to_execute=imports )

    #loop through the catalogues and add a columns with n(z)
    if( args.parallel == False ):  #if: parallel
        for fn in args.ifname:  #file name loop
            #substitute n(z)
            assign_noz( fn, interp_noz, **vars(args) ) 
    #run the script using the IPython parallel environment 
    else:    #if: parallel
        engines_id = lview.client.ids  #get the id of the engines_id
        initstatus = lview.queue_status()  #get the initial status

        #submit the jobs and save the list of jobs
        import os
        runs = [ lview.apply( assign_noz, os.path.abspath(fn.name), interp_noz,
            **vars(args) ) for fn in args.ifname ]

        if args.verbose :   #if some info is required
            IPp.advancement_jobs( lview, runs, engines_id, update=args.update,
                    init_status=initstatus )
        else:   #if no info at all is wanted
            lview.wait( jobs=runs )  #wait for the end

        #just check for any error
        results = [r.result for r in runs]
    #end if: parallel

    exit()
