#!/usr/bin/python
# -*- coding: utf-8 -*-
"""rescale the x, y and z positions of object in a catalogue by a given amount """

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

    description = """Subtract from the first three columns of the input files the
    three floats. If they are not provided they will be computed from the files
    themselves. """

    p = ap.ArgumentParser(description=description,
            formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'),
            help="Input file name(s), containing z in one of the columns")

    p = apc.version_verbose( p, '1' )

    p, group = apc.insert_or_replace1(p)
    p, group = apc.overwrite_or_skip(p)

    p.add_argument("--rescale", action="store", type=float, nargs=3,
            help="""Values to subtract from the first three columns""")
    p.add_argument("-o", "--offset", action="store", type=float, default=0,
            help="""Extra offset to add to 'rescale'""")

    p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+',
            help="Format of the output files")

    description = """Parameters related to the parallel computation"""
    p, parallel = apc.parallel_group( p, description=description )

    return p.parse_args(args=argv)  
#end def parse(argv)


def get_minimum( f ):
    """Parameters
    ----------
    f: file object or string
        file containing the catalogue
    output
    ------
    mincat: 3 element array

    """
    if( type(f) == file ):  #if f is a file object
        fname = f.name  #get the file name
    else:  #it's alread the file name
        fname = f

    return np.loadtxt( f, usecols=[0,1,2] ).min( axis=0 )  
#end def get_minimum

def subtract_fromfile( f, absmin, **kwargs):
    """read file 'f', subtract absmin+offset from the first three columns and
    save the file
    Parameters
    ----------
    f: file object or string
        file containing the catalogue
    absmin: list or array of three floats
        values to subtract from the first 3 columns
    output
    ------
    none

    accepted kwargs that affects the function
    +verbose: verbose mode [True|False] 
    +replace: replace string *replace[0]* with *replace[1]* in f.name
    +insert: insert string *insert[0]* before *insert[1]* in f.name
    +skip: existing file names skipped [True|False]
    +overwrite: existing file names overwritten [True|False]
    +offset: extra offset to add to absmin
    +fmt: format of the output file
    """
    ofile = mf.create_ofile_name(f, **kwargs) # create the output file name

    cat = np.loadtxt( f )  #read the input catalogue
    cat[:,:3] -= absmin+kwargs['offset']  #subtract the minimum

    np.savetxt(ofile, cat, fmt=kwargs['fmt'], delimiter='\t')
    return None
#end def subtract_fromfile

if __name__ == "__main__":   # if is the main

    import sys
    args = parse(sys.argv[1:])
    args.rescale = np.array( args.rescale )  #convert to numpy array rescale

    #if parallel computation required, check that Ipython.parallel.Client 
    #is in installed and that the ipycluster has been started
    if args.parallel :
        import os
        from ipython_parallel import Load_balanced_view as Lbv
        parallel_env = Lbv()  #initialize the object with all my parallen stuff
        args.parallel = parallel_env.is_parallel_enabled()

    if args.parallel :
        imports = [ 'import numpy as np', 'import my_functions as mf', ]
        parallel_env.exec_on_engine( imports )

    #search for the absolute minimum if not given
    if(args.rescale == None):
        if args.verbose : 
            print("Finding the absolute maximum and minimum")
        if( args.parallel == False ):  #if: parallel
            mincat = []  #store the minima
            for fn in args.ifname:  #file name loop
                mincat.append( get_minimum( fn ) )
        #run the script using the IPython parallel environment 
        else:    #if: parallel
            initstatus = parallel_env.get_queue_status()  #get the initial status
            #submit the jobs and save the list of jobs
            runs = [parallel_env.apply( get_minimum, os.path.abspath(fn.name) ) 
                for fn in args.ifname]

            if args.verbose :   #if some info is required
                parallel_env.advancement_jobs(runs, update=args.update,
                        init_status=initstatus)
            else:   #if no info at all is wanted
                parallel_env.wait(jobs=runs)  #wait for the end

            #get the minimum
            mincat = [r.result for r in runs]
        #endif: parallel
        args.rescale = np.min(mincat, axis=0)  #get the absolute minimum

        if args.verbose : 
            print("Absolute minimum: ", args.rescale)

    #rescale
    if args.verbose : 
        print("Subtracting the minimum")
    if( args.parallel == False ):  #if: parallel
        for fn in args.ifname:  #file name loop
            subtract_fromfile(fn, args.rescale, **vars(args))
    else:    #if: parallel
        initstatus = parallel_env.get_queue_status()  #get the initial status
        #submit the jobs and save the list of jobs
        runs = [parallel_env.apply(subtract_fromfile, os.path.abspath(fn.name),
            args.rescale, **vars(args)) for fn in args.ifname]

        if args.verbose :   #if some info is required
            parallel_env.advancement_jobs(runs, update=args.update,
                    init_status=initstatus)
        else:   #if no info at all is wanted
            parallel_env.wait(jobs=runs)  #wait for the end
        #clear the variable in the parallel environment to avoid filling up memory
        parallel_env.clear_cache()  
    #endif: parallel

    exit()
