#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Move or swap some column in a (list of) file(s)
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

    description = "Move or swap some column in a (list of) file(s)."

    p = ap.ArgumentParser(description=description, )

    p.add_argument("from", action="store", type=apc.int_or_list, 
            help="""Column(s) to copy. (if more than one, give them as a string
            of integers separated by spaces)""")
    p.add_argument("to", action="store", type=apc.int_or_list,
            help="""Where to copy the column(s). (if more than one, give them
            as a string of integers separated by spaces)""")
    p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'),
            help="""Input file name(s), containing ra and dec in the first two
            columns""")

    p = apc.version_verbose(p, '0.1')

    p.add_argument("-s", "--swap", action="store_true", 
            help="Swap 'from' with 'to'")

    p, group = apc.insert_or_replace1(p, print_def=True)
    p, group = apc.overwrite_or_skip(p)

    p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+', 
        help="Format of the output files. (default: %(default)s)")

    description = """Parameters related to the parallel computation"""
    p, parallel = apc.parallel_group( p, description=description )

    return p.parse_args(args=argv)

def cselect( f, col, constraint, **kwargs):
    """read file 'f', substitute a columns with noz(z), with z in 'f' itself, and
    save in a file.
    Parameters
    ----------
    f: file object or string
        file containing the catalogue
    col: integer
        column to use for the selection
    constr: string
        selection criterion
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

    col = cat[:, col]
    col = eval( constraint )

    np.savetxt( ofile, cat[col, :], fmt=kwargs['fmt'], delimiter='\t' )
    if(args.verbose == True):
        print( "File '{0}' saved".format(ofile) )

    return None
#end cselect( f, col, constr, **kwargs):

if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])

    #if parallel computation required, check that Ipython.parallel.Client 
    #is in installed and that the ipycluster has been started
    if args.parallel :
        from ipython_parallel import Load_balanced_view as Lbv
        parallel_env = Lbv()  #initialize the object with all my parallen stuff
        args.parallel = parallel_env.is_parallel_enabled()

    #loop through the catalogues and add a columns with n(z)
    if( args.parallel == False ):  #if: parallel
        for fn in args.ifname:  #file name loop
            #substitute n(z)
            cselect( fn, args.column, args.constr, **vars(args) )
    #run the script using the IPython parallel environment 
    else:    #if: parallel
        imports = [ 'import numpy as np', 'import my_functions as mf', ]
        parallel_env.exec_on_engine( imports )

        initstatus = parallel_env.queue_status()  #get the initial status

        #submit the jobs and save the list of jobs
        import os
        runs = [ parallel_env.apply( cselect, os.path.abspath(fn.name), args.column,
            args.constr, **vars(args) ) for fn in args.ifname ]

        if args.verbose :   #if some info is required
            parallel_env.advancement_jobs( runs, update=args.update,
                    init_status=initstatus )
        else:   #if no info at all is wanted
            parallel_env.wait( jobs=runs )  #wait for the end

        #just check for any error
        results = [r.result for r in runs]
        #clear the variable in the parallel environment to avoid filling up memory
        parallel_env.clear_cache()  
    #end if: parallel

    exit()
