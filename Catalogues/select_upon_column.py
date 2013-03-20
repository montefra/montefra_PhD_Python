#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Read one or more files and extract a part of the file according to the
constraints given for a column
"""

import my_functions as mf
import numpy as np
import pandas as pd
import re

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

    description = """  Read one or more files and extract a part of the file according to the constraints given.
        Examples: 
        %(prog)s 'c3 < 4' filename
        %(prog)s '(c3 < 4) | (c3 > 8)' filename
        %(prog)s '(c4 > 8) & (c4 < 4)' filename
        WARNING: The constraint is evaluated with 'eval' and no check is done upon local or global variables."""

    p = ap.ArgumentParser(description=description, formatter_class=ap.RawDescriptionHelpFormatter)

    #p.add_argument("column", action="store", type=int, help="Column to check")
    p.add_argument("constr", action="store", 
        help="""Constraints to be applied to the desired column""")
    p.add_argument("ifname", nargs='+', action=apc.file_exists(),
        help="Input file name(s)")

    p = apc.version_verbose( p, '1.0' )

    p, group = apc.insert_or_replace(p, print_def=True)
    p, group = apc.overwrite_or_skip(p)
    
    p, pandas = apc.pandas_group(p)
 
    p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+', 
        help="Format of the output files. (default: %(default)s)")

    description = """Parameters related to the parallel computation"""
    p, parallel = apc.parallel_group(p, description=description)

    return p.parse_args(args=argv)

def cselect(f, constraint, **kwargs):
    """read file 'f', substitute a columns with noz(z), with z in 'f' itself, and
    save in a file.
    Parameters
    ----------
    f: file object or string
        file containing the catalogue
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
    +pandas: use pandas for the input
    +chunks: chunksize in pandas.read_table
    +fmt: format of the output file
    """
    ofile = mf.create_ofile_name(f, **kwargs) # create the output file name
    if kwargs['verbose']:
        print("Processing file '{}'".format(f))

    if kwargs['pandas']:
        if kwargs['chunks'] is None:
            cat = pd.read_table(f, header=None,
                    skiprows=mf.n_lines_comments(f), sep='\s')
            cat = cat[eval(constraint)]
            np.savetxt(ofile, cat, fmt=kwargs['fmt'], delimiter='\t')
        else:
            chunks = pd.read_table(f, header=None, sep='\s',
                    skiprows=mf.n_lines_comments(f), chunksize=kwargs['chunks'])
            with open(ofile, 'w') as of:
                for cat in chunks:
                    cat = cat[eval(constraint)]
                    np.savetxt(of, cat, fmt=kwargs['fmt'], delimiter='\t')
    else:
        cat = np.loadtxt(f)
        cat = cat[eval(constraint),:]
        np.savetxt(ofile, cat, fmt=kwargs['fmt'], delimiter='\t')
##end cselect( f, col, constr, **kwargs):

if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])

    #if parallel computation required, check that Ipython.parallel.Client 
    #is in installed and that the ipycluster has been started
    if args.parallel :
        from ipython_parallel import Load_balanced_view as Lbv
        parallel_env = Lbv()  #initialize the object with all my parallen stuff
        args.parallel = parallel_env.is_parallel_enabled()

    # create the string to evaluate
    pattern = re.compile(r"c(\d+?)")  #re pattern with the columns name
    if args.pandas:  #expression to evaluate to execute the operation
        to_evaluate = pattern.sub("cat[\\1]", args.constr)
    else:
        to_evaluate = pattern.sub("cat[:,\\1]", args.constr)

    #loop through the catalogues and add a columns with n(z)
    if(args.parallel == False):  #if: parallel
        for fn in args.ifname:  #file name loop
            #substitute n(z)
            cselect(fn, to_evaluate, **vars(args))
    #run the script using the IPython parallel environment 
    else:    #if: parallel
        imports = ['import numpy as np', 'import my_functions as mf', 
                "import pandas as pd"]
        parallel_env.exec_on_engine(imports)

        initstatus = parallel_env.get_queue_status()  #get the initial status

        #submit the jobs and save the list of jobs
        import os
        runs = [parallel_env.apply(cselect, os.path.abspath(fn), to_evaluate,
            **vars(args)) for fn in args.ifname]

        if args.verbose :   #if some info is required
            parallel_env.advancement_jobs(runs, update=args.update,
                    init_status=initstatus)
        else:   #if no info at all is wanted
            parallel_env.wait(jobs=runs)  #wait for the end

        #just check for any error
        results = [r.result for r in runs]
        #clear the variable in the parallel environment to avoid filling up memory
        parallel_env.clear_cache()  
    #end if: parallel

    exit()
