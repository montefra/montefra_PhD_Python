#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Move or swap some column in a (list of) file(s)
"""

import my_functions as mf
import numpy as np
import pandas as pd

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

    p.add_argument("from_cols", action="store", type=apc.int_or_list, 
            help="""Column(s) to copy. (if more than one, give them as a string
            of integers separated by spaces)""")
    p.add_argument("to_cols", action="store", type=apc.int_or_list,
            help="""Where to copy the column(s). (if more than one, give them
            as a string of integers separated by spaces)""")
    p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'),
            help="""Input file name(s)""")

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

def move_columns(f, from_columns, to_columns, **kwargs):
    """
    Read file 'f', substitute the content of columns 'to_cols' with the one of
    'from_cols'
    Parameters
    ----------
    f: file object or string
        file containing the catalogue
    from_columns: list of ints
        list of columns to copy
    to_columns: list of ints
        list of columns where to copy
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
    ofile = mf.create_ofile_name(f, **kwargs) # create the output file name

    cat = pd.read_table(f, header=None, skiprows=mf.n_lines_comments(f), sep='\s') 

    cat[to_columns] = cat[from_columns]

    np.savetxt(ofile, cat, fmt=kwargs['fmt'], delimiter='\t')
    if(kwargs['verbose'] == True):
        print("File '{0}' saved".format(ofile))
#end def move_columns(f, from_cols, to_cols, **kwargs):

def swap_columns(f, from_columns, to_columns, **kwargs):
    """
    Read file 'f', swap the content of columns 'from_cols' with the one of
    'to_cols'
    Parameters
    ----------
    f: file object or string
        file containing the catalogue
    from_columns: list of ints
        first set of columns
    to_columns: list of ints
        second set of columns
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
    ofile = mf.create_ofile_name(f, **kwargs) # create the output file name

    cat = pd.read_table(f, header=None, skiprows=mf.n_lines_comments(f), sep='\s') 

    # swap the colums
    temp_from_cols, temp_to_cols = from_columns[:], to_columns[:]
    temp_from_cols.extend(to_columns)
    temp_to_cols.extend(from_columns)
    cat[temp_to_cols] = cat[temp_from_cols]

    np.savetxt(ofile, cat, fmt=kwargs['fmt'], delimiter='\t')
    if kwargs['verbose'] == True:
        print("File '{0}' saved".format(ofile))
#end def swap_columns(f, from_cols, to_cols, **kwargs):

if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])

    if len(args.from_cols) != len(args.to_cols):
        print("""It is not possible to copy (swap) {0} columns into (with) {1}
                columns""".format(len(args.from_cols), len(args.to_cols)))

    #if parallel computation required, check that Ipython.parallel.Client 
    #is in installed and that the ipycluster has been started
    if args.parallel:
        from ipython_parallel import Load_balanced_view as Lbv
        parallel_env = Lbv()  #initialize the object with all my parallen stuff
        args.parallel = parallel_env.is_parallel_enabled()

    #loop through the catalogues and add a columns with n(z)
    if(args.parallel == False):  #if: parallel
        if args.swap:
            for fn in args.ifname:  #file name loop
                swap_columns(fn, args.from_cols, args.to_cols, **vars(args))
        else:
            for fn in args.ifname:  #file name loop
                move_columns(fn, args.from_cols, args.to_cols, **vars(args))
    #run the script using the IPython parallel environment 
    else:    #if: parallel
        imports = [ 'import numpy as np', 'import my_functions as mf', 
                'import pandas as pd']
        parallel_env.exec_on_engine(imports)

        initstatus = parallel_env.get_queue_status()  #get the initial status

        #submit the jobs and save the list of jobs
        import os
        if args.swap:
            runs = [parallel_env.apply( swap_columns, os.path.abspath(fn),
                args.from_cols, args.to_cols, **vars(args)) for fn in args.ifname]
        else:
            runs = [parallel_env.apply( move_columns, os.path.abspath(fn),
                args.from_cols, args.to_cols, **vars(args)) for fn in args.ifname]

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
