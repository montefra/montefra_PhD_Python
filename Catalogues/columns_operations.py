#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
execute operations between columns and save in a given column
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

    description = """Execute operations between columns and save in a column.
    Example: c3+c6-1: add column 3 to column 6 and subtract 1."""

    p = ap.ArgumentParser(description=description,)

    p.add_argument("operation", action="store", 
            help="""Operation to execute between columns.""")
    p.add_argument("to_col", action="store", type=int,
            help="""Column where the operation is saved. If it is larger than
            the number of columns, the result is appended""")
    p.add_argument("ifname", nargs='+', action=apc.file_exists(),
            help="""Input file name(s)""")

    p = apc.version_verbose(p, '0.1')

    p.add_argument("-s", "--substitute", action="store", type=float,
            help="""If given, substitutes the content of the columns used for
            the operation with '%(dest)s'. This is executed before copying the
            result of the operation to the desired column""")

    p, group = apc.insert_or_replace1(p, print_def=True)
    p, group = apc.overwrite_or_skip(p)

    p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+', 
        help="Format of the output files.")

    description = """Parameters related to the parallel computation"""
    p, parallel = apc.parallel_group( p, description=description )

    return p.parse_args(args=argv)

def columns_operations(f, operations, to_column, **kwargs):
    """
    read file, performe the desired operations between columns, save it and
    write to a file.

    Signature: 
        columns_operations("file.dat", 'c3+c6-1', 4) 
        #read file, add content of column 3 and column 6, subtract 1 and save
        the result in column 4. 

    Parameters
    ----------
    f: file object or string
        file containing the catalogue
    operations: string
        columns with operations to perform
    to_columns int
        save there the result of operations
    output
    ------
    none

    accepted kwargs that affects the function
    +substitute: substitute the content of the columns involved with operations
        with this value, if not None
    +verbose: verbose mode [True|False] 
    +replace: replace string *replace[0]* with *replace[1]* in f.name
    +insert: insert string *insert[0]* before *insert[1]* in f.name
    +skip: existing file names skipped [True|False]
    +overwrite: existing file names overwritten [True|False]
    +fmt: format of the output file
    """
    ofile = mf.create_ofile_name(f, **kwargs) # create the output file name

    # read the input catalogue
    cat = pd.read_table(f, header=None, skiprows=mf.n_lines_comments(f)) 

    pattern = re.compile(r"c(\d+?)")  #re pattern with the columns name
    new_column = eval(pattern.sub("cat[\\1]", operations)) # execute the operation

    # substitute the columns used for the operations with the given value
    if kwargs['substitute'] is not None:
        columns_read = [int(s) for s in pattern.findall(operations)]
        cat[columns_read] = kwargs['substitute']
    cat[to_column] = new_column #copy the result of the operation

    # save the converted catalogue
    np.savetxt(ofile, cat, fmt=kwargs['fmt'], delimiter='\t')

if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])
    
    #if parallel computation required, check that Ipython.parallel.Client 
    #is in installed and that the ipycluster has been started
    if args.parallel :
        from ipython_parallel import Load_balanced_view as Lbv
        parallel_env = Lbv()  #initialize the object with all my parallen stuff
        args.parallel = parallel_env.is_parallel_enabled()

    if(args.parallel == False):  # if: parallel
        for fn in args.ifname:  #file name loop
            columns_operations(fn, args.operation, args.to_col, **vars(args))
    else: # else: parallel
        imports = ['import numpy as np', 'import my_functions as mf', 
                "import pandas as pd"]
        parallel_env.exec_on_engine(imports)

        initstatus = parallel_env.get_queue_status()  #get the initial status

        #submit the jobs and save the list of jobs
        import os
        runs = [parallel_env.apply(columns_operations, os.path.abspath(fn), args.operation,
            args.to_col, **vars(args)) for fn in args.ifname]

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
