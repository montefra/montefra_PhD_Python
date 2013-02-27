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
    p.add_argument("ifname", action=apc.file_exists(), nargs='+', 
            help="""Input file name(s)""")

    p = apc.version_verbose(p, '0.1')

    p.add_argument("-s", "--swap", action="store_true", 
            help="Swap 'from' with 'to'")
    p.add_argument("--substitute", action="store", type=float,
            help="""If given, substitutes the content of the columns
            'from_cols' not contained in 'to_cols' with '%(dest)s'. This is
            executed only if 'swap' is *False* and after moving""")

    p, group = apc.insert_or_replace(p, print_def=True)
    p, group = apc.overwrite_or_skip(p)

    p, pandas = apc.pandas_group(p)
    
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
    +swap: swap the columns instead of just moving
    +pandas: use pandas for the input
    +chunks: chunksize in pandas.read_table
    +substitute: value to substitute in 'from_columns' that are not in 'to_columns'
    +fmt: format of the output file
    """
    ofile = mf.create_ofile_name(f, **kwargs) # create the output file name
    if kwargs['verbose']:
        print("Processing file '{}'".format(f))
    #get the columns to move or swap
    substitute = None
    temp_from_cols, temp_to_cols = from_columns[:], to_columns[:]
    if kwargs['swap']: # swap the colums
        temp_from_cols.extend(to_columns)
        temp_to_cols.extend(from_columns)
    else:
        if kwargs['substitute'] is not None: #find the columns that need substitution
            substitute = [i for i in temp_from_cols if i not in temp_to_cols]
            if len(substitute) == 0:
                substitute = None #set back to None if there are no columns with value to substitute

    if kwargs['pandas']:
        use_pandas(f, ofile, temp_from_cols, temp_to_cols, substitute, **kwargs)
    else:
        cat = np.loadtxt(f)
        cat[:,temp_to_cols] = cat[:,temp_from_cols]
        if substitute is not None:
            cat[:, substitute] = kwargs['substitute']
        np.savetxt(ofile, cat, fmt=kwargs['fmt'], delimiter='\t')
#end def move_columns(f, from_cols, to_cols, **kwargs):

def use_pandas(fin, fout, from_columns, to_columns, sub, **kwargs):
    """
    Reads the file using pandas read_table and save it using savetxt. The rest is like 'convert_save'
    Parameters
    ----------
    f: file object or string
        file containing ra, dec and z
    fout: string
        output file name
    from_columns: list of ints
        list of columns to copy
    to_columns: list of ints
        list of columns where to copy
    sub: list
        where to substitute in the catalogue
    kwargs: keyword arguments
    output
    ------
    none

    used kwargs that affects the function
    +pandas: use pandas for the input
    +chucks: read, elaborate and print file in chunks of size 'chunks'
    +substitute: value to substitute in 'from_columns' that are not in 'to_columns'
    +fmt: format of the output file
    """

    if kwargs['chunks'] is None:  #read the whole file in one go
        cat = pd.read_table(fin, header=None, sep='\s', skiprows=mf.n_lines_comments(fin))
        cat[to_columns] = cat[from_columns]
        if sub is not None:
            cat[sub] = kwargs['substitute']
        np.savetxt(fout, cat, fmt=kwargs['fmt'], delimiter='\t')

    else:  #read the file in chuncks
        chunks = pd.read_table(fin, header=None, sep='\s',
                skiprows=mf.n_lines_comments(fin), chunksize=kwargs['chunks'])
        with open(fout, 'w') as fo: #open the output file
            for cat in chunks:  #loop over the chunks
                cat[to_columns] = cat[from_columns]
                if sub is not None:
                    cat[sub] = kwargs['substitute']
                np.savetxt(fo, cat, fmt=kwargs['fmt'], delimiter='\t')
#end def use_pandas(fin, fout, distance, **kwargs)

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
        for fn in args.ifname:  #file name loop
            move_columns(fn, args.from_cols, args.to_cols, **vars(args))
    #run the script using the IPython parallel environment 
    else:    #if: parallel
        import os # the absolute path and file name of this script
        path, fname = os.path.split(os.path.abspath(sys.argv[0]))
        imports = [ 'import numpy as np', 'import my_functions as mf', 
                'import pandas as pd',
                #add the script directory to the python path
                'import sys', 'sys.path.append("{0}")'.format(path),     
                #import the desired function in the namespace
                'from {0} import *'.format(os.path.splitext(fname)[0])]  
        parallel_env.exec_on_engine(imports)

        initstatus = parallel_env.get_queue_status()  #get the initial status

        #submit the jobs and save the list of jobs
        import os
        runs = [parallel_env.apply(move_columns, os.path.abspath(fn),
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
