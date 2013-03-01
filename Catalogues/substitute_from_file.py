#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Replace in a (list of) file(s) a column according to the content of 'replacement'
file.
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

    description = """Read two columns c1, c2 from 'replacement' file, then
    replace in 'column' of the input file(s) c1 -> c2."""

    p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument("column", action="store", type=int, help="Column to substitute")
    p.add_argument("replacement", action=apc.file_exists(),
        help="""File containing the substitution to do.""")
    p.add_argument("ifname", nargs='+', action=apc.file_exists(),
        help="Input file name(s)")

    p = apc.version_verbose(p, '1.0')

    p, group = apc.insert_or_replace(p, print_def=True)
    p, group = apc.overwrite_or_skip(p)
    
    p.add_argument("-c", "--repl-columns", nargs='+', action=apc.required_length(1,2), 
            type=int, default=[0,1], help="Columns to read from 'replacement' file")

    p, pandas = apc.pandas_group(p)
 
    p.add_argument("--fmt", default="%7.6e", action=apc.store_fmt, nargs='+', 
        help="Format of the output files. (default: %(default)s)")

    description = """Parameters related to the parallel computation"""
    p, parallel = apc.parallel_group(p, description=description)

    return p.parse_args(args=argv)

def read_replacement(fname, usecols, verbose=False):
    """Read the replacement file and return the two columns
    Parameters
    ----------
    fname: string
        file to read
    usecols: 2 item list
        columns to read in from the input file
    verbose: bool (option)
    
    output
    ------
    replacement: N*2 ndarray
    """
    if verbose:
        print("Reading file '{}'".format(fname))

    class Replacement(object):
        """
        read the file containing the replacement and provide a __call__
        method to get the content of the file
        Parameters
        ----------
        fname: string
            file to read
        usecols: list with one or two elements
            columns to read from the input file
        """
        def __init__(self, fname, usecols):
            self.replacement = np.loadtxt(fname, usecols=usecols)
            self.ncols = len(usecols)

        def _one_column(self, x):
            """
            called by __call__ when only one column is read
            """
            assert x.dtype == int, "'x' must be of integer type"
            return self.replacement[x]
        def _two_columns(self, x):
            """
            called by __call__ when two columns are read
            x casted to float and copied to local variable
            """
            _x = x.astype(np.float)
            for (r0, r1) in self.replacement:
                _x[x==r0] = r1
            return _x
            

        def __call__(self, x):
            """
            Returns the content of the input file at point x.
            If only one column read from the file, x must be an index (a list
            of indeces) of the read array. If the read columns are 2, returns 
            the values in the second one where 'x' matches the first one. 'x'
            is casted to float to avoid problems if passing from indeces to
            float values.
            Parameters
            ----------
            x: int, float or list/ndarray
                value(s) 
            output
            ------
            ndarray with the same shape of x
            """
            if self.ncols == 1:
                return self._one_column(x)
            else:
                return self._two_columns(x)
    # end class Replacement(object):
    
    return Replacement(fname, usecols)
# end def read_replacement(fname, usecols, verbose=False):

def substitute_from_file(fname, col, replacementc, **kwargs):
    """
    Read file 'fname' and replace 'column' matching 'replacement_arr[:,0]' with
    'replacement_arr[:,1]'

    Parameters
    ----------
    fname: file object or string
        file containing the catalogue
    col: integer
        column to use
    replacementc: callable
        function that execute the replacement
        
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
    if kwargs['verbose']:
        print("Processing file '{}'".format(fname))
    ofile = mf.create_ofile_name(fname, **kwargs) # create the output file name

    if kwargs['pandas']:
        if kwargs['chunks'] is None:
            cat = pd.read_table(fname, header=None, skiprows=mf.n_lines_comments(fname), sep='\s')
            cat[col] = replacementc(cat[col])
            np.savetxt(ofile, cat, fmt=kwargs['fmt'], delimiter='\t')
        else:
            chunks = pd.read_table(fname, header=None, sep='\s',
                    skiprows=mf.n_lines_comments(fname), chunksize=kwargs['chunks'])
            with open(ofile, 'w') as fo: #open the output file
                for cat in chunks:  #loop over the chunks
                    cat[col] = replacementc(cat[col])
                    np.savetxt(fo, cat, fmt=kwargs['fmt'], delimiter='\t')

    else:
        cat = np.loadtxt(fname)
        cat[:,col] = replacementc(cat[:,col])
        np.savetxt(ofile, cat, fmt=kwargs['fmt'], delimiter='\t')
#end def substitute_from_file(fname, col, replacementc, **kwargs):

if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])

    if args.parallel:
        raise RuntimeError("The class 'Replacement' can't be sent to the engines.")

    #if parallel computation required, check that Ipython.parallel.Client 
    #is in installed and that the ipycluster has been started
    if args.parallel :
        from ipython_parallel import Load_balanced_view as Lbv
        parallel_env = Lbv()  #initialize the object with all my parallen stuff
        args.parallel = parallel_env.is_parallel_enabled()

    #read in the file with the replacement
    replacement = read_replacement(args.replacement, args.repl_columns,
            verbose=args.verbose)

    #loop through the catalogues and add a columns with n(z)
    if(args.parallel == False):  #if: parallel
        for fn in args.ifname:  #file name loop
            #substitute n(z)
            substitute_from_file(fn, args.column, replacement, **vars(args))
    #run the script using the IPython parallel environment 
    else:    #if: parallel
        imports = ['import numpy as np', 'import my_functions as mf', 
                "import pandas as pd"]
        parallel_env.exec_on_engine(imports)

        initstatus = parallel_env.get_queue_status()  #get the initial status

        #submit the jobs and save the list of jobs
        import os
        runs = [parallel_env.apply(substitute_from_file, os.path.abspath(fn), args.column,
            replacement, **vars(args)) for fn in args.ifname]

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
