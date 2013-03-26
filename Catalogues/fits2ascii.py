#!/usr/bin/python
# -*- coding: utf-8 -*-
# read fits files and save them to ascii tables using 
# fitsio by Erin Scheldon. Case sensitive

import fitsio 
import my_functions as mf
import numpy as np
import re

# operations permitted between columns of the fits. If extended add the
# corresponding pattern in compilation of the regular expression
permitted_operations = ['+', '-', '/', '*', '**']

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

    description = """Give a list of fits files, read some of the
    columns into a numpy array and then save them into ascii files. Operations
    between colums permitted but numbers are interpreted as such and not
    columns. Permitted operations: {0}
    """.format(permitted_operations)

    p = ap.ArgumentParser(description=description,
            formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument("ifname", nargs='+', action=apc.file_exists(),
            help="Input file name(s), containing z in one of the columns")

    p = apc.version_verbose( p, '1' )

    p, group = apc.insert_or_replace(p)
    p, group = apc.overwrite_or_skip(p)


    finfo = p.add_argument_group( title="file information", description="""Get
            information from the fits files""")
    finfo.add_argument("--show-columns", action="store_true", help="""Print on
            screen the info for the first file in 'ifname'""") 
    finfo.add_argument("--all-files", action="store_true", help="""Show the info of all
            the files. Ignored if '--show-columns' not used.""")
    finfo.add_argument("--tofile", action="store", type=ap.FileType("w"),
            help="""Write the output to file '%(dest)s' instead than to
            stdout. Ignored if '--show-columns' not used.""")

    p.add_argument("-c", "--columns", nargs='+', default=["RA", "DEC", "Z"],
            help="""Read the given column(s). Can be: integer, string and lists.
            Operations permitted for string entry like 'a+b-1': 'a' and 'b' are column
            name and '1' a number.""")

    p.add_argument("--fmt", action=apc.StoreFmt, nargs='+', default="%7.6e",
            help="Format of the output files.")

    description = """Parameters related to the parallel computation"""
    p, parallel = apc.parallel_group( p, description=description )

    return p.parse_args(args=argv)  
# end def parse(argv)

def get_fileinfo(fnames, allfiles=False, tofile=None):
    """
    Get a list of file names and print the information stored in the fits
    (HDU and corresponding info)
    """
    if tofile is not None:
        def printto(string):
            "print string to file"
            tofile.write(str(string) + "\n")
    else:
        def printto(string):
            "print string to stdout"
            print(string)

    with fitsio.FITS(fnames[0]) as fits:
        printto(fits[:])
    if allfiles and len(fnames) > 1:
        for fn in fnames[1:]:
            with fitsio.FITS(fn) as fits:
                printto(fits[:])
# end def get_fileinfo( fnames, allfiles=False, tofile=None ):

def fits2ascii(fname, selected_columns, operations=None, **kwargs):
    """
    Read some columns a fits file 'fname' and save the output into an ascii table.
    Operations between columns permitted. The output file has the order given
    in columns

    Signature: 
        fits2ascii( "file.fits", 0 )  #read and save the given column
        fits2ascii( "file.fits", [0,4,5,3] )  #read and save the given columns (in that order)
        fits2ascii( "file.fits", ['a', 'b', 'n', 'c'] ) #as before but with column name
        fits2ascii( "file.fits", ['a', 'b+n-1', 'c'] ) #save a 3 columns ascii table with 'a' in the first,
        the result of the fitsio.read_column('b')+fitsio.read_column('n')-1 in the second and 'c' in the third

    Parameters
    ----------
    fname: string
        file containing the catalogue
    selected_columns: integer, string, list of integers or strings.
        columns to read. If strings, operation permitted: in this case
        the *name* of the column, *not* the number, must be given
    operations: regular expression pattern
        pattern containing the list of accepted operators. If None no operation
        done

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
    if kwargs['verbose']:
        print("extracting columns '{0}' from file '{1}'".format(selected_columns, fname)) 
    ofile = mf.create_ofile_name(fname, **kwargs) # create the output file name

    # check if there are operations to execute
    check_operations = []
    if operations is None:
        check_operations = False
    else:
        for c in selected_columns:
            try:
                n_matches = len(operations.findall(c))
            except TypeError:  # The number of the column given
                check_operations.append(False)
            else:  #If a string is given
                if n_matches > 0:
                    check_operations.append(True)
                else:
                    check_operations.append(False)
    # If there are no operations to perform, open the file, read the columns
    # and save them
    if sum(check_operations) == 0:
        cat = fitsio.read(fname, columns=selected_columns)[selected_columns]
    # If there are operations involved
    else:
        with fitsio.FITS(fname) as fits:
            cat = read_with_operations(fits[1], selected_columns, operations)
    np.savetxt(ofile, cat, fmt=kwargs['fmt'], delimiter='\t')
# end def fits2ascii(f, columns, operations=None, **kwargs):

def read_with_operations(fitstable, columns, operations):
    """
    Parameters
    ----------
    fitstable: FITS HDU (from fitsio)
    columns: integer, string, list of integers or strings.
        columns to read. If strings, operation permitted: in this case
        the *name* of the column, *not* the number, must be given
    operations: regular expression pattern
        pattern containing the list of accepted operators. If None no operation
        done

    output
    ------
    catalogue: numpy array
        catalogue with the extracted catalogues from the fits file  
    """
    catalogue = []
    for c in columns:
        if isinstance(c, (int, long)):  # if the number of the column given
            catalogue.append(fitstable.read_column(c))
        else:  # if it is a string
            split_on_operators = [s.strip() for s in operations.split(c)]
            if len(split_on_operators) == 1:  # no operation required
                catalogue.append(fitstable[c][:])
            else:  # do the operations
                # get only the names of the columns
                columns_name = [col for col in split_on_operators 
                        if not is_float(col)]
                # regex for the columns name
                columns_name_pattern =  \
                    re.compile("({0})".format("|".join(columns_name)))
                # substitute the column name with the command to read it
                to_execute = columns_name_pattern.sub(r"fitstable['\1'][:]", c)
                # read the colums and do the operations required
                catalogue.append(eval(to_execute))
    return np.array(catalogue).T # return the catalogue as numpy array
# end def read_with_operations(fitstable, columns, operations)

def is_float(string):
    """
    Check fi the string is a float
    Parameters
    ----------
    string:
        input string to evaluate

    output
    ------
    is_float: bool
    """
    try:
        float(string)
        return True
    except ValueError:
        return False
# end def is_float(string)

if __name__ == "__main__": 

    import sys
    args = parse(sys.argv[1:])

    if args.show_columns:
        get_fileinfo( args.ifname, allfiles = args.all_files, 
                tofile = args.tofile )
        exit()

    #if parallel computation required, check that Ipython.parallel.Client 
    #is in installed and that the ipycluster has been started
    if args.parallel :
        from ipython_parallel import Load_balanced_view as Lbv
        parallel_env = Lbv()  #initialize the object with all my parallen stuff
        args.parallel = parallel_env.is_parallel_enabled()

    # regular expression pattern for the permitted operations
    pattern = re.compile(r"\+|\-|\/|\*{1,2}")

    columns = args.columns

    if( args.parallel == False ):  # if: parallel
        for fn in args.ifname:  # file name loop
            fits2ascii(fn, columns, operations=pattern, **vars(args))
    else: # else: parallel
        # execute some import on all engines
        import os
        # the absolute path and file name of this script
        path, fname = os.path.split(os.path.abspath(sys.argv[0]))
        module = os.path.splitext(fname)[0]
        # command to run on all the engines
        imports = ['import numpy as np', 'import my_functions as mf',
            # add this script directory to the python path
            'import sys', 'sys.path.append("{0}")'.format(path),     
            # import the functions in this file in the namespace
            'from {0} import {1}'.format(module, 'read_with_operations'),
            'from {0} import {1}'.format(module, 'is_float')]
        parallel_env.exec_on_engine(imports)

        initstatus = parallel_env.get_queue_status()  #get the initial status

        #submit the jobs and save the list of jobs
        runs = [parallel_env.apply(fits2ascii, os.path.abspath(fn), columns,
            operations=pattern, **vars(args)) for fn in args.ifname]

        if args.verbose :   #if some info is required
            parallel_env.advancement_jobs(runs, update=args.update,
                    init_status=initstatus)
        else:   #if no info at all is wanted
            parallel_env.wait(jobs=runs)  #wait for the end

        #clear the variable in the parallel environment to avoid filling up memory
        parallel_env.clear_cache()  
    #end if: parallel

    exit()

