#!/usr/bin/python
# -*- coding: utf-8 -*-
#read fits files and save them to ascii tables using 
#fitsio by Erin Scheldon. Case sensitive

import fitsio 

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

    description = """Give a list of fits files, read the data into a numpy
    array and then save them into ascii files """

    p = ap.ArgumentParser(description=description,
            formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument("ifname", nargs='+', type=ap.FileType('r'), 
            action = apc.close_file,
            help="Input file name(s), containing z in one of the columns")

    p = apc.version_verbose( p, '1' )

    p, group = apc.insert_or_replace1(p)
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

    p.add_argument("-c", "--columns", nargs='+', default=["RA", "DEC", "Z"], help="""Read the given
            columns. Can be: integer, string and lists. Operations permitted
            for string entry like 'a+b-1': 'a' and 'b' are column name and '1'
            a number.""")

    p.add_argument("--fmt", action=apc.store_fmt, nargs='+',
            help="Format of the output files.")

    description = """Parameters related to the parallel computation"""
    p, parallel = apc.parallel_group( p, description=description )

    return p.parse_args(args=argv)  
#end def parse(argv)

def get_fileinfo( fnames, allfiles=False, tofile=None ):
    """
    Get a list of file names and print the information stored in the fits
    (HDU and corresponding info)
    """
    if tofile is not None:
        def printto( string ):
            "print string to file"
            tofile.write( str(string)+"\n")
    else:
        def printto( string ):
            "print string to stdout"
            print( string )

    with fitsio.FITS( fnames[0] ) as fits:
        printto( fits[:] )
    if allfiles and len(fnames) > 1:
        for fn in fnames[1:]:
            with fitsio.FITS( fn ) as fits:
                printto( fits[:] )
#end def get_fileinfo( fnames, allfiles=False, tofile=None ):

def fits2ascii( f, columns=None, **kwargs ):
    """
    read a fits file 'f', or some columns if cols is not None and save the
    output into an ascii table.  Operations between columns permitted. The
    output file has the order given in columns

    Signature: 
        fits2ascii( "file.fits" ) #print the whole file
        fits2ascii( "file.fits", columns= [0,4,5,3] )  #read and save the given columns (in that order)
        fits2ascii( "file.fits", columns= ['a', 'b', 'n', 'c'] ) #as before but with column name
        fits2ascii( "file.fits", columns= ['a', 'b+n-1', 'c'] ) #save a 3 columns ascii table with 'a' in the first,
        #the result of the fitsio.read_column('b')+fitsio.read_column('n')-1 in the second and 'c' in the third

    Parameters
    ----------
    f: file object or string
        file containing the catalogue
    columns: None, integer, string, list of integers or strings.
        columns to read. If None the whole file read and saved to an ascii
        table (using np.savetxt) If strings, operation permitted: in this case
        the *name* of the column, *not* the number, must be given

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



if __name__ == "__main__":   #if it's the main

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

    if( args.parallel == False ):  #if: parallel
        pass
    else: #else: parallel
        pass
    #end if: parallel

    exit()

