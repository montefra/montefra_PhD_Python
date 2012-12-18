#!/usr/bin/python
# -*- coding: utf-8 -*-
#read fits files and save them to ascii tables using 
#pyfits

import pyfits as pf

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

    p.add_argument("ifname", action="store", nargs='+', type=ap.FileType('r'),
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

