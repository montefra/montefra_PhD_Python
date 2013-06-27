#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Read a list of files containing power spectra and adjust the normalisation
and/or shot noise 
"""
import argparse as ap
import io_custom as ioc
import my_functions as mf
import numpy as np
from pk_utils import PS_header
import ps_mean_err as ps_me

ps_me.ns_choices.remove("ran_nosh")
ps_me.ns_choices.remove("dat_nosh")

class read_n(ap.Action):
    """If the file name and the elements to read are given, read then right ahead"""

    def __call__(self, parser, namespace, values, option_string=None):
        try:
            usecols = [int(i) for i in values[1:]]  #convert the second and third element in integers
        except ValueError:
            msg="The second and third element of '{}' must be integers".format(option_string)
            raise ap.ArgumentTypeError(msg)
        n_tr = np.loadtxt(values[0], usecols=usecols)  #number tot and redshift
        if n_tr.ndim == 2: #if more than one line read
            n_tr = n_tr[0]
        # as an ascii file is read, no sense to check for higher dimensionalities
        setattr(namespace, self.dest, n_tr) #set the two numbers in the namespace

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
    import argparse_custom as apc
    
    description = """    Change the mean and/or the amplitude of the input
    files
    from the same FFT grid and in the wavenumbers. The structure of the input 
    files must be like the following:
        # 	 sum(w) 	 sum(w^2n(z)) 	 sum(w^2)
        #data	 ###             ###             ###   
        #random	 ###             ###             ###   
        #	k	P(k)    	P(k)+noise      n_modes
        [...]
    The output files will have the same structure
    """
    p = ap.ArgumentParser(description=description, formatter_class=apc.RawDescrArgDefHelpFormatter)

    p.add_argument("ifname", action=apc.file_exists(), nargs='+', 
            help="""Input file name(s).""")

    p = apc.version_verbose(p, '1')

    p.add_argument("--before", action="store", help="""Insert the string of
            selected cases and a dot before '%(dest)s' instead of the end of the output
            file name""")

    p, group = apc.overwrite_or_skip(p)
    
    p.add_argument( "--norm-sh", nargs="+", choices=ps_me.ns_choices,
            default=ps_me.ns_choices[0], metavar='CHOICES', help="""Cases to work on.
            The part before the '_' regards the power spectrum normalisation,
            the one after the shot noise. The choices are CHOISES={0}. The input
            power spectra are for case '{0[0]}'. If 'all' given, all the cases
            done.""".format(ps_me.ns_choices))

    p.add_argument("--correct-sh", nargs=3, action=read_n,
            help="""'%(dest)s[0]' contains one line (header ignored) whose
            '%(dest)s[1]' and '%(dest)s[2]' elements are 'n_tot' and 'n_redshift',
            used to correct the shot noise estimated from the randoms substituting
            'P_SN' with 'n_tot/n_redshift * P_SN'. If given, the file is read
            directly and the two numbers are returned. Any line besides the first
            id ignored""")

    p.add_argument("--fmt", default="%7.6e", action=apc.StoreFmt, nargs='+',
            help="Format of the output files")

    return p.parse_args(args=argv)
# end def parse(argv)

def read_file(fname):
    """
    read the file 'fname' saving the full header, the power spectrum table and the content of the header.
    Parameters
    ----------
    fname: string
        file name
    output
    ------
    full_header: list
        full header from the file
    pk: ndarray
        power spectrum table
    header: PS_header object
        actual content of the header
    """
    #read the file name and the header
    with open(fname, 'r') as f:
        header_lines = mf.n_lines_comments(f) #get the size of the header
        full_header = [f.readline() for i in xrange(header_lines)] #save the header as it is
        pk = np.loadtxt(f)  #read the power spectrum
        header = PS_header(f) # get the interesting parts of from the header

    return full_header, pk, header

def read_resc_print(fname, choises, sh_correct=None, before=None,
        overwrite=False, skip=False, fmt='%7.6e'):
    """
    Read the file and its header, do the required rescaling and save the result
    Parameters
    ----------
    fname: string
        file name
    choises: list of strings
        choises for the normalisation and 
    sh_correct: 2 floats
        n_tot and n_redshift used to correct the shot noise   
    before: string
        when creating the output file name, insert the choise and a dot before .before
    overwrite, skip: bool
        overwrite or skip the current file if existing
    fmt: string or list of strings
        format for np.savetxt
    """
    
    full_header, pk, header = read_file(fname)

    SN_resc = list(ps_me.rescale_sn([header,]))  #to add to change shot noise
    N_resc = list(ps_me.rescale_norm([header,])) #to multiply to rescale amplitude
    if sh_correct is not None:
        n_tot, n_redshift = sh_correct
        #to substract to correct the shot noise from the randoms
        SN_corr = list(ps_me.correct_sn([header,], [n_tot,],[n_redshift,]))

    for c in choises:
        #create the output file name
        if before is None:
            ofile = fname+'.'+c
        else:
            ofile = fname.replace(before, c+'.'+before)
        #check if output file exists
        if ioc.file_exists(ofile):
            if overwrite is False:
                print("The output file '{}' already exists, aborting".format(ofile))
                exit()
            elif skip is True:
                continue  #go to the next iteration of the choises loop

        #create the power spectrum 
        opk = np.copy(pk) #copy the array into the output
        if c == 'ran_ran':
            pass  #do nothing
        if c == 'ran_dat':
            opk[:,1] += SN_resc
        if c == 'dat_ran':
            opk[:,1] *= N_resc
        if c == 'dat_dat':
            opk[:,1] = (opk[:,1]+SN_resc)*N_resc
        if c == 'ran_rancorr':
            opk[:,1] -= SN_corr
        if c == 'dat_rancorr':
            opk[:,1] = (opk[:,1]-SN_corr)*N_resc
        
        with open(ofile, 'w') as of:
            of.writelines(full_header)
            np.savetxt(of, opk, delimiter='\t', fmt=fmt)


if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])

    #check the choises
    args.norm_sh = ps_me.check_choises(args.norm_sh, args.correct_sh)
    if args.verbose:
        print("Choises checked")

    for fn in args.ifname:  # loop over the file names
        if args.verbose:
            print("process file '{}'".format(fn))
        read_resc_print(fn, args.norm_sh, sh_correct=args.correct_sh,
                before=args.before, overwrite=args.overwrite, skip=args.skip,
                fmt=args.fmt)

    exit()

