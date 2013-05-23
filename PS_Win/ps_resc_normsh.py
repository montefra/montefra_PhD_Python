#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Read a list of files containing power spectra and adjust the normalisation
and/or shot noise 
"""
import itertools as it
import numpy as np
import ps_mean_err as ps_me

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

    p = apc.insert_or_replace(p)

    p = apc.overwrite_or_skip(p)
    
    p.add_argument( "--norm_sh", nargs="+", choices=ps_me.ns_choices,
            default=ps_me.ns_choices[0], metavar='CHOICES', help="""Cases to work on.
            The part before the '_' regards the power spectrum normalisation,
            the one after the shot noise. The choices are CHOISES={0}. The input
            power spectra are for case '{0[0]}'. If 'all' given, all the cases
            done.""".format(ps_me.ns_choices))

    p.add_argument("--correct_sh", nargs=3, help="""Corrects the shot noise
            estimated from the randoms (i.e. for the cases '*_ran')
            substituting 'P_SN' with 'n_tot/n_redshift * P_SN': 'n_tot' and
            'n_redshift' are the number of objects that should have redshift
            and the number of objects with measured redshifts. Those two
            numbers are in columns '%(dest)s[1]' and '%(dest)s[2]' of file
            '%(dest)s[0]'. The number of lines in '%(dest)s[0]' must be the
            same as the number of input files and the order assumed to be the
            same as the ordered list of files.""")

    return p.parse_args(args=argv)
# end def parse(argv)


if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])

    #check the choises
    args.norm_sh = check_choises(args.norm_sh, args.correct_sh)


    exit()

