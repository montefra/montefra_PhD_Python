#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Read a list of files containing power spectra and compute mean, standard deviation
and, if required, the covariance
"""
import itertools as it
import numpy as np

# norm_sh argument choices
ns_choices = ["ran_ran", "ran_dat", "dat_ran", "dat_dat", "ran_nosh", "dat_nosh", "all"]
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
    
    description = """    Computes the mean, standard deviation and, if reqired, the covariance 
    from a list of power spectra. The power spectra are assumed to be computed 
    from the same FFT grid and in the wavenumbers. The structure of the input 
    files must be like the following:
        # 	 sum(w) 	 sum(w^2n(z)) 	 sum(w^2)
        #data	 ###             ###             ###   
        #random	 ###             ###             ###   
        #	k	P(k)    	P(k)+noise      n_modes
        [...]
    The output files with the mean and standard deviation will have the following structure
        #       	 sum(w) 	 sum(w^2n(z)) 	 sum(w^2)
        #mean_data	 ###             ###             ###  
        #mean_random	 ###             ###             ###  
        #stddev_data	 ###             ###             ###  
        #stddev_random	 ###             ###             ###  
        #	k	mean P(k)    	stddev P(k)      n_modes
        [...]
    """
    p = ap.ArgumentParser(description=description, formatter_class=apc.RawDescrArgDefHelpFormatter)

    p.add_argument("ofname", action="store", help="Output file name")

    p.add_argument("ifname", action=apc.file_exists(), nargs='+', 
            help="Input file name(s)")

    p = apc.version_verbose(p, '0.1')
    
    p.add_argument("-o", "--overwrite", action="store_true", 
            help="Overwrite the output file name")

    p.add_argument("-c", "--covariance", action='store', 
            help="Computes and saves the covariance matrix to file '%(dest)s'") 

    p.add_argument( "--norm_sh", nargs="+", action=apc.required_length(1, 4), choices=ns_choices,
            default=ns_choices[0], metavar='CHOICES', help="""Cases to work on.
            The part before the '_' regards the power spectrum normalisation,
            the one after the shot noise. The choices are CHOISES={0}. The input
            power spectra are for case '{0[0]}'. If 'all' given, all the cases
            done.""".format(ns_choices))

    return p.parse_args(args=argv)
# end def parse(argv)

def rescale_norm(headers):
    """
    Give a list of pk headers, returns an iterator of ratios of the randoms and data
    normalisation to multiply to each pk
    """
    for h in headers:
        yield h.get_N2randoms()/h.get_N2data() 

def rescale_sn(headers):
    """
    Give a list of pk headers, returns the an iterator of differences of the shot
    noise from the randoms and the data, normalised with the randoms, to add to each pk
    """
    for h in headers:
        yield (h.get_SNrandoms() - h.get_SNdata())/h.get_N2randoms()

def read_ps(fnlist, choises):
    """
    Reads and saves into lists the header, the second and third colums for the
    input files, and returns them
    Parameters
    ----------
    fnlist: list of strings
        files  to read
    choises: list of strings
        choises for the normalisation and 
    output
    ------
    pk: dict of list of 1D numpy arrays
        power spectra
    headers: list of PS_header instances
    """
    from pk_utils import PS_header
    pk= {}  # power spectra
    headers = []  #PS_header instances
    temppk, temppknsh = [],[]

    #read files and headers
    for fn in fnlist:
        headers.append(PS_header(fn))
        tpk, tpknsh = np.loadtxt(fn, usecols=[1,2]).T
        temppk.append(tpk)
        temppknsh.append(tpk)

    # ns_choices = ["ran_ran", "ran_dat", "dat_ran", "dat_dat", "ran_nosh", "dat_nosh"]
    # normalisation: randoms; shot noise: random 
    if ns_choices[0] in choises:
        pk[ns_choices[0]] = temppk[:]
    # normalisation: randoms; shot noise: data
    if ns_choices[1] in choises:
        pk[ns_choices[1]] = [p+sn for p,sn in it.izip(temppk, rescale_sn(headers))]
    # normalisation: data; shot noise: random
    if ns_choices[2] in choises:
        pk[ns_choices[2]] = [p*rn for p,rn in it.izip(temppk, rescale_norm(headers))]
    # normalisation: data; shot noise: data
    if ns_choices[3] in choises:
        pk[ns_choices[3]] = [(p+sn)*rn for p,rn, sn in it.izip(temppk, rescale_norm(headers), rescale_sn(headers))]
    # normalisation: randoms; shot noise: none
    if ns_choices[4] in choises:
        pk[ns_choices[4]] = temppknsh[:]
    # normalisation: data; shot noise: none
    if ns_choices[5] in choises:
        pk[ns_choices[5]] = [p*rn for p,rn in it.izip(temppknsh, rescale_norm(headers))]

    return pk, headers
# end def read_ps(fnlist, choises):

def mean_std_headers(headers):
    """
    computes the mean and std deviation of the sums contained in the headers of the power spectra
    Parameters
    ----------
    headers: list of PS_header instances
    output
    ------
    mean_std: list
        mean of data and randoms, standard deviation of data and randoms. 
        Each element is a size 3 numpy array
    """
    mean_data = np.mean([h.data for h in headers], axis=0)
    std_data = np.std([h.data for h in headers], axis=0)
    mean_random = np.mean([h.random for h in headers], axis=0)
    std_random = np.std([h.random for h in headers], axis=0)

    return [mean_data, mean_random, std_data, std_random]

def mean_std_pks(pk):
    """
    do the mean and standard deviation of the power spectra in the input
    dictionaries
    Parameters
    ----------
    pk: dict of lists of numpy arrays
        power spectra
    output
    ------
    mean: dict of numpy arrays
    stddev: dict of numpy arrays
    dictionaries with the same keys as the input one
    """
    mean, stddev = {}, {}

    return mean, stddev

def covariance(pk):
    pass

if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])


    # check the choices 
    if 'all' in args.norm_sh:
        args.norm_sh = ns_choices[:-1]
    else:  #get unique elements
        args.norm_sh = list(set(args.norm_sh))

    # check if output file exists
    if not args.overwrite:
        from io_custom import file_exists
        if file_exists(args.ofname):
            print("File '{}' exists. Delete or rename it, or use the '--overwrite' option".format(args.ofname))

    #read the files
    if args.verbose:
        print("Reading the power spectra")
    pk, headers = read_ps(args.ifname)
    
    #read k and number of modes
    k, n_modes = np.loadtxt(args.ifname[0], usecols=[0,-1]).T

    #do the mean and standard deviation of the headers
    header_mean_std = mean_std_headers(headers)

    #do the mean and standard deviation of the power spectra

    exit()
        
