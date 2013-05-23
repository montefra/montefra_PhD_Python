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

    p.add_argument("ofname", action="store", 
            help="""Output file name. If more that one element in '--norm_sh'
            is given, the string associated to the various cases is inserted
            after '%(dest)s.'""")

    p.add_argument("ifname", action=apc.file_exists(), nargs='+', 
            help="Input file name(s). The files are ordered, to be able to use 'correct_sh'")

    p = apc.version_verbose(p, '0.1')
    
    p.add_argument("-o", "--overwrite", action="store_true", 
            help="Overwrite the output file name")

    p.add_argument("--before", action="store", help="""Insert the string of
            selected cases and a dot before '%(dest)s' instead of the end of the output
            file name""")

    p.add_argument("-c", "--covariance", action='store', 
            help="""Computes and saves the covariance matrix to file '%(dest)s'.
            If more that one element in '--norm_sh' is given, the string
            associated to the various cases is inserted after '%(dest)s.'""")

    p.add_argument( "--norm-sh", nargs="+", choices=ns_choices,
            default=ns_choices[0], metavar='CHOICES', help="""Cases to work on.
            The part before the '_' regards the power spectrum normalisation,
            the one after the shot noise. The choices are CHOISES={0}. The input
            power spectra are for case '{0[0]}'. If 'all' given, all the cases
            done.""".format(ns_choices))

    p.add_argument("--correct-sh", nargs=3, help="""Corrects the shot noise
            estimated from the randoms (i.e. for the cases '*_ran')
            substituting 'P_SN' with 'n_tot/n_redshift * P_SN': 'n_tot' and
            'n_redshift' are the number of objects that should have redshift
            and the number of objects with measured redshifts. Those two
            numbers are in columns '%(dest)s[1]' and '%(dest)s[2]' of file
            '%(dest)s[0]'. The number of lines in '%(dest)s[0]' must be the
            same as the number of input files and the order assumed to be the
            same as the ordered list of files.""")

    p.add_argument("-t", "--total-number", type=int, help="""If given and
            smaller than the number of input file names, randomly remove
            enought file, and rows in 'correct_sh', in order to get the the
            desired number""")

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

def correct_sn(headers, n_obj, n_red):
    """
    Given a list of pk headers, and two arrays with the total number of objects
    and the number of objects with redshifts, returns the correction to
    subtract from each pk
    """
    for h, no, nz in it.izip(headers, n_obj, n_red):
        yield (no/nz-1) * h.get_SNrandoms()/h.get_N2randoms()

def read_ps(fnlist, choises, sh_correct, total_number):
    """
    Reads and saves into lists the header, the second and third colums for the
    input files, and returns them
    Parameters
    ----------
    fnlist: list of strings
        files  to read
    choises: list of strings
        choises for the normalisation and 
    sh_correct: None or a file name and two integers
        Files and column in the files to use to correct the shot noise
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

    #sort the file list
    fnlist.sort()

    # remove some of the files if required
    if total_number is not None and total_number < len(fnlist):
        to_remove = np.random.randint(0, high=len(fnlist), #random indeces to remove
                size=len(fnlist)-total_number)
        to_remove.sort() #sort it
        for tr in to_remove[::-1]:  #remove the correspoind indeces
            fnlist.pop(tr)       #from the bottom of the list

    #read files and headers
    for fn in fnlist:
        headers.append(PS_header(fn))
        tpk, tpknsh = np.loadtxt(fn, usecols=[1,2]).T
        temppk.append(tpk)
        temppknsh.append(tpknsh)

    #read the file to correct the shot noise if needed
    if sh_correct is not None:
        n_tot, n_redshift = np.loadtxt(sh_correct[0], usecols=[int(i) for i in sh_correct[1:]]).T
        # try to remove lines if some file has been removed
        try:
            n_tot = np.delete(n_tot, to_remove)
            n_redshift = np.delete(n_redshift, to_remove)
        except NameError: #if 'to_remove' is not defined, do nothing
            pass
        if n_tot.size != len(fnlist):
            print("""The number of lines in file {} must be the same as the
                    number of input power spectra""".format(sh_correct[0]))
            exit()
        

    # ns_choices = ["ran_ran", "ran_dat", "dat_ran", "dat_dat", "ran_nosh", "dat_nosh", 'ran_rancorr', 'dat_rancorr']
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
    if 'ran_rancorr' in choises:
        pk['ran_rancorr'] = [p-sn for p,sn in it.izip(temppk, correct_sn(headers, n_tot, n_redshift))] 
    if 'dat_rancorr' in choises:
        pk['dat_rancorr'] = [(p-sn)*rn for p,sn, rn in it.izip(temppk, correct_sn(headers, n_tot, n_redshift), rescale_norm(headers))] 

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
    dictionary
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
    for key, value in pk.iteritems():
        mean[key] = np.mean(value, axis=0)
        stddev[key] = np.std(value, axis=0, ddof=1)
    return mean, stddev

def covariance(pk):
    """
    do the covariance matrix of the power spectra in input
    dictionary
    Parameters
    ----------
    pk: dict of lists of numpy arrays
        power spectra
    output
    ------
    covmat: dict of numpy arrays
        dictionary with the same keys as the input one
    """
    covmat = {}
    for key, value in pk.iteritems():
        covmat[key] = np.cov(value, rowvar=0)
    return covmat

def check_ofiles(ofname, cases, overwrite=False, covfname=None, before=None):
    """
    Create and check the output file names for the power spectra and, if give,
    the covariance matrix
    Parameters
    ----------
    ofname: string
        file name of the output mean power spectrum, if len(cases)==1, base for
        the output file names otherwise.
    cases: list of strings
        cases for which the mean, standard deviation (and covariance) is to be
        computed
    overwrite: bool (optional)
        overwrite the output files if existing
    covfname: string
        file name of the output mean power spectrum, if len(cases)==1, base for
        the output file names otherwise.
    before: string (optional)
        if len(cases)>1, the output file names are created inserting the cases
        before 'before' instead of appenting to the file name

    output
    ------
    ofnames: dict
        dictionary of output file names for the mean power spectra
    covfnames: dict
        only if covfname is not None: dictionary of output file names for the covariances
    """
    from io_custom import file_exists
    class OutputException(Exception):
        def __init__(self, fname):
            """get the file name and create a message with it"""
            string = "File '{}' exists. Delete or rename it, or use the '--overwrite' option"
            self.string = string.format(fname)
        def __str__(self):
            return repr(self.string)
    def _create_fnames(ofname, cases, before):
        if len(cases) == 1: # if only one case is provided
            return {cases[0] : ofname}
        elif before is None:  # if the case is to be appended to the file name
            return {c: ofname+c for c in cases}
        else:
            return {c: ofname.replace(before, c+'.'+before) for c in cases}

    # create the list of output file names
    ofnames = _create_fnames(ofname, cases, before)
    if covfname is not None:
        covfnames = _create_fnames(covfname, cases, before)
    if not overwrite:  # check if all the files exist
        all_files = ofnames.values() if (covfname is None) else ofnames.values()+covfnames.values()  # create a unique list of file names
        for fn in all_files: 
            if file_exists(fn):
                raise OutputException(ofname)
    if covfname is None:
        return ofnames
    else:
        return ofnames, covfnames
# end def check_ofiles(ofname, cases, overwrite=False, covfname=None, before=None):

def print_pk(ofnames, k, mean, std, n_modes, header=None):
    """
    Save the power spectra 'pk' in files 'ofnames'.
    Parameters
    ----------
    ofnames: dict
        output file names
    mean: dict
        mean power spectra
    std: dict
        standard deviations
    k: 1D numpy array
        wavenumbers of the pk
    n_modes: 1D numpy array
        number of modes per k shell
    header: list of 4 1D numpy arrays (optional)
        mean and standard deviation of the data and random sums
    """
    header_pre = "#\t\tsum(w)\tsum(w^2n(z))\tsum(w^2)\n"
    header_post= "#\tk\tmean P(k)\tstddev P(k)\tn_modes\n"
    header_line_start = ["mean_data", "mean_random", "stddev_data", "stddev_random"]
    header_template = "#{0}\t{1[0]:.4f}\t{1[1]:.4f}\t{1[2]:.4f}\n"
    fmt = ['%7.6e', '%7.6e', '%7.6e', '%d' ]
    for key, fvalue in ofnames.iteritems(): #iter through the file names
        with open(fvalue, 'w') as f:
            if header is not None: #write the header if not given
                f.write(header_pre)
                for hls, h in it.izip(header_line_start, header):
                    f.write(header_template.format(hls, h))
            f.write(header_post)
            np.savetxt(f, np.vstack([k, mean[key], std[key], n_modes]).T, delimiter='\t', fmt=fmt)
# end def print_pk(ofnames, k, mean, std, n_modes, header=None):

def print_cov(ofnames, cov, k=None):
    """
    Save to files the covariance matrices.
    Parameters
    ----------
    ofnames: dict
        output file names
    cov: dict
        covariance matrices
    k: 1D numpy array (optional)
        values of k for the covariance matrix. Written as header if provided
    """
    for key, fvalue in ofnames.iteritems():
        with open(fvalue, 'w') as f:
            f.write("#k: "+"\t".join(['{0:7.6e}'.format(kk) for kk in k])+"\n")
            np.savetxt(f, cov[key], delimiter='\t', fmt='%7.6e')

if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])

    # check the choices 
    if isinstance(args.norm_sh, str):
        args.norm_sh = [args.norm_sh,]
    if 'all' in args.norm_sh:
        args.norm_sh = ns_choices[:-1]
    else:  #get unique elements
        args.norm_sh = list(set(args.norm_sh))
    # if correct_sh given, add ran_rancorr and/or dat_rancorr if ran_ran and/or dat_ran are present
    if args.correct_sh is not None:
        if 'ran_ran' in args.norm_sh:
            args.norm_sh.append('ran_rancorr')
        if 'dat_ran' in args.norm_sh:
            args.norm_sh.append('dat_rancorr')

    # create and check the output files
    ofnames = check_ofiles(args.ofname, args.norm_sh, overwrite=args.overwrite,
            covfname=args.covariance, before=args.before)
    if args.covariance is not None:
        ofnames, covfnames = ofnames

    # read the files
    if args.verbose:
        print("Reading the power spectra")
    pk, headers = read_ps(args.ifname, args.norm_sh, args.correct_sh, args.total_number)
    
    # read k and number of modes
    k, n_modes = np.loadtxt(args.ifname[0], usecols=[0,-1]).T

    if args.verbose:
        print("Compute the mean and standard deviations")
    # do the mean and standard deviation of the headers
    header_mean_std = mean_std_headers(headers)
    # do the mean and standard deviation of the power spectra
    meanpk, stdpk = mean_std_pks(pk)
    if args.verbose:
        print("Print the mean and standard deviations")
    print_pk(ofnames, k, meanpk, stdpk, n_modes, header=header_mean_std)

    # do the covariance
    if args.covariance is not None:
        if args.verbose:
            print("Compute the covariance matrices")
        covmat = covariance(pk)
        if args.verbose:
            print("Print the covariance matrices")
        print_cov(covfnames, covmat, k=k)


    exit()
        
