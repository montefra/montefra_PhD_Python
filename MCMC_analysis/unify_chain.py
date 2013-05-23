#!/usr/bin/python
# -*- coding: utf-8 -*-

import glob   #allows search of files with unix wildcards
import itertools as it
import my_statistic as ms
import numpy as np    # import numpy
import optparse as op  #import optsparse: allows nice command line option handling
import os    #contain OS dependent stuffs: it helps with sistem portability  
import re   #regular expressions
import sys    #mondule sys


def options(p):
    """
    This function accept the option parser istance and fill it with options.

    The version must be set when creating the optparse istance

    Parameters
    ----------
    p: option parser istance

    output: tuple with the options and the arguments
    ---------
    """

    p.set_usage("""%prog [options] file_root
    This program reads the files of the COSMOMC containing the chain results and the parameter name definition.
    It discard the first part of each chain and merge the remaining parts of the chains in a unique file.
    If required returns a file with the 1D marginalized constraint for all the parameters.
    The output files has the same root name as the input ones.
    The file names is based upon the COSMOMC style (http://cosmologist.info/cosmomc).
    The common part of the file names 'file_root' must be given as first argument. 'file_root' must contain the full absolute or relative file path.
    The chain have a structure like: n_times likelihood parameters
    """)   # basic informations of the program

    p.add_option("-v", "--verbose", action="store_true", dest="verbose",
            help="Produces more output.")  # verbose option
    p.add_option("-f", "--fraction", action="store", dest="fraction",
            type="float", default="0.5", help="""The first 'FRACTION' of the
            chain that is going to be discarted. 'FRACTION' must be [0,1].
            [default: %default]""")
    p.add_option("-m", "--marginalized", action="store_true", dest="marg",
            help="""If given, the program reads the file containing the
            parameter names and writes a file with the 1D marginalized
            constraints.""")
    p.add_option("-p", "--percentile", action="append", type="float",
            dest="percentile", help="""It allow to set the percentile for the
            output 1D marginalized contraints. If called more then once allows
            to give more than one percentile value. 'PERCENTILE' must be
            [0,100]. If PERCENTILE=68, gives back the lower and upper 1-sigma
            constraints.  [default: 68.0]. To provide more the one percentile,
            call the option more than once""")
    p.add_option("--extension", action="store", dest="ext", type="string",
            default="txt", help="""Extension of the input and output files.
            [default: %default]""")
    p.add_option("--fmt", action="store", dest="fmt", type="string",
            default="%8.7e", help="""Formatting of the output chain [default:
            %default]""")
    p.add_option("-d", "--discard", action="store", dest="discard",
            type="string", help="""Discard files containing a certain string""")
    p.add_option("-o", "--output-dir", action="store", dest="outdir",
            type="string", help="""If an output directory given, the input path
            is substitutes with the output one.""")
    p.add_option("-c", "--convergence", action="store_true", dest="conv",
            default=False, help="""If selected, the Gelman & Rubin convergency
            test run. The value of R for every parameters will be printed in
            the file 'file_root.conv'""")

    return p.parse_args()

def check_optarg(options, args):
    """
    Check if the input options and parameters are ok.
    If the percentile scores are not given it sets the default.
    
    Parameters
    ----------
    options: list
        list of options returned by optparse
    args: list
        list of arguments returned by optparse

    output: list
        list of options
    """

    if(len(args) < 1):   #file_root required
        print """One argument is required by the program
        Type './%s -h' for more informations""" % os.path.basename(sys.argv[0])
        sys.exit(10)
    if(options.fraction < 0. or options.fraction > 1.):
        print "The fraction of the chain to be discarted must be a number in the range [0,1]"
        sys.exit(11)
    if(options.marg == True):
        if(options.percentile == None):
            options.percentile = [68.]
            if(options.verbose == True):
		print "percentile has been set to its default value %f" % options.percentile[0] 
        for i in options.percentile:
            if(i<=0 or i>=100):
		print "Percentile must be (0,100)"
		sys.exit(12)
        options.percentile.sort() # sort the percentile from the input
    else:
        if(options.percentile != None):
            print """The option -p or --percentile requires that -m or --marginalized is 'True'
            Type './%s -h' for more informations""" % os.path.basename(sys.argv[0])
            sys.exit(13)

    return options

def get_file_names(options, args):
    """
    Get the file names and remove the output files if existing and the ones to be discarded
    
    Parameters
    ----------
    options: list
        list of options returned by optparse
    args: list
        list of arguments returned by optparse

    output: list
        list of file names
    """

    flist = glob.glob(args[0]+"*."+options.ext) #list all the files that contain chains

    tbr=[]   #list of to be removed files
    exp = re.compile(os.path.basename(args[0])+".*_total\."+options.ext)  #template of the output file name for the full chain
    for f in flist:   #search for output files in the file list
        if(exp.search(f) != None):
            tbr.append(f)
    if(options.discard != None):   #discard files that contains a given string
        for f in flist:
            if(f.find(options.discard) != -1):
                tbr.append(f)
    for r in tbr:  #discard files to be discarded avoiding errors if there are repetitions in tbr
        try:
            flist.remove(r)
        except:
            pass

    if len(flist) == 0:
        print "No files with root '%s'." %args[0]
        sys.exit(20)
    else:
        return flist


def read_cut(fname, frac=0.5, verbose=False):
    """It reads the table in file 'fname' and store it into a numpy array.
    Then returns the array with the first fraction 'frac' cut.
    'frac' must be [0,1] (no check for this in the function)
    if 'verbose == True' print out more information """
    if(verbose == True):
        print "reading file %s" % i
    chain = np.loadtxt(fname)   #read the ith chain

    if(verbose == True):
        print "cutting file %s" % i
    discarded = np.sum(chain[:,0])*frac  #count the number steps of the chain and optain the number of steps to discard
    for j in np.arange(chain.shape[0]):  #get the number of lines to be discarded 
        part_sum = np.sum(chain[:j,0])
        if(part_sum > discarded):
            break
    return chain[j:,:]

def read_paramnames(fname, verbose = False):
    """
    Read the parameter names list from file 'fname' and return a list of strings.

    Parameters
    ----------
    fname: string
        input file name
    verbose: bool (optional)
        verbose mode

    output: list of strings
        list of strings containing the lines of the parameter name file
    """
    if(verbose == True):
        print "Reading the file with the parameters' names: %s." % fname
    pf = open(fname, 'r')   #open the file
    paramnames = pf.readlines()    #read saving the lines as a list of strings
    pf.close()              #close
    return paramnames

def printGR(fname, conv, paramnames):
    """
    Save the list of R values for each parameter, computed with the Gelman & Rubin criterion.
    
    Parameters
    ----------
    fname: string
        output file name
    conv: list
        list of R values computed
    paramnames: list of string
        list of strings containing the parameter names with the format "name \t latex code"
    """

    out = open(fname, 'w')    #open file
    out.write("R values from the Gelman and Rubin criterion:\n")   #Header
    for c,p in zip(conv,paramnames):
        if(np.isnan(c) == False):
            out.write("%s: R-1 = %3.2e \n" %(p.strip().split("\t")[1],c-1))

    out.write("\n Best and worse R-1 values: %3.2e, %3.2e\n" %(np.nanmin(conv)-1, np.nanmax(conv)-1))
    out.close()
    pass


def marg(chain, sd=True, percentile=None, verbose=False):
    """ 
    Given a MCMC chain 'chain' and a list of percentiles, return the mean, 
    the standard devition and the values corresponding to the percentiles required

    Parameters
    ----------
    chain: 2D array
        array containing the chain with the following format: weight, likelihood, parameters
    sd: bool (optional)
        if True compute the standard deviation, if False, doesn't
    percentile: list (optional)
        if not None compute the percentiles provided 
    verbose: bool (optional)
        verbose mode

    output: 1D and 2D array
        mean, stddev (if computed) and scores at the given percentiles (if computed)
    """

    if(verbose == True):
        print "Starting with the 1D marginalization."

    mean = np.average(chain[:,2:], weights=chain[:,0], axis=0)  #mean of all the columns containing parameters
    if(sd == True):
        stddev = ms.stddev(chain[:,2:], weights=chain[:,0], axis=0)      #stddev of all the columns containing parameters
    if(verbose == True):
        print "Mean and standard deviation computed."

    if(percentile != None):
        tchain = ms.append_weights(chain[:,2:], chain[:,0])   #replicate the rows of the total chain according to the weights (done only once). WARNING: the first two columns are discarded

        perc = (100. - np.array(percentile))/2.    # convert perc to an array containing the lower and higher franction corresponding the all the percentiles ordered from smaller to larger
        perc = np.append(perc, 100.-perc)
        params_perc = ms.percentile(tchain, perc=perc)
        if(verbose == True):
            print "Percentiles computed."

    if(sd == True and percentile != None):
        return mean, stddev, params_perc
    elif(sd == True and percentile == None):
        return mean, stddev
    elif(sd == False and percentile != None):
        return mean, params_perc
    else:
        return mean


def print_marg(fname, mean, stddev, percentile, scores, paramnames):
    """
    Writes a file 'fname' with a latex table containing the mean, 
    standard deviation and percentile scores for the parameters

    Parameters
    ----------
    fname: string
        output file name
    mean: 1D array
        mean of the parameters
    stddev: 1D array
        standard deviation of the parameters
    percentile: list
        list of the percentiles to be considered
    scores: 2D array
        scores of the percentiles. The first column are the percentiles to consider
        For each percentile, the two scores at (percentile/2.) and (100 - percentile/2.) must be listed
    paramnames: list of strings
        list of strings containing the parameter names with the format "name \t latex code"
    """

    list_perc = range(len(percentile)-1,-1,-1)  #create a reversed iterator to write the percentile scores

    out = open(fname, 'w')    #open file
    out.write("\\begin{table} \n  \\centering \n  \\begin{minipage}{80mm} \n")   #LaTex table
    out.write("    \\caption{Table from fileroot \verb=%s=} \\label{tab:%s} \n" % (os.path.splitext(fname)[0], os.path.basename(os.path.splitext(fname)[0])) )
    out.write("    \\begin{tabular}{ c c ")   # table type
    for i in percentile:
        out.write("c ")
    out.write("} \n    \\hline \n")
    out.write("      parameter name \t & \t mean $\pm$ stddev")   #write the header
    for i in percentile:
        out.write("\t & \t mean $\pm %3.1f\\%%$" % i)
    out.write("\\\\ \n      \\hline \n")
    for i,p,m,s in it.izip(it.count(), paramnames, mean, stddev):  #write the parameters 'mean+-perc'
        if(np.abs(s)>1e-7 or np.sum(np.abs(scores[:,i+1]-m)>1e-7)>0):
            out.write("      $%s$" % p.strip().split("\t")[1])   #parameter name
            out.write("\t & \t $%3.2e \pm %3.2e$" % (m, s))
            for j in list_perc:
		out.write("\t & \t $%3.2e^{%3.2e}_{%3.2e}$" % (m, scores[j,i+1]-m, scores[-j-1,i+1]-m) )
            out.write("\\\\ \n")
    out.write("      \\hline \n    \\end{tabular} \n  \\end{minipage} \n\\end{table}")   #LaTex table

    out.close()  #close file

    pass


if __name__ == "__main__":   # if is the main
    """This routine reads all the files from a COSMOMC run, throw away the first part,
    put them together and save the result in a new file"""

    (options, args) = options(op.OptionParser(version="%prog version 1.1"))   #create the object optparse

    options = check_optarg(options, args)   #check if options and arguments are ok

    flist = get_file_names(options, args)  #get the file list

    if(options.outdir != None):    #if the outdir is given substitute the directory
        outroot = os.path.join(options.outdir, os.path.basename(args[0]))
    else:
        outroot = args[0]
    outfile = outroot+"."+str(options.fraction)+"_total."+options.ext  # outfile name for the total chain
    if(options.marg == True):
        outfile1D = outroot+"."+str(options.fraction)+".margestats" # outfile for the 1D marginalized 
    if(options.conv == True):
        outfileGR = outroot+"."+str(options.fraction)+".conv" # outfile for the results of the convergency test

    if(options.verbose == True):
        print "The program will work on the following files %s." % flist

    tot_chain = [] #initialise list that contains the arrays from the various input files
    for i in flist:
        tot_chain.append(read_cut(i, frac=options.fraction, verbose=options.verbose))

    if(options.conv == True or options.marg == True):
        paramnames = read_paramnames(args[0]+".paramnames", verbose = options.verbose)  #read the parameter names
        if(len(paramnames) != tot_chain[0].shape[1]-2):  #check if there are errors
            print "The number of parameters in the chain and in the parameter names file do not correspond." 
            sys.exit(100)

    if(options.conv == True):
        convlist = ms.GR_criterion(tot_chain)
        if(options.verbose == True):
            print "Saving the values from the Gelman & Rubin convergence criterion to file '%s'." %outfileGR
        printGR(outfileGR, convlist, paramnames)  #print the R values

    tot_chain = np.vstack(tot_chain)   #convert the list of arrays in a single array

    if(options.verbose == True):
        print "Total chain created. It will be saved in the file '%s." % outfile
    np.savetxt(outfile, tot_chain, delimiter='\t', fmt=options.fmt)   #output chain saved

    if(options.marg == True):   #1D marginalization required
        mean, stddev, params_perc = marg(tot_chain, percentile=options.percentile, verbose=options.verbose)  #do the mean, standard deviation and percentile scores

        if(options.verbose == True):
            print "Writing to the file %s the 1D statistics." % outfile1D
        print_marg(outfile1D, mean, stddev, options.percentile, params_perc, paramnames)   #print out the 1D marginalised file

    exit()

