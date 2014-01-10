#!/usr/bin/python3
# -*- coding: utf-8 -*-

import itertools as it
import my_statistic as ms
import numpy as np    # import numpy
import os    #contain OS dependent stuffs: it helps with sistem portability  
import sys    #mondule sys
from warnings import warn

class MyWarning(UserWarning):
    pass

def parse(argv):
    """
    This function accept a list of strings, create and fill a parser instance 
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
    import textwrap as tw

    description = tw.dedent("""\
    Given a list of file roots and of corresponding tags, reads the parameter
    name files and the chain files. Then computes the mean and either the
    standard deviation (default) or the given percentile of the parameters for
    each chain (chains assumed to have the following structure:
    chain[:,0]=weight; chain[:,1]=likelihood]; chain[:,2:]=parameters).
    The mean and the stddev or percentile are printed in a latex table with the following format:
    --------------------------------------------------------------------------------------
                       tag1                   tag2          ...        tagn
    --------------------------------------------------------------------------------------
    parameter1    m1+-s1(p1)(file1)    m1+-s1(p1)(file2)    ...    m1+-s1(p1)(filen)
    parameter2    m2+-s2(p2)(file1)    m2+-s2(p2)(file2)    ...    m2+-s2(p2)(filen)
    ...                 ...                  ...            ...             ...
    parameterm    mm+-sm(pm)(file1)    mm+-sm(pm)(file2)    ...    mm+-sm(p1)(filen)
    --------------------------------------------------------------------------------------
    """)

    p = ap.ArgumentParser(description=description,
            formatter_class=apc.RawDescrArgDefHelpFormatter)

    p.add_argument("froottag", metavar="froot-tag", nargs='+', action=apc.multiple_of(2,
            reshape=True), help="Input file name(s) and header tags. They must be given in pairs")

    p = apc.version_verbose(p, '2')

    # columns 
    p.add_argument("-c", "--columns", action="store", nargs='+', 
            help="""Name of the columns as they appear in the first column of
            each parameter file name.""")

    class What(ap.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if values[0] not in ['s', 'p']:
                msg = "Only 's' or 'p' accecpted as first argument"
                raise ap.ArgumentTypeError(msg)
            else:
                try:
                    values[1] = float(values[1])
                except ValueError:
                    msg = "The second argument must be a float"
                    raise ap.ArgumentTypeError(msg)
            setattr(namespace, self.dest, values)
    p.add_argument("-w", "--what", nargs=2, default=["s", 1.], action=What,
            help="""Decide what to do: '%(dest)s[0]= s' (standard deviation) or 'p'
            (percentile); if 's' the '%(dest)s[1]*sigma' values saved, if 'p' the
            '%(dest)s[1] percentile values around the mean computed and saved.""")

    p.add_argument("-r", "--rescale", nargs='+', default=[], action=apc.multiple_of(2,
            reshape=True), help="""Rescale the values of variable '%(dest)s[0]' by
            '%(dest)s[1]'. The same rescaling factor is shown in the parameter label.  The
            variable is identified with the short name from the parameter file.  Multiple
            rescaling can be drawn providing couple of variable-rescaling.""")

    p.add_argument("--comment", action="store_true", help="""Comment rows of
            the matrix when all the elements have null variance""")

    # file related options
    pf = p.add_argument_group(description='Input-output file options')

    pf.add_argument("--ext-chain", action="store", default=".0.5_total.dat",
            help="""Extension to create the chain file names""")
    pf.add_argument("--ext-parmn", action="store", default=".paramnames",
            help="""Extension of the parameter file names""")

    pf.add_argument("-o", "--output", default="latextable.tex", 
            help="Output file name")

    pf.add_argument("-f", "--format", default="%7.6f",
            help="Formatter for the numbers in the output table.")

    return p.parse_args(args=argv)

def do_std_perc(what):
    """
    find if standard deviation or percentile needs to be computed
    Parameters
    ----------
    what: list
        what[0] = ['s', 'p'], what[1] = float
    output
    ------
    do_std: bool
        do the standard deviation
    perc_lev: float
        percentile score
    """
    if what[0] == 's':
        return True, None
    else:
        return False, what[1]

def do_mean_std(froots, **kwargs):
    """
    Read the files and compute the mean and standard deviation/percentile
    Parameters
    ----------
    froots: list
        list of file roots
    kwargs: dictionary
        +verbose: bool
        +columns: list or None
            name of the columns to read
        +ext_paramn: string
            extension of the parameter name file
        +ext_chain: string
            extension of the chains file
        +what: list
            do the standard deviation or the percentile
    output
    ------
    keylist: list
        ordered lists of parameter short names
    longnames: dictionary
        key: parameter short name; value: long name
    mean: dictionary
        key: parameter short name; value: list of mean values with 
        len(mean[k]) == len(froots)
    std_perc: dictionary
        key: parameter short name; value: list of standard deviation or
        percentile with len(mean[k]) == len(froots)
    """
    import contour_func as cf   #functions used in contour plot
    import unify_chain as uc

    # if froots is a string convert to a list
    is_string, _ = cf._check_type(froots)
    if is_string:
        froots = [froots]
    paramnames = cf.get_paramnames(froots, params=kwargs.get('columns'),
            ext=kwargs.get('ext_parmn'), verbose=kwargs.get('verbose', False))

    # get the column numbers and the chains
    cols = [[j[0]+2 for j in i] for i in paramnames]
    chains = cf.get_chains(froots, cols, ext=kwargs.get('ext_chain'),
            verbose=kwargs.get('verbose', False))

    # create a list of keywords and the dictionary of long names
    keylist = [j[1] for j in paramnames[0]]
    longnames = {j[1]: j[2] for j in paramnames[0]}
    for pn in paramnames[1:]: # loop through the remaining paramnames
        for i,p in enumerate(pn):
            if p[1] not in keylist:
                prev_ind = keylist.index(keylist[i-1])
                keylist.insert(prev_ind+1, p[1])
                longnames[p[1]] = p[2]

    do_std, perc_lev = do_std_perc(kwargs['what'])
    # compute mean and std or percentile
    mean = {k:[] for k in keylist}
    std_perc = {k:[] for k in keylist}
    for pn, ch in zip(paramnames, chains): # loop the files
        # compute mean and stddev or percentile
        lmean, lstd_perc = uc.marg(ch, sd=do_std, percentile=perc_lev,
                verbose=kwargs.get('verbose', False), has_like=False)
        loc_kl = [j[1] for j in pn]  # list of short names of the current loop
        # loop the general key list
        for k in keylist:  
            try: # if k is in the list
                kind = loc_kl.index(k)
                mean[k].append(lmean[kind])
                if do_std:
                    std_perc[k].append(lstd_perc[kind])
                else:
                    std_perc[k].append(lstd_perc[:,kind+1]-lmean[kind])
            except ValueError:
                mean[k].append(np.nan)
                std_perc[k].append(np.nan)

    return keylist, longnames, mean, std_perc

def do_rescale(mean, std_perc, labels, rescale):
    """
    Rescale mean and standard deviation 
    Parameters
    ----------
    mean: dictionary
        key: parameter short name; value: list of mean values with 
        len(mean[k]) == len(froots)
    std_perc: dictionary
        key: parameter short name; value: list of standard deviation or
        percentile with len(mean[k]) == len(froots)
    labels: dictionary
        key: parameter short name; value: labels for the table rows
    rescale: list of len(2) lists 
        new_labels[i] = [key, rescale]
    output
    ------
    mean, std_perc, labels: dictionary
        input with the variables rescaled
    """
    for (k, r) in rescale:
        if k not in mean:
            warn("Key '{}' is not in mean and std dictionaries".format(k), MyWarning)
        else:
            r = float(r)
            mean[k] = [m*r for m in mean[k]]
            std_perc[k] = [sp*r for sp in std_perc[k]]
            if(r.is_integer()==True):
                strr = str(int(r))
            else:
                strr = str(r)
            labels[k] = strr+r'\times '+labels[k] 
    return mean, std_perc, labels

def make_table(mean, std_perc, labels, **kwargs):
    """
    create the table 
    Parameters
    ----------
    mean: dictionary
        key: parameter short name; value: list of mean values with 
        len(mean[k]) == len(froots)
    std_perc: dictionary
        key: parameter short name; value: list of standard deviation or
        percentile with len(mean[k]) == len(froots)
    labels: dictionary
        key: parameter short name; value: long name
    kwargs: dictionary
        +verbose: bool
        +what: list
            do the standard deviation or the percentile
        +comment: bool
            don't show empty lines

    output
    ------
    table: dictionary
        key: parameter short name; value: line to save in the table
    """
    if kwargs.get('verbose', False):
        print("Making the table")
    do_std, _ = do_std_perc(kwargs['what'])

    # cell input
    from string import Template
    if do_std:
        cell = Template(r"$${0:$fmt} \pm {1:$fmt}$$")
    else:
        cell = Template(r"$${0:$fmt}_{{{1[0]:+$fmt}}}^{{{1[1]:+$fmt}}}$$")
    cell = cell.substitute(fmt=kwargs['format'].strip('%'))

    table = {}
    for k,l in labels.items():
        line = ["$"+l+"$"]
        for m,s in zip(mean[k], std_perc[k]):
            if np.isnan(m): # if the mean is nan, the variable does not exist
                line.append("--")
            else:
                line.append(cell.format(m, s))
        if line.count("--") != len(mean[k]) or kwargs.get('comment', False):
            table[k] = " & ".join(line)+r"\\[0.5mm]"

    return table

def save_table(fname, tags, keylist, table):
    """
    Save the table into file fname
    Parameters
    ----------
    fname: string
        file name of the latex table
    tags: list
        tags of the column headers
    keylist: list
        list of keys of the table. Used to print the table ordered
    table: dictionary
        key: parameter short name; value: line to save in the table
    """
    
    otable  = [r'\begin{table*}']
    otable += [r'  \centering']
    otable += [r'  \begin{minipage}{160mm}']
    otable += [r'    \caption{<+Caption text+>}']
    otable += [r'    \label{tab:<+label+>}']
    otable += [r'    \begin{tabular}{l'+'c'*len(tags)+'}']
    otable += [r'      \hline']
    otable += [r'       & '+' & '.join(tags)+r" \\"]
    otable += [r'      \hline']

    for k in keylist:
        otable += ['      '+table[k]]

    otable += [r'      \hline']
    otable += [r'    \end{tabular}']
    otable += [r'  \end{minipage}']
    otable += [r'\end{table*}']

    with open(fname, 'w') as f:
        f.writelines('\n'.join(otable))


#############################################
###                 Main                  ###
#############################################

def main(argv):
    args = parse(argv)

    # separate the file roots from the tags
    args.froot, args.tag = [], []
    for ft in args.froottag: 
        args.froot.append(ft[0])
        args.tag.append(ft[1])

    # compute the mean and standard deviation or percentile
    keylist, longnames, mean, std_perc = do_mean_std(args.froot, **vars(args))

    # rescale some value
    mean, std_perc, longnames = do_rescale(mean, std_perc, longnames, args.rescale)

    # create the latex table
    dtable = make_table(mean, std_perc, longnames, **vars(args))
    # and write it to file
    save_table(args.output, args.tag, keylist, dtable) 


if __name__ == "__main__":   # if is the main

    import sys
    main(sys.argv[1:])
    exit()
