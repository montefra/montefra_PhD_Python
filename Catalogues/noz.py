#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Computes N(z) and n(z) from one column in a (set of) file(s) and save it (them)
for each file and/or computes the mean and save it. Weights can be applied. 
The routine can also be used to create simple histograms from single columns of files
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

    description = """
    Computes N(z) and n(z) from one column in a (set of) file(s) and save it
    (them) for each file and/or computes the mean and save it. Weights can be
    applied. 
    When the 'area' is not specified, a standard histogram is produced.
    The output files have the following structure:
        bin_center  bin_lower  bin_upper  hist (std)
    std is added only if the mean is required
    """
    p = ap.ArgumentParser(description=description,
            formatter_class=apc.RawDescrArgDefHelpFormatter)

    int_or_oper = """If integer, the given column is considered. Operations
    between columns can be performed preceding the column number with a c: e.g.
    c3+c6-1 add columns 3 and 6 and subtract 1 from the result. only after the
    histogram is created"""

    p.add_argument("column", action="store", type=apc.int_or_str, help="""Columns
            containing the redshift (or whatever variable for the histogram).
            """+int_or_oper)
    p.add_argument("ifname", nargs='+', action=apc.file_exists(),
            help="""Input file name(s)""")
    
    p, group = apc.insert_or_replace(p)
    p, group = apc.overwrite_or_skip(p)

    p = apc.version_verbose(p, '0.1')

    p.add_argument("-w", "--weight", action="store", type=apc.int_or_str, 
            help="""Weight for the histogram. Operation permitted as for 'column'
            If boolean expressions like '(c3==1)&(c6==0)', the weight is interpreted 
            as 1 if *True*, as 0 if *False*""")

    p.add_argument("-n", "--nbins", action="store", type=int, default='50', 
            help="Number of bins per histogram.")
    p.add_argument("--range", action="store", nargs=2, type=float, 
            help="""Lower and upper range of the bins. Computed from the first
            file if any of the following are required: 
                i) the mean; 
                ii) n(z) [so the area needs to be computed];
                iii) pandas and chunks
            """)

    p.add_argument("--mean", action="store", type=apc.outfile,
            help="""Save the mean in file '%(dest)s'""")
    p.add_argument("--mean-only", action="store_true", 
            help="Only the mean of the histograms is saved")

    p, pandas = apc.pandas_group(p)

    description="""
    If the area of the survey is given, the effective volume per bin is
    computed and the histogram of n(z) is returned. 'h0=1' is equivalent to
    requiring the volume to be in (Mpc/h)^3. For performance reasons, if the
    range is not given it will be computed as for when the mean is requested"""
    p, cosmo = apc.cosmology_group(p, description=description, h0_def=1)
    cosmo.add_argument("-a", "--area", type=float, action="store", 
            help="Area of the survey in steradians")

    p.add_argument("--fmt", default="%7.6e", action=apc.StoreFmt, nargs='+',
            help="Format of the output files")

    description = """Parameters related to the parallel computation"""
    p, parallel = apc.parallel_group( p, description=description )

    return p.parse_args(args=argv)

#---------
# create a dictionary of pieces of docstring
# and a decorator that substitute docstring using string.format method
hist_doc ={
    'common_params': """Parameters
    ----------
    f: string
        file containing the catalogue
    readcols: list
        columns to read from the input file
    hdata: None, int or string
        use column *int* or evaluate *string* to get the data to feed to the
        histogram. If *None* only one column has been read in""",

    'return_values':"""hist: array
        The values of the histogram. 
    bin_centers: array of dtype float 
        Return the bin centers (length(hist))
    bin_edges: array of dtype float
        Return the bin edges (length(hist)+1).""",

    'common_kwargs': """accepted kwargs that affects the function:
    +weight: weights to use when creating the histogram. Can be a int or a
        string with the operation to execute
    +verbose: verbose mode [True|False] 
    +nbins: number of bins for the histogram
    +range: lower and upper limit of the histogram
    +pandas: use pandas for the input
    +chunks: chunksize in pandas.read_table""",


    'kwargs_file': """+replace: replace string *replace[0]* with *replace[1]* in f.name
    +insert: insert string *insert[0]* before *insert[1]* in f.name
    +skip: existing file names skipped [True|False]
    +overwrite: existing file names overwritten [True|False]
    +fmt: format of the output file
    """
}
def format_docstring(dic_hd):
    def wrapper(func):
        doc = func.__doc__
        doc = doc.format(**dic_hd)
        func.__doc__ = doc
        return func
    return wrapper

#-------
#get the desired column(s) of the input catalogue
def get_colums(cat, hdata, weight, use_pandas):
    """ gets the (piece of the) catalogue and returns the desired columns 

    Parameters
    ----------
    cat: numpy ndarray or pandas dataframe
        catalogue from which to extract
    hdata: None, int or string
        use column *int* or evaluate *string* to get the data to feed to the
        histogram. If *None* only one column has been read in
    weight: weights to use when creating the histogram. Can be a int or a
        string with the operation to execute
    use_pandas: bool
        use pandas for the input

    Output
    ------
    values: 1D numpy array or pandas Series
        values to use to create the histogram
    weights: 1D numpy array or pandas Series
        values to use as weights in the histogram

    """
    if hdata is None:
        values = cat[:]
    elif isinstance(hdata, (int, long)):
        if use_pandas:
            values = cat[hdata]
        else:
            values = cat[:,hdata]
    else:
        values = eval(hdata)

    if weight is None:
        weights = None
    elif isinstance(weight, (int, long)):
        if use_pandas:
            weights = cat[weight]
        else:
            weights = cat[:,weight]
    else:
        weights = eval(weight)

    # if a boolean expression is give to estimate the weights, 
    # convert the bool array to int (True=1, False=0)
    if weights is not None and weights.dtype is np.dtype('bool'):
        weights = weights.astype(int)

    return values, weights

@format_docstring(hist_doc)
def get_range(f, readcols, hdata, use_pandas, nchunks):
    """get the minimum and maximum from the columns of the file
    {common_params}
    use_pandas: bool
        use pandas for the input
    nchunks: int
        chunksize in pandas.read_table
    output
    ------
    crange: list
        minimum and maximum in the column
    """
    if use_pandas:
        if nchunks is None:
            cat = pd.read_table(f, header=None, skiprows=mf.n_lines_comments(f),
                    sep='\s')[readcols]
            values, _ = get_colums(cat, hdata, None, use_pandas)
            return [np.min(values), np.max(values)]
        else:
            chunks = pd.read_table(f, header=None, sep='\s',
                    skiprows=mf.n_lines_comments(f), chunksize=nchunks)
            cmin, cmax = [],[]  #contain minimum and maximum of the 
            for cat in chunks:  #loop over the chunks
                values, _ = get_colums(cat[readcols], hdata, None, use_pandas)
                cmin.append(np.min(values))
                cmax.append(np.max(values))
            return [np.min(cmin), np.max(cmax)]
    else:
        cat = np.loadtxt(f, usecols=readcols)
        values, _ = get_colums(cat, hdata, None, use_pandas)
        return [np.min(values), np.max(values)]

#-------
# create and return and/or save histograms
@format_docstring(hist_doc)
def hist_return(f, readcols, hdata, **kwargs):
    """
    Read the required 'col' from file 'f' and return the histogram
    {common_params}
    output
    ------
    {return_values}

    {common_kwargs}
    """

    if kwargs['verbose']:
        print("Process file '{}'".format(f))

    use_pandas = kwargs.get('pandas', False)
    hweights = kwargs.get('weight', None)

    if use_pandas:
        nchunks = kwargs.get('chunks', None)
        if nchunks is None:
            cat = pd.read_table(f, header=None, skiprows=mf.n_lines_comments(f),
                    sep='\s')[readcols]
            values, weights = get_colums(cat, hdata, hweights, use_pandas)
            hist, bin_edges = np.histogram(values, bins=kwargs['nbins'],
                    range=kwargs['range'], weights=weights)
        else:
            chunks = pd.read_table(f, header=None, sep='\s',
                    skiprows=mf.n_lines_comments(f), chunksize=nchunks)
            for cat in chunks:  #loop over the chunks
                values, weights = get_colums(cat[readcols], hdata, hweights, use_pandas)
                temphist, bin_edges = np.histogram(values, bins=kwargs['nbins'],
                        range=kwargs['range'], weights=weights)
                #sum all the histogram as they are created
                try:
                    hist += temphist
                except NameError:
                    hist = np.copy(temphist)
    else:
        cat = np.loadtxt(f, usecols=readcols)
        values, weights = get_colums(cat, hdata, hweights, use_pandas)
        hist, bin_edges = np.histogram(values, bins=kwargs['nbins'],
                range=kwargs['range'], weights=weights)

    return hist, (bin_edges[1:]+bin_edges[:-1])/2., bin_edges 

@format_docstring(hist_doc)
def hist_save_return(f, readcols, hdata, **kwargs):
    """
    Read the required 'col' from file 'f', do the histogram 
    save it to file and return the it 
    {common_params}
    output
    ------
    {return_values}
    
    {common_kwargs}
    {kwargs_file}
    """
    ofile = mf.create_ofile_name(f, **kwargs) # create the output file name
    hist, bin_centers, bin_edges = hist_return(f, readcols, hdata, **kwargs)
    # create the output 2d array 
    outhist = np.vstack([bin_centers, bin_edges[:-1], bin_edges[1:], hist]).T
    np.savetxt(ofile, outhist, fmt=kwargs['fmt'], delimiter='\t')
    return hist, bin_centers, bin_edges

@format_docstring(hist_doc)
def hist_save(f, readcols, hdata, **kwargs):
    """
    Read the required 'col' from file 'f', do the histogram 
    save it to file 
    {common_params}
    output
    ------
    None
    
    {common_kwargs}
    {kwargs_file}
    """
    hist_save_return(f, readcols, hdata, **kwargs)
    return None

@format_docstring(hist_doc)
def hist_ndensity_return(f, readcols, hdata, volume_z, **kwargs):
    """
    Read the required 'col' from file 'f' and return the histogram of the number density
    {common_params}
    volume_z: list
        effective volume in each bin of the histogram. 
        Warning: kwargs['range'] cannot be *None* and 
        kwargs['nbins'] == len(volume_z). No check performed
    output
    ------
    {return_values}

    {common_kwargs}
    """
    hist, bin_centers, bin_edges = hist_return(f, readcols, hdata, **kwargs)
    hist = hist / volume_z
    return hist, bin_centers, bin_edges

@format_docstring(hist_doc)
def hist_ndensity_save_return(f, readcols, hdata, volume_z, **kwargs):
    """
    Read the required 'col' from file 'f', do the histogram of the number
    density, save it to file and return it
    {common_params}
    volume_z: list
        effective volume in each bin of the histogram. 
        Warning: kwargs['range'] cannot be *None* and 
        kwargs['nbins'] == len(volume_z). No check performed
    output
    ------
    {return_values}

    {common_kwargs}
    {kwargs_file}
    """
    ofile = mf.create_ofile_name(f, **kwargs) # create the output file name
    hist, bin_centers, bin_edges = hist_ndensity_return(f, readcols, hdata,
            volume_z, **kwargs)
    # create the output 2d array 
    outhist = np.vstack([bin_centers, bin_edges[:-1], bin_edges[1:], hist]).T
    np.savetxt(ofile, outhist, fmt=kwargs['fmt'], delimiter='\t')
    return hist, bin_centers, bin_edges

@format_docstring(hist_doc)
def hist_ndensity_save(f, readcols, hdata, volume_z, **kwargs):
    """
    Read the required 'col' from file 'f', do the histogram of the number
    density and save it to file
    {common_params}
    volume_z: list
        effective volume in each bin of the histogram. 
        Warning: kwargs['range'] cannot be *None* and 
        kwargs['nbins'] == len(volume_z). No check performed
    {common_params}
    output
    ------
    None
    
    {common_kwargs}
    {kwargs_file}
    """
    hist_ndensity_save_return(f, readcols, hdata, volume_z, **kwargs)
    return None


def select_hist_function(mean, mean_only, redshift_volumes):
    """
    Return the function to use to create the histogram
    Parameters
    ----------
    mean: bool
        mean of the histograms required
    mean_only: bool
        only the mean required (not every histogram)
    redshift_volumes: None or list of float
        if not None, the number density required (not the simple histogram)
    output
    ------
    func: function
        reference to the function to use
    """
    if redshift_volumes is None:  # simple histogram
        if mean is not None: # if the mean is required
            if mean_only:  # only the mean to file
                return hist_return
            else: # both single histograms and mean to files
                return hist_save_return
        else: # no mean histograms to files
            return hist_save
    else: # cosmological volume per redshift bin
        if mean is not None: # if the mean is required
            if mean_only:  # only the mean to file
                return hist_ndensity_return
            else: # both single histograms and mean to files
                return hist_ndensity_save_return
        else: # no mean histograms to files
            return hist_ndensity_save

def effective_volume(z_edges, area, cosmology):
    """
    Compute the effective volume in (Mpc/h)^3 or Mpc^3 for each redshift bin 
    Parameters
    ----------
    z_edges: 1d array
        edges of the redshift intervals
    area: float
        area of the footprint in steradians
    cosmology: list
        Omega_matter (om), Omega_k (ok), w_DE (wde) and h0 
    output
    ------
    volume: array
        effective volume per redshift bin 
    """
    # create the cosmology and distance instance
    import cosmologydir as c
    cos = c.Setcosmology(om=cosmology[0], ok=cosmology[1], wde=cosmology[2],
            h0=cosmology[3])
    dis = c.Distance(cos)
    # compute the effective volumes
    return dis.effective_volume_sr_z(z_edges) * area

def parse_operations(column, weight, use_pandas):
    """
    parse the column and weight entry, return a list of columns to read and,
    if any, the operations to perform for the entries and weights of the hists
    Parameters
    ----------
    column: int or string
        column to use of the histogram, or operations among columns to create
        the input of the histogram
    weight: int or string
        same as before, but for the weights
    use_pandas: bool
        use pandas or numpy
    output
    ------
    usecols: list of int
        list of columns to read in ordered with first the *histdata* and then
        the *histweight* columns 
    histdata: None, int or string
        None if *weigth* is None and *columns* int
        histdata=0 if *columns* int
        string to *eval* if *columns* is string 
    histweight: None, int, string
        None if *weight* is None
        int if *weight* is int with the position of the weight in *usecols*
        string to *eval* if *weight* is string 

    example:
        parse_operations(1, None, True/False) --> [1,], None, None #only one column to read: data

        parse_operations(1, 3, False) --> [1,3], 0, 1 # data and weight to read in.
            # the '0' entry in the read array is data, the '1' is weight
        parse_operations(1, 3, True) --> [1,3], 1, 3 # data and weight to read in.

        parse_operations('c1+c2-1', None, False) --> [1,2], 'cat[:,0]+cat[:,1]-1', None
        parse_operations('c1+c2-1', None, True) --> [1,2], 'cat[:,1]+cat[:,2]-1', None

        parse_operations('c1+c2-1', 3, False) --> [1,2,3], 'cat[:,0]+cat[:,1]-1', 2
            # three columns to read, the data are the operation in the string,
            # the weight the third entry in the read array
        parse_operations('c1+c2-1', 3, True) --> [1,2,3], 'cat[:,1]+cat[:,2]-1', 3

        parse_operations('c1+c2-1', 'c4*3+c5', False) --> 
                [1,2,4,5], 'cat[:,0]+cat[:,1]-1', 'cat[:,2]*3+cat[:,4]'
            # four columns to read, the data and weight are the operation in
            # the first and second strings, with respect to the read array
        parse_operations('c1+c2-1', 'c4*3+c5', True) --> 
                [1,2,4,5], 'cat[:,1]+cat[:,2]-1', 'cat[:,4]*3+cat[:,5]'
    """
    # column in input file is cat[#] if using pandas and cat[:,#] if numpy
    if use_pandas:
        cataloguestr = "cat[\\1]"
    else:
        cataloguestr = "cat[:,\\1]"

    import re
    pattern = re.compile(r"c(\d+)")  # re pattern with the columns name

    #get the desired column(s) to do the histogram
    if isinstance(column, (int, long)):
        usecols=[column,]
        if use_pandas:
            histdata = column
            after_read_cols = usecols[:] #index of the columns once read in
        else:
            histdata = 0
            after_read_cols = [histdata,]
    else:
        datamatches = pattern.findall(column)
        usecols = list(np.unique([int(m) for m in datamatches])) # unique indeces to read
        if use_pandas:
            after_read_cols = usecols[:]
        if not use_pandas:
            after_read_cols = range(len(usecols))
            # replace the columns required with the ones read in
            for o,n in zip(usecols, after_read_cols):
                column = column.replace('c'+str(o), 'c'+str(n))
        # change the from the column number in file to the column number in the 
        # read numpy array or pandas DataFrame
        histdata = pattern.sub(cataloguestr, column)

        
    # get the desired column(s) to do the weights
    if weight is None:
        histweight = None
        if isinstance(histdata, (int, long)):
            histdata = None
    elif isinstance(weight, (int, long)):
        #if the weight column has already been included into the data
        if weight in usecols:
            windex = usecols.index(weight)
            if use_pandas:
                histweight = usecols[windex]
            else:
                histweight = after_read_cols[windex]
        else:  #new colum
            usecols.append(weight)
            if use_pandas:
                histweight = weight
            else:
                histweight = after_read_cols[-1]+1
    else:
        weightmatches = pattern.findall(weight)
        wusecols = list(np.unique([int(m) for m in weightmatches])) # unique indeces to read
        #check if columns already included ore not
        if use_pandas:  #here simply add the new ones
            for wuc in wusecols:
                if wuc not in usecols:
                    usecols.append(wuc)
        else: #if numpy used is a bit more complex
            for wuc in wusecols:
                if wuc in usecols:
                    windex = usecols.index(wuc)
                    tempw= after_read_cols[windex]
                else:
                    usecols.append(wuc)
                    after_read_cols.append(after_read_cols[-1]+1)
                    tempw = after_read_cols[-1]
                #substitute the input column with the one in the numpy array
                weight = weight.replace('c'+str(wuc), 'c'+str(tempw))

        # change the from the column number in file to the column number in the 
        # read numpy array or pandas DataFrame
        histweight = pattern.sub(cataloguestr, weight)

    return usecols, histdata, histweight


if __name__ == "__main__":   #if it's the main

    import sys
    args = parse(sys.argv[1:])
    
    #if parallel computation required, check that Ipython.parallel.Client 
    #is in installed and that the ipycluster has been started
    if args.parallel:
        from ipython_parallel import Load_balanced_view as Lbv
        parallel_env = Lbv()  #initialize the object with all my parallen stuff
        args.parallel = parallel_env.is_parallel_enabled()

    # parse the column and weight arguments and return the list of columns
    # to read, the index or string to create the data and weight for the 
    # histograms
    usecols, histdata, args.weight= parse_operations(args.column, args.weight, args.pandas)

    # if the mean or the volume of the redshift bins is required 
    # get the range from the first file, if not passed in the command line
    if args.range is None:
        if args.mean or args.area is not None or (args.pandas and args.chunks is not None):
            if args.verbose:
                print("Extract the range from the first file")
            args.range = get_range(args.ifname[0], usecols, histdata, args.pandas, args.chunks)

    # if the number density required, compute the effective volume for each
    # redshift bin
    if args.area is not None:
        if args.verbose:
            print("compute the effective volume per redshift bin")
        z_edges = np.linspace(args.range[0], args.range[-1], num=args.nbins+1)
        redshift_volumes = effective_volume(z_edges, args.area, 
                [args.om, args.ok, args.wde, args.h0])
    else:
        redshift_volumes = None

    # chose which function to use to create the histograms
    do_hist = select_hist_function(args.mean, args.mean_only, redshift_volumes)

    if(args.parallel == False):  # if: parallel
        if args.area is None:
            hists = [do_hist(fn, usecols, histdata, **vars(args)) for fn in args.ifname]
        else:
            hists = [do_hist(fn, usecols, histdata, redshift_volumes, **vars(args))
                    for fn in args.ifname]
    else: # else: parallel
        import os
        #the absolute path and file name of this script
        path, fname = os.path.split(os.path.abspath(sys.argv[0]))
        #command to run on all the engines
        imports = ['import numpy as np', 'import my_functions as mf',
#                'import pandas as pd', 
                #add the script directory to the python path
                'import sys', 'sys.path.append("{0}")'.format(path),     
                #import the desired function in the namespace
                'from {0} import *'.format(os.path.splitext(fname)[0])]  
        parallel_env.exec_on_engine(imports)

        initstatus = parallel_env.get_queue_status()  #get the initial status

        #submit the jobs and save the list of jobs
        if args.area is None:
            runs = [parallel_env.apply(do_hist, os.path.abspath(fn), usecols,
                histdata, **vars(args)) for fn in args.ifname]
        else:
            runs = [parallel_env.apply(do_hist, os.path.abspath(fn), usecols,
                histdata, redshift_volumes, **vars(args)) for fn in
                args.ifname]

        if args.verbose :   #if some info is required
            parallel_env.advancement_jobs(runs, update=args.update,
                    init_status=initstatus)
        else:   #if no info at all is wanted
            parallel_env.wait(jobs=runs)  #wait for the end

        #get the result
        hists = [r.result for r in runs]
        #clear the variable in the parallel environment to avoid filling up memory
        parallel_env.clear_cache()  
    #end if: parallel

    if args.mean is not None:
        if args.verbose:
            print("compute and save the mean")
        # copy bin centers and egdes out of the list of lists
        bin_centers = hists[0][1][:]
        bin_edges = hists[0][2][:]
        hists = [h[0] for h in hists] # replace the hists with the histograms alone

        mean_hists = np.mean(hists, axis=0)
        sigma_hist = np.std(hists, axis=0)
        outhist = np.vstack([bin_centers, bin_edges[:-1], bin_edges[1:],
            mean_hists, sigma_hist]).T
        np.savetxt(args.mean, outhist, fmt=args.fmt, delimiter='\t')

    exit()
