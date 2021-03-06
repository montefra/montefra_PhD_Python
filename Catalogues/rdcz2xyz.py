#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Converts catalogues with ra, dec, redshift into positions x y z in Mpc/h assuming a cosmology

The output file has the structure
x	y	z	w	bias	n(x,y,z)	n_b(z)	M	redfhist
(The columns that are not present in the input file are filled with 1, except the redshift

"""

#import constants as const
import my_functions as mf
import numpy as np
import pandas as pd
#from memory_profiler import profile

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

    description = """ Convert the file from ra, dec, redshift into cartesian coordinates assuming a cosmology.
    The name of the output file name is derived from the input one and the modification 
    can be set with the options. ra, dec and redshift are by default assumed to be in the first three columns,
    and all the following columns are by default copied after x, y and z. If the number of columns in 
    the input file (or in the ones read) is more than 8 the exceding ones are cut.
    Redshift is copied in the last column. 
    The output file has the structure
    x	y	z	w	bias	n(x,y,z)	n_b(z)	M	redshift
    (The columns that are not present in the input file are filled with 1, except the redshift)"""
    p = ap.ArgumentParser(description=description,
            formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument("ifname", nargs='+', action=apc.file_exists(),
            help="""Input file name(s), containing ra and dec in the first two
            columns""")

    p = apc.version_verbose( p, '1' )

    p, group = apc.insert_or_replace(p)
    p, group = apc.overwrite_or_skip(p)

    p, pandas = apc.pandas_group(p)

    description="""Cosmology to use to convert ra, dec and z in distances. h0=1
    is equivalent to have the distance in Mpc/h with any other value of h0"""
    p, cosmo = apc.cosmology_group( p, description=description, h0_def=1. )
    cosmo.add_argument("--zrange", action="store", nargs=2, type=float,
            help="""Lower and upper limit for the redshift. If this option is
            given the distance is computed only ones and then interpolated in
            the values in the files""")
    cosmo.add_argument("--nbins", action="store", type=int, default='500',
            help='Number of bins in redshift.')

    p.add_argument('--negative-z', action='store', choices=[None, 'skip', 'tozero'],
            default=None, help="""If *None* no check on negative redshifts,
            otherwise either skip the corresponding lines or set to 0""")

    p.add_argument("--fmt", default="%7.6e", action=apc.StoreFmt, nargs='+',
            help="Format of the output files")

    p.add_argument("--usecols", action="store", nargs="+", type=int,
            help="""Read the selected columns. By default read all the
            columns.If thi option is used, make sure the the first three
            columns read in are ra, dec and redshift.""") 

    description = """Parameters related to the parallel computation"""
    p, parallel = apc.parallel_group( p, description=description )

    return p.parse_args(args=argv)  
#end def parse(argv)

def comoving_distance( om, ok, wde, h0, zrange=None, nbins=None):
    """Compute the comoving distance given the cosmology and returns a function
    reference. If *zrange* is not *None*, the the comoving distance is evaluated
    in *nbins* in *zrange* and then a UnivariateSpline is returned
    Parameters
    ----------
    om, ok, wde, h0: floats
        omega matter, omega curvature, dark energy equation of state and reduced
        hubble parameter (h0=1 is equivalent of getting the distance in Mpc/h) 
    zrange: 2 element list
        lower and upper limit in the redshift range
    nbins: int
        number of bins in redshift 
    output
    ------
    dis: function
        function that evaluates the comoving distance at given redshift(s)

    Examples
    --------
    d = distance( 0.27, 0, -1, 1)
    z = np.linspace(0, 1, num=50)
    comdis = d(z)
    """

    import cosmologydir as c
    cos = c.Setcosmology( om=om, ok=ok, wde=wde, h=h0 ) #set the cosmology
    d = c.Distance( cos )   #create the distance object
    if( args.zrange == None ):  #If no interpolation required
        dis = d.comoving_distance_z #reference the function that compute the distance
    else:   #otherwise compute the comoving distance 
        import scipy.interpolate as spi  #scipy interpolation routines
        z = np.linspace( zrange[0], zrange[1], num=nbins ) #find the redshifts
        D_c = d.comoving_distance_z(z)    #get the distance
        dis = spi.InterpolatedUnivariateSpline( z, D_c )   #create the spline function
    return dis   
#end def comoving_distance( ... )

def rdz2xyz(rdz, dis):
    """Get an array with ra, dec and redshift and return an array with x, y and z in Mpc/h 

    ----------
    Parameters
    rdz: 2D numpy array
        array with ra, dec, redshift
    dis: function
        function that returns the distance at given redshifts

    Output
    xyz: 2D numpy array
        comoving coordinates: x, y, z
    """
    xyz = np.empty_like(rdz)  #create the output array
    rdz[:,0] *= np.pi/180     #convert ra from deg to rad
    rdz[:,1] = np.pi/2. - rdz[:,1] * np.pi/180     #convert dec from deg to rad
    rdz[:,2] = dis( rdz[:,2] )       #transform redshift in distance in Mpc/h
    xyz[:,0] = rdz[:,2] * np.sin(rdz[:,1]) * np.cos(rdz[:,0])   #fill the x coordinates
    xyz[:,1] = rdz[:,2] * np.sin(rdz[:,1]) * np.sin(rdz[:,0])   #fill the y coordinates
    xyz[:,2] = rdz[:,2] * np.cos(rdz[:,1])    #fill the z coordinates

    return xyz   
#end def rdz2xyz(rdz, dis):

#@profile
def convert_save(f, distance, **kwargs ):
    """
    Read file *f*, converts ra, dec, z into cartesian coordinates, computing the
    comoving distance at redshift z stored in *distance*, and save to a new file
    Parameters
    ----------
    f: file object or string
        file containing ra, dec and z
    distance: function
        function that evaluates the comoving distance at given redshift(s)
    kwargs: keyword arguments
    output
    ------
    max, min: lists
        maximum and minimum values of x, y and z
        If kwargs['skip'] == True and the output file name already exists, a *None*
        is returned

    accepted kwargs that affects the function
    +verbose: verbose mode [True|False] 
    +replace: replace string *replace[0]* with *replace[1]* in f.name
    +insert: insert string *insert[0]* before *insert[1]* in f.name
    +skip: existing file names skipped [True|False]
    +overwrite: existing file names overwritten [True|False]
    +pandas: use pandas for the input
    +chunks: chunksize in pandas.read_table
    +usecols: columns to read from the input files. the first three must be ra,
        dec and redshift
    +negative_z: check or not for negative redshifts and perform action [None, 'skip', 'tozero']
    +fmt: format of the output file
    """
    if kwargs['verbose']:
        print("Processing file '{}'".format(f))
    ofile = mf.create_ofile_name(f, **kwargs) # create the output file name

    if kwargs['pandas']:
        return use_pandas(f, ofile, distance, **kwargs)
    else:
        return use_numpy(f, ofile, distance, **kwargs)

#end def convert_save(f, distance, **kwargs ):

def use_numpy(fin, fout, distance, **kwargs):
    """
    Reads the file using numpy loadtxt and save it using savetxt. The rest is like 'convert_save'
    Parameters
    ----------
    f: file object or string
        file containing ra, dec and z
    fout: string
        output file name
    distance: function
        function that evaluates the comoving distance at given redshift(s)
    kwargs: keyword arguments
    output
    ------
    max, min: lists
        maximum and minimum values of x, y and z
        If kwargs['skip'] == True and the output file name already exists, a *None*
        is returned

    used kwargs that affects the function
    +usecols: columns to read from the input files. the first three must be ra,
        dec and redshift
    +negative_z: check or not for negative redshifts and perform action [None, 'skip', 'tozero']
    +fmt: format of the output file
    """
    cat = np.loadtxt(fin, usecols=kwargs['usecols'])

    if kwargs['negative_z'] is not None:
        negz = cat[:,2] < 0
        if kwargs['negative_z'] == 'skip':
            cat = cat[~negz,:]
        else:
            cat[negz,2] = 0.

    (nrows, ncolumns) = cat.shape
    if ncolumns > 9:
        cat = cat[:,:9]
    elif ncolumns < 9:
        cat = np.hstack([cat, np.ones([nrows, 9-ncolumns])])

    cat[:,8] = cat[:,2]   #save the redshift in the last column of the output file

    cat[:,:3] = rdz2xyz(np.copy(cat[:,:3]), distance)   #convert ra, dec, red in x,y,z in Mpc/h

    # save the converted catalogue
    np.savetxt(fout, cat, fmt=kwargs['fmt'], delimiter='\t')
    return np.amax(cat[:,:3], axis=0), np.amin(cat[:,:3], axis=0)
#end def use_numpy(fin, fout, distance, **kwargs):

def use_pandas(fin, fout, distance, **kwargs):
    """
    Reads the file using pandas read_table and save it using savetxt. The rest is like 'convert_save'
    Parameters
    ----------
    f: file object or string
        file containing ra, dec and z
    fout: string
        output file name
    distance: function
        function that evaluates the comoving distance at given redshift(s)
    kwargs: keyword arguments
    output
    ------
    max, min: lists
        maximum and minimum values of x, y and z

    used kwargs that affects the function
    +pandas: use pandas for the input
    +chucks: read, elaborate and print file in chunks of size 'chunks'
    +usecols: columns to read from the input files. the first three must be ra,
        dec and redshift
    +negative_z: check or not for negative redshifts and perform action [None, 'skip', 'tozero']
    +fmt: format of the output file
    """

    if kwargs['chunks'] is None:  #read the whole file in one go
        cat = pd.read_table(fin, header=None, sep='\s', skiprows=mf.n_lines_comments(fin))
        if kwargs['usecols'] is not None:
            cat = cat[kwargs['usecols']]

        if kwargs['negative_z'] is not None:
            cat = set_negative_z_pandas(cat, kwargs['negative_z'])
        cat = create_out_pandas(cat, distance) # convert to output array
        np.savetxt(fout, cat, fmt=kwargs['fmt'], delimiter='\t')

        return np.array(cat[range(3)].max()), np.array(cat[range(3)].min())

    else:  #read the file in chuncks
        chunks = pd.read_table(fin, header=None, sep='\s',
                skiprows=mf.n_lines_comments(fin), chunksize=kwargs['chunks'])
        cmin, cmax = [],[]  #contain minimum and maximum of the 
        with open(fout, 'w') as fo: #open the output file
            for cat in chunks:  #loop over the chunks

                if kwargs['usecols'] is not None:
                    cat = cat[kwargs['usecols']]

                if kwargs['negative_z'] is not None:
                    cat = set_negative_z_pandas(cat, kwargs['negative_z'])
                cat = create_out_pandas(cat, distance) # convert to output array
                np.savetxt(fo, cat, fmt=kwargs['fmt'], delimiter='\t')

                cmin.append(np.array(cat[range(3)].min()))
                cmax.append(np.array(cat[range(3)].max()))

        return np.amax(cmax, axis=0), np.amin(cmin, axis=0)
#end def use_pandas(fin, fout, distance, **kwargs)

def set_negative_z_pandas(c, negative_z):
    """check negative z in the catalogue 'c' and to what required"""
    negz = c[2] < 0
    if negative_z == 'skip':
        c = c[~negz]
    else:
        c[2][negz] = 0.
    return c
#end def set_negative_z_pandas(c, negative_z):

def create_out_pandas(cat, distance):
    """from cat, create the output. 'distance' is used to convert ra, dec
    and redshift to cartesian x, y and z"""
    ncolumns = cat.columns.size
    nrows = cat.index.size
    # set the number of columns to 9
    if ncolumns > 9:
        cat = cat.ix[:,:9]
    elif ncolumns < 9:
        ones = np.ones(nrows)
        for i in range(ncolumns, 9):
            cat[i] = ones
    # copy redshift (column 3) into the last one
    cat[8] = cat[2]

    cat[range(3)] = pd.DataFrame(rdz2xyz(np.array(cat[range(3)]), distance)) 
    return cat
#end def create_out_pandas(cat, distance):


if __name__ == "__main__":   # if is the main

    import sys
    args = parse(sys.argv[1:])

    #compute the comoving distance for the given cosmology
    dis = comoving_distance( args.om, args.ok, args.wde, args.h0,
            zrange=args.zrange, nbins=args.nbins)
    if args.verbose:
        print("Distance vs redshift function done")

    #if parallel computation required, check that Ipython.parallel.Client 
    #is in installed and that the ipycluster has been started
    if args.parallel :
        from ipython_parallel import Load_balanced_view as Lbv
        parallel_env = Lbv()  #initialize the object with all my parallen stuff
        args.parallel = parallel_env.is_parallel_enabled()

    #run the script in serial mode
    if(args.parallel == False):  #if: parallel
        #initialise the list of maxima and minima in the output file
        maxi, mini = [], []
        for fn in args.ifname:  #file name loop
            #convert the coordinates and return maxima and minima
            temp = convert_save(fn, dis, **vars(args) ) 
            if(temp != None):
                maxi.append(temp[0])
                mini.append(temp[1])
    #run the script using the IPython parallel environment 
    else:    #if: parallel
        #execute some import on all engines
        import os
        #the absolute path and file name of this script
        path, fname = os.path.split(os.path.abspath(sys.argv[0]))
        #command to run on all the engines
        imports = ['import numpy as np', 'import my_functions as mf',
                'import pandas as pd',
                #add the script directory to the python path
                'import sys', 'sys.path.append("{0}")'.format(path),     
                #import the desired function in the namespace
                'from {0} import *'.format(os.path.splitext(fname)[0])]  
        parallel_env.exec_on_engine(imports)

        initstatus = parallel_env.get_queue_status()  #get the initial status

        #submit the jobs and save the list of jobs
        runs = [parallel_env.apply(convert_save, os.path.abspath(fn),
            dis, **vars(args)) for fn in args.ifname]

        if args.verbose :   #if some info is required
            parallel_env.advancement_jobs(runs, update=args.update,
                    init_status=initstatus)
        else:   #if no info at all is wanted
            parallel_env.wait( jobs=runs )  #wait for the end

        #get the maxima and minima from the computations excluding the None
        maxi = [r.result[0] for r in runs if r.result is not None]
        mini = [r.result[1] for r in runs if r.result is not None]

        #clear the variable in the parallel environment to avoid filling up memory
        parallel_env.clear_cache()  
    #end if: parallel

    #compute absolute maximum and minimum
    absmax = np.max( maxi, axis=0)
    absmin = np.min( mini, axis=0)

    maxstring_lenght = len("Difference")
    string_template = "{:<{}}  "+"  ".join(["{:^9}",]*3)
    float_template = "{:<{}}: "+", ".join(["{:>+9.3f}",]*3)

    print(string_template.format(" ", maxstring_lenght, "x", "y", "z"))
    print(float_template.format("Maximun", maxstring_lenght, *absmax))
    print(float_template.format("Minimum", maxstring_lenght, *absmin))
    print(float_template.format("Difference", maxstring_lenght, *absmax-absmin))

    exit()
