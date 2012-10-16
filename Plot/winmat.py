#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Plot rows of the window matrix"""

import itertools as it
import matplotlib.pyplot as plt
import myplotmodule as mpm
import numpy as np
import sys

ls = it.cycle( ['-' , '--' , '-.' , ':'] )

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
  
  import argparse as ap  #import optsparse: allows nice command line option handling
  import argparse_custom as apc

  description = """Plot the rows of the window matrices """
  p = ap.ArgumentParser(description=description, formatter_class=ap.ArgumentDefaultsHelpFormatter)

  p.add_argument("ifroot", action="store", nargs='+', 
      help="Root of the file name(s) containing the window matrix and the values of k for the plots")

  p = apc.version_verbose( p, '1' )

  p.add_argument("-w", "--win-ext", default=".Wij.dat", action="store", help="""Extension for the window matrix file with dimension NxM. 
      The full name must be 'ifroot+%(dest)s'.""")
  p.add_argument("--k-ext", default=".kj.dat", action="store", help="""Extension for the file with the k with dimension M. 
      The full name must be 'ifroot+%(dest)s'.""")

  p.add_argument("-k", "--krange", nargs=2, default=[-1,-1], type=float, action="store", help="""Range of k to plot. 
      If -1 minimum and/or maximum from the file of k.""")

  p.add_argument("--rows", nargs=3, default=[-1,-1,1], type=int, action="store", help="""Rows of the window matrix to plot: 
      begin, end and stride. If %(dest)s[0]<0, set to 0; if %(dest)s[1]<0 or >N, set to 0.""")
  p.add_argument("--shift", type=int, default=0, action="store", help="""Shift the starting point of the rows for 
      each of the input window matrices""")

  p.add_argument("--lw", type=float, default=1.7, action="store", help="Line width")

  p.add_argument("-l", "--legend", action="store", nargs="+", help="""Enable the legene in position '%(dest)s[0]'
      with tags '%(dest)s[1:]'. Must be len(%(dest)s)=len(ifroot)+1.""")

  p.add_argument("--cmap", action=apc.store_cycle, default="Paired", nargs="+", help="""Colormap(s) to use cyclically. 
      http://matplotlib.sourceforge.net/examples/pylab_examples/show_colormaps.html """)

  p.add_argument("-y", "--y-lim", action="store", type=float, nargs=2, default=[0,0.3], help="Limits on the y axis")

  p.add_argument("-o", "--ofile", action="store", type=ap.FileType('w'), help="Output file name. If not given, plot shown.")

  return p.parse_args(args=argv)

if __name__ == "__main__":   # if is the main

  args = parse(sys.argv[1:])

  if( args.legend and len(args.legend) != len(args.ifroot)+1):
    print( """The number of tags for the legend must be the same as the number of input files""" )

  if not isinstance(args.cmap, it.cycle) :
    args.cmap = it.cycle( args.cmap )
  fig = plt.figure(1)   #figure
  ax = fig.add_subplot(111)   #axis

  lines=[]
  for i, fr in enumerate(args.ifroot):
    if( args.verbose ):
      print("Read and plot the window matrix from file {0}{1}".format(fr, args.win_ext))
    Wij = np.loadtxt( "{0}{1}".format(fr, args.win_ext) )
    kj = np.loadtxt( "{0}{1}".format(fr, args.k_ext) )

    #set the range in k used for the plot
    krange = args.krange[:]
    if( krange[0] < 0 ): 
      krange[0] = 0
    else:
      krange[0] = (kj<krange[0]).sum()
    if( krange[1] < 0 ):
      krange[1] = kj.size
    else:
      krange[1] = (kj<=krange[1]).sum()

    #set the range of the rows to plot
    rows = np.array( args.rows, dtype=int )  #convertion to numpy array
    if( rows[0] <0 ):
      rows[0] = 0
    if( rows[1] <0 or rows[1] >Wij.shape[0] ):
      rows[1] = Wij.shape[0]

    kj = kj[ krange[0]:krange[1] ]
    ax, l2 = mpm.plot_winmat_rows(ax, Wij[:, krange[0]:krange[1] ], kj, ext=rows[:2]+i*args.shift, cmap=args.cmap.next(), 
	spacing=rows[2], lw=args.lw, ls=ls.next(), legend=True)
    lines.append(l2)

  if( args.legend ):
    ax.legend( lines, args.legend[1:], loc=args.legend[0], frameon=False) 

  ax.set_xlim( kj[0], kj[-1])
  ax.set_ylim( args.y_lim )

  ax.set_xlabel( 'k [h/Mpc]' )
  ax.set_ylabel( r'$W_{ij}$' )

  plt.tight_layout()

  if( args.ofile == None ):
    plt.show()
  else:
    fig.savefig( args.ofile )

     
  exit()

