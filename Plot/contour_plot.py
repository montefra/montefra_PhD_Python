#!/usr/bin/python
# -*- coding: utf-8 -*-

import colmaps as cm
import contour_func as cf
import itertools as it  #optimized iteration tools
import matplotlib.pyplot as plt  #plotting stuff
import myplotmodule as mpm   # my model used to plot
import numpy as np  #numpy
import os    #contain OS dependent stuffs: it helps with sistem portability  
import sys    #mondule sys

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

    description = """Draw contourplots for a list of MCMC chains.
    File name and structure is base on CosmoMC.
    ifroot+ext-chain must have the following structure:
    'counter, likelihood, parameters'.
    ifroot+ext-paramn: are the paramname files
    """

    p = ap.ArgumentParser(description=description,
            formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument("columns", action="store", nargs=2, help="""Name of the
            two colums as appears in the first column of each parameter file
            name.""")

    p.add_argument("ifroot", action="store", nargs='+', 
            help="Input file name(s), containing z in one of the columns")

    p = apc.version_verbose(p, '2')

    # file related options
    pf = p.add_argument_group(description='Input-output file options')

    pf.add_argument("--ext-chain", action="store", default=".0.5_total.dat",
            help="""Extension to create tht chain file names""")
    pf.add_argument("--ext-parmn", action="store", default=".paramnames",
            help="""Extension of the parameter file names""")

    pf.add_argument('-o', '--output', help="""Save the contour plot into file
            '%(dest)s' instead of plotting it""")

    # plot options
    pp = p.add_argument_group(description="Plot related options")

    pp.add_argument("--figure-size", nargs=2, default=[10.,10.], type=float,
            action=apc.Cm2Inch, help="Figure size in cm")

    pp.add_argument('-x', '--xlim', type=float, nargs=2, help="""Plot x
            limits""")
    pp.add_argument('-y', '--ylim', type=float, nargs=2, help="""Plot y
            limits""")

    pp.add_argument("-a", "--alpha", action="store_true", help="""Turn on
            transparency.""")

    pp.add_argument("-w", "--line-width", action="store", type=float,
            default=2, help="Line width for the contour plots")

    pp.add_argument("--no-fill", action="store_false", help="""Don't draw
            filled contours""")

    pp.add_argument("--horizontal", type=float, 
            help="Draw a horizontal line at y='%(dest)s'")
    pp.add_argument("--vertical", type=float, 
            help="Draw a vertical line at x='%(dest)s'")
    pp.add_argument("--diagonal", type=float, nargs=2, help="""Draw a straight
            line at y='%(dest)s[0]'*x + '%(dest)s[1]'""")
    pp.add_argument("--line", type=float, nargs=3, help="""Draw a line
            y='%(dest)s[0]'*x^'%(dest)s[1]' + '%(dest)s[2]'""")

    # 2D histograms options
    p2D = p.add_argument_group(description="""Histogram creation options""")

    p2D.add_argument("-b", "--num_bins", action="store", type=int, default="80", 
            help="Number of bins for the 2D histogram") 

    p2D.add_argument("-s", "--smooth", action="store_true", help="""Turn on the
            smothing of the 2D histogram""")

    p2D.add_argument("-l", "--level", type=float, nargs="+",
            help="""Confidence levels for the contour plot. '%(dest)s' must be
            [0,100].""")

    #text and font options
    pfo = p.add_argument_group(description="Text and font options")

    pfo.add_argument("--font-size", type=float, default=20, 
            help="Axis font size")

    pfo.add_argument("-t", "--text", nargs='+', action=apc.multiple_of(3,
            reshape=True), help="""Writes the text '%(dest)s[2]' at coordinates
            x='%(dest)s[0]', y='%(dest)s[1]'. Multiple text can be given providing,
            giving triplet of x,y, and text""")
    pfo.add_argument("--text-fsize", type=float, help="""Texts font size.
            Defaults to axis font size""")
    pfo.add_argument("--text-bkgr", help="Text background")

    pfo.add_argument("--xlabel", help="Custom x label")
    pfo.add_argument("--ylabel", help="Custom y label")

    #legend options 
    pl = p.add_argument_group(description='Legend options')

    pl.add_argument("--legend", nargs='+', help="""Legend tags. If given, the
        number of elements must be the same as the number of input root and in
        the same order""")

    pl.add_argument("--loc", type=apc.int_or_str, default=0,
            help='Legend location (see matplotlib legend help)')

    pl.add_argument("--legend-fsize", type=float, help="""Legend tags font
            size. Defaults to axis font size""")

    pl.add_argument("--color-text", action="store_true", help="""Legend text
            with the same color as the lines in the contour plots.""")

    pl.add_argument("--invert-legend", action="store_true", 
            help="Invert the legend entries")

    pl.add_argument("--legend-frame", help="""Set the legend frame to the given
            color. Default no frame""")

    return p.parse_args(args=argv)  
#end def parse(argv)

#############################################
###                 Main                  ###
#############################################

if __name__ == "__main__":   # if is the main

    import sys
    args = parse(sys.argv[1:])

    # command line parameter checks
    if args.level is not None:
        args.level = np.array(args.level)/100.
        if not np.all((args.level>=0)&(args.level<=1)):
            print("Option '--level' must be between 0 and 100")
            exit()
        n_levels = args.level.size
    if args.legend is not None and len(args.legend)!=len(args.ifroot):
        print("The number of legend tags must be the same as the number of file roots")
        exit()
        
    #line style twicks uncomment for final files
    t0,t1 = mpm.linestyles[:2]   #invert the first and second elements of the line style
    mpm.linestyles[:2] = t1,t0
    #t0,t1 = mpm.linestyles[2:4]   #invert the third and fourth elements of the line style
    mpm.linestyles[2:4] = t1,t0

    # get parameter names
    paramnames = cf.get_paramnames(args.ifroot, args.columns,
            verbose=args.verbose, ext=args.ext_parmn) 

    #extract the columns to read from 'paramnames'
    usecols = [[i[0]+2 for i in pn] for pn in paramnames]
    #get the chains 
    totchain = cf.get_chains(args.ifroot, usecols, verbose=args.verbose,
            ext=args.ext_chain)

    if args.verbose:
        print("compute x and y range")
    xr, yr = [np.infty,-np.infty], [np.infty,-np.infty]
    for ch in totchain:
        xr = [min(xr[0], np.amin(ch[:,1])), max(xr[1], np.amax(ch[:,1]))]
        yr = [min(yr[0], np.amin(ch[:,2])), max(yr[1], np.amax(ch[:,2]))]

    if args.verbose:
        print "Converting the MCMC chains in 2d histograms"
    h2D, xe, ye = cf.hist2D(totchain, xr=xr, yr=yr, bins=args.num_bins,
            smooth=args.smooth)

    #select the colors and make the custom colors
    if args.level is not None:  #if the confidence levels are given
        # return a list of arrays containing the ampliudes relative to the
        # maxima corresponding to the confidence levels
        amplitude = cf.h2D2conflev(h2D, args.level)
        n_shades = 1   #if alpha shades wanted give the number of them
        if args.alpha:
            n_shades = len(args.ifroot)
        #create authomatic colors
        colors = cm.cols_shades(len(args.ifroot), n_levels, alpha=n_shades)
        #colors = cm.custom_colors()
    
    fig=plt.figure(1, figsize=args.figure_size)  #open the figure
    subpl = fig.add_subplot(1,1,1)  #create a subplots

    #do the contour plots
    lines = []
    for i,h,x,y in it.izip(it.count(), h2D, xe, ye):   #plot the filled contours plots
        if(args.level != None):  #if the confidence levels are given
            col = (tuple(colors[i,j,:]) for j in range(n_levels))
            alphas = colors[i,0,3]
            if args.no_fill:
                #filled contours
                subpl.contourf(x, y, h, colors=col, alpha=alphas,
                        levels=amplitude[i], origin='lower')
            #plot the empty contours and get the outermost line
            lines.append(subpl.contour(x, y, h,
                colors=(tuple(colors[i,-1,:3]),), levels=amplitude[i],
                linestyles=mpm.linestyles[i%mpm.numls],
                linewidths=args.line_width, origin='lower').collections[0])
        else:
            subpl.contourf(x, y, h, origin='lower')
            #plot the contours and get the outermost line
            lines.append(subpl.contour(x, y, h,
                linestyles=mpm.linestyles[i%mpm.numls],
                linewidths=args.line_width, origin='lower').collections[0])

    #plot all the extra lines
    if(args.horizontal != None):
        subpl.axhline(y=args.horizontal, color="k", linestyle=':', lw=1.5)
    if(args.vertical != None):
        subpl.axvline(x=args.vertical, color="k", linestyle=':', lw=1.5)
    if(args.diagonal != None):
        subpl.plot(x, args.diagonal[0]*x+args.diagonal[1], color="k",
                linestyle=':', lw=1.5)
    if(args.line != None):
        subpl.plot(x, args.line[0]*np.power(x,args.line[1])+args.line[2],
                color="k", linestyle=':', lw=1.5)

    # write the text
    if(args.text != None):
        if(args.text_fsize == None):
            args.text_fsize = args.font_size
        for (x,y,t) in args.text:
            txt = subpl.text(float(x),float(y),t, fontsize=args.text_fsize)
        if(args.text_bkgr != None):
            txt.set_backgroundcolor(args.text_bkgr)

    #set the axis size
    if(args.xlim != None):
        subpl.set_xlim(args.xlim)
    if(args.ylim != None):
        subpl.set_ylim(args.ylim)

    #set xlabel and font size
    if(args.xlabel == None):
        args.xlabel = "$"+paramnames[0][0][2]+"$"
    if(args.ylabel == None):
        args.ylabel = "$"+paramnames[0][1][2]+"$"
    subpl.set_xlabel(args.xlabel, fontsize=args.font_size)
    subpl.set_ylabel(args.ylabel, fontsize=args.font_size)
    for label in subpl.get_xticklabels():
        label.set_fontsize(args.font_size)
    for label in subpl.get_yticklabels():
        label.set_fontsize(args.font_size)

    #legend
    if(args.legend != None):  #write the legend if required
        if args.invert_legend:
            lines.reverse()
            args.legend.reverse()
        l = subpl.legend(lines, args.legend, args.loc, numpoints=1, borderpad=0.3,
                columnspacing=0.8, handlelength=1.4, borderaxespad=0.1,
                labelspacing=0.2)
        if(args.legend_frame is None):
            l.draw_frame(False)
        else:
            l.get_frame().set_edgecolor(args.legend_frame)
            l.get_frame().set_facecolor(args.legend_frame)

        txt = l.get_texts()
        if(args.legend_fsize == None):
            args.legend_fsize = args.font_size
        plt.setp(txt, fontsize=args.legend_fsize)
        if args.color_text:
            for i,line in enumerate(lines):
                txt[i].set_color(list(line.get_color()[0]))

    #save or display the figure
    fig.tight_layout()
    if args.output is None:
        plt.show()
    else:
        fig.savefig(args.output)

    exit()
