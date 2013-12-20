#!/usr/bin/python3
# -*- coding: utf-8 -*-

import contour_func as cf   #functions used in contour plot
import itertools as it  #optimized iteration tools
import matplotlib.pyplot as plt  #plotting stuff
import myplotmodule as mpm   # my model used to plot
import numpy as np  #numpy
from warnings import warn

class MyException(Exception):
    pass
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

    description = """Draw 1D histograms for parameters from a list of MCMC chains.
    File name and structure is base on CosmoMC.
    ifroot+ext-chain must have the following structure:
    'counter, likelihood, parameters'.
    ifroot+ext-paramn: are the paramname files
    """

    p = ap.ArgumentParser(description=description,
            formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument("ifroot", action="store", nargs='+', 
            help="Input file name(s)")

    p = apc.version_verbose(p, '2')

    # columns 
    p.add_argument("-c", "--columns", action="store", nargs='+', 
            default=['omegam', 'omegal', 'H0'],
            help="""Name of the columns as they appear in the first column of
            each parameter file name.""")


    # file related options
    pf = p.add_argument_group(description='Input-output file options')

    pf.add_argument("--ext-chain", action="store", default=".0.5_total.dat",
            help="""Extension to create the chain file names""")
    pf.add_argument("--ext-parmn", action="store", default=".paramnames",
            help="""Extension of the parameter file names""")

    pf.add_argument('-o', '--output', help="""Save the contour plot into file
            '%(dest)s' instead of plotting it""")

    # plot options
    pp = p.add_argument_group(description="Plot related options")

    pp.add_argument("--set-MNRAS", action='store_true', 
            help='Use MNRAS presets for plotting')

    pp.add_argument("--figure-figsize", nargs=2, default=[10.,10.], type=float,
            action=apc.Cm2Inch, help="Figure size in cm")

    pp.add_argument("-n", "--n-cols", type=int, default=3, 
            help='Number of subplot columns')

    pp.add_argument("--bounds", action="store", type=float, nargs=4,
            default=[0, 0, 1, 1], help="""(left, bottom, right, top) in the
            normalized figure coordinate passed to 'plt.tight_layout'""")

    pp.add_argument("--x-lim", nargs='+', action=apc.multiple_of(3, reshape=True), 
            default=[], help="""Set x axis range in subplot '%(dest)s[0]' to
            '%(dest)s[1:3]'.  The subplot is identified with the short name
            from the parameter file.""")

    pp.add_argument("--x-label", nargs='+', action=apc.multiple_of(2, reshape=True), 
            default=[], help="""Set x axis label in subplot '%(dest)s[0]' to
            '%(dest)s[1]'.  The sublot is identified with the short name from
            the parameter file. Multiple labels can be drawn providing couple
            of subplot-label.""")
    pp.add_argument("--y-label", default="$P/P_{max}$", help='x axis label')

    pp.add_argument("-w", "--lines-linewidth", action="store", type=float,
            help="Line width for the contour plots")

    pp.add_argument("--vertical", nargs='+', action=apc.multiple_of(2,
        reshape=True), default=[], help="""Plot a vertical line in subplot
        '%(dest)s[0]' at x='%(dest)s[1]'. The subplot is identified with the
        short name from the parameter file. Multiple lines can be drawn
        providing couple of subplot-line.""") 

    # 1D histograms options
    p1D = p.add_argument_group(description="""Histogram creation options""")

    p1D.add_argument("-b", "--num-bins", action="store", type=int, default="80", 
            help="Number of bins for the 2D histogram") 

    #text and font options
    pfo = p.add_argument_group(description="Text and font options")

    pfo.add_argument("--axes-labelsize", type=float, help="Axis font size")

    pfo.add_argument("-t", "--text", nargs='+', action=apc.multiple_of(4,
            reshape=True), default=[], help="""Writes the text '%(dest)s[3]' at
            coordinates x='%(dest)s[1]', y='%(dest)s[2]' in axes '%(dest)s[0]'.
            Multiple text can be given providing, giving triplet of x,y, and
            text""")
    pfo.add_argument("--font-size", type=float, help="""Texts font size.""")
    pfo.add_argument("--text-bkgr", help="Text background")


    #legend options 
    pl = p.add_argument_group(description='Legend options')

    pl.add_argument("--legend", nargs='+', default=[], help="""Legend tags. If
            given, the number of elements must be the same as the number of input root
            and in the same order""")

    pl.add_argument("--legend-plot", action="store", 
            help="""Put the legend in plot '%(dest)s'. The subplot is identified
            with the short name from the parameter file. If '%(dest)s' is not
            in the list of parameters, the legend is drawn in the first empty
            subplot, if any, otherwise the figure legend is drawn. If not
            given, the legend is drawn in a plot (depends on the dictionary
            hashing)""")
    pl.add_argument("--loc", type=apc.int_or_str, default=0,
            help='Legend location (see matplotlib legend help)')

    pl.add_argument("--legend-fontsize", type=float, help="""Legend tags font
            size.""")

    pl.add_argument("--color-text", action="store_true", help="""Legend text
            with the same color as the lines in the contour plots.""")

    pl.add_argument("--invert-legend", action="store_true", 
            help="Invert the legend entries")

    pl.add_argument("--legend-frameon", action="store_true", default=False,
            help='Legend frame')

    return p.parse_args(args=argv)  

def get_hists(froots, ext_paramname, ext_chains, cols, **kwargs):
    """
    find the give columns in the parameter file names, read them from the chain
    files and create the 1D histogramsz
    Parameters
    ----------
    froots: list
        file name roots
    ext_parmn: string
        extension of the parameter file names
    ext_chain: string
        extension of the chain file names
    columns: list
        short name of the columns to read
    kwargs: dictionary
        +verbose: bool
        +x_lim: list of lists
            list of x limits. Each element is a list containing the short parameter
            name and the max and minimum limit
        +num_bins: int 
            number of bins used to compute the histogram

    Output
    ------
    hists: dictionary
        key: the parameter short name; value: a list of numpy arrays with the
        normalised histograms in the same order as the froots
    bin_edges: dictionary
        key: the parameter short name; value: bin edges as returned by
        `numpy.histogram`
    x_labels: dictionary
        key: the parameter short name; value: x labels from the parameter files.
    """
    if kwargs.get('verbose', False):
        print('Reading the paramname files')
    # read the paramname files
    paramindex = cf.get_paramnames(froots, params=cols, ext=ext_paramname,
        verbose=kwargs.get('verbose', False))

    # dictionary containing the x limits. The short column name is the key.
    ranges = {i:None for i in cols}
    bin_edges = {i:None for i in cols}
    nbins = kwargs.get('num_bins', 10)
    for k, v1, v2 in kwargs.get('x_lim', []):
        if k not in ranges:
            warn("Key '{}' is not in column list".format(k), MyWarning)
        else:
            ranges[k] = [float(v1), float(v2)]
            bin_edges[k] = np.linspace(*ranges[k], num=nbins+1)
    # dictionary containing the x labels
    x_labels = {pi[1]: pi[2] for pi in paramindex[0]}
    # dictionary containing the histograms. The short column name is the key.
    hists = {i:[] for i in cols}

    # read the columns.
    if kwargs.get('verbose', False):
        print('Reading the chains')
    for fr, pindex in zip(froots, paramindex):
        # get the column numbers
        usecols = [i[0]+2 for i in pindex]
        chains = np.loadtxt(fr+ext_chains, usecols=[0]+usecols)
        # create the histograms
        for i, pi in enumerate(pindex):
                # compute the 1D histogram
                hist, bin_edge = np.histogram(chains[:,i+1], bins=nbins,
                        range=ranges[pi[1]], weights=chains[:,0])
                # normalise to the maximum and save in hists
                hists[pi[1]].append(hist/float(hist.max()))
                if ranges[pi[1]] is None:
                    ranges[pi[1]] = bin_edge[[0,-1]]
                    bin_edges[pi[1]] = bin_edge
 
    return hists, bin_edges, x_labels 

def make_figure(key_list, xlabels, n_cols):
    """
    Create the figure. Save the subplots into a
    dictionary and make the x axis labels
    Parameters
    ----------
    xlabels: dictionary
        list of keys of the following dictionary: short names from the
        parameter files
    n_cols: int
        number of columns of subplots, the number or rows is computed
        internally
    output
    ------
    fig: matplotlib figure
    axs_dic: dic
        dictionary of subplot objects with keys from key_list
    """
    # compute the number of rows needed
    if len(key_list)%n_cols == 0:
        n_rows = len(key_list)//n_cols
    else:
        n_rows = len(key_list)//n_cols + 1
        if len(key_list)%n_cols < n_cols/2:
            warn("There are too many empty subplots", MyWarning)

    fig, axs = plt.subplots(nrows=n_rows, ncols=n_cols, sharey=True, squeeze=False)
    axs = axs.flatten()

    # move the axs to a dictionary with the elements of key_list as key.
    # Set the x labels and remove the y labels except in the last plot
    axs_dic={}
    for i, k, ax in zip(it.count(), key_list, axs):
        if ax.is_first_col():
            ax.set_ylabel('$\mathcal{L}/\mathcal{L}_{max}$')
        ax.set_xlabel("${}$".format(xlabels[k]))
        axs_dic[k] = ax

    # save empty subplots
    if i+1 < axs.size:
        axs_dic['empty'] = axs[i+1:]

    return fig, axs_dic

def plot_hists(axs_dic, hists, bin_edges, labels):
    """
    Plot the histograms into the axes
    Parameters
    ----------
    axs_dic: dictionary
        key: parameter short name; value: axes/subplot object
    hists: dictionary
        key: the parameter short name; value: a list of numpy arrays with the
        normalised histograms in the same order as the froots
    bin_edges: dictionary
        key: the parameter short name; value: bin edges as returned by
        `numpy.histogram`
    labels: list
        legend labels. If None a cycle of None created
    """
    if not labels: # if labels is an empty list or None
        labels = it.cycle([None])

    for k, ax in axs_dic.items():
        if k == 'empty': 
           continue
        be = bin_edges[k]
        be = np.r_[be[0], np.repeat(be[1:-1], 2), be[-1]]
        ls = it.cycle(mpm.linestyles)
        col = it.cycle(mpm.colors)
        for h,l in zip(hists[k], labels):
            ax.plot(be, np.repeat(h, 2), ls=next(ls), c=next(col), label=l)

#############################################
###                 Main                  ###
#############################################

def main(argv):
    args = parse(argv)

    # command line parameter checks
    if args.legend and len(args.legend)!=len(args.ifroot):
        raise MyException("""The number of legend tags must be the same as the
                             number of file roots""")

    mpm.set_rcParams(**vars(args))

    if args.verbose:
        print("Reading the input files")
    # read the input files and create the histograms
    hists, bin_edges, x_labels = get_hists(args.ifroot, args.ext_parmn, args.ext_chain,
            args.columns, **vars(args))

    if args.verbose:
        print("Make the figure and axes")
    mpm.subsitute_labels(args.x_label, x_labels)
    # create the figure and axes
    fig, daxs = make_figure(args.columns, x_labels, args.n_cols)

    # set the x and y ranges
    xlimits = [[k, be[0], be[-1]] for k, be in bin_edges.items()]
    ylimits = [[k, 0, 1.05] for k, _ in bin_edges.items()]
    mpm.change_xylim(daxs, x_lim=xlimits, y_lim=ylimits)

    # plot the histograms
    if args.verbose:
        print("Plotting the histograms")
    plot_hists(daxs, hists, bin_edges, args.legend)

    # do the legend
    if args.legend:
        if args.verbose:
            print("Drawing the legend")
        if args.legend_plot is None:
            args.legend_plot = list(daxs.keys())[0]
        legend_ = mpm.draw_legend(fig, daxs, args.legend_plot, args.loc)

    # do vertical lines
    mpm.plot_horiz_vert(daxs, vert=args.vertical)

    # write the text
    dtexts = mpm.write_text(daxs, args.text, args.text_bkgr)

    # clear empty axes 
    for ax in daxs.get('empty', []):
        ax.set_axis_off()

    # show or save the figure
    plt.tight_layout(h_pad=0, pad=0.2, rect=args.bounds)
    if args.output is None:
        plt.show()
    else:
        fig.savefig(args.output)

if __name__ == "__main__":   # if is the main

    import sys
    main(sys.argv[1:])
    exit()
