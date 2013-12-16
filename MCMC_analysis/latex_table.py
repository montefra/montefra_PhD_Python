#!/usr/bin/python
# -*- coding: utf-8 -*-

import contour_func as cf   #functions used in contour plot
import itertools as it
import my_statistic as ms
import numpy as np    # import numpy
import optparse as op  #import optsparse: allows nice command line option handling
import os    #contain OS dependent stuffs: it helps with sistem portability  
import sys    #mondule sys
import unify_chain as uc

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

  p.set_usage("""
  %prog [options] file_root1 tag1 [file_root2 tag2 ... file_rootn tagn]
  Given a list of file roots and of corresponding tags, reads the parameter name files and the chain files.
  Then computes the mean and either the standard deviation (default) or the given percentile
  of the parameters for each chain (chains assumed to have the following structure: chain[:,0]=weight;
  chain[:,1]=likelihood]; chain[:,2:]=parameters).
  The mean and the stddev or percentile are printed in a latex table with the following format:
  --------------------------------------------------------------------------------------
  		tag1			tag2			...	tagn
  --------------------------------------------------------------------------------------
  parameter1	m1+-s1(p1)(file1)	m1+-s1(p1)(file2)	...	m1+-s1(p1)(filen)
  parameter2	m2+-s2(p2)(file1)	m2+-s2(p2)(file2)	...	m2+-s2(p2)(filen)
  ...		...			...			...	...
  parameterm	mm+-sm(pm)(file1)	mm+-sm(pm)(file2)	...	mm+-sm(p1)(filen)
  --------------------------------------------------------------------------------------
  """)

  p.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Produces more output.")  # verbose option

  p.add_option("-w", "--what", action="store", dest="what", nargs=2, default=["s", "1"], help="Decide what to do: what[0]= 's' (standard deviation) or 'p' (percentile); if 's' the what[1]*sigma values saved, if 'p' the what[1] percentile values around the mean comuted and saved. [Default: %default]")
  p.add_option("-f", "--format", action="store", dest="fmt", default="%7.6f", help="Formatter for the numbers in the output table. [Default: %default]")

  p.add_option("-r", "--rescale", action="append", dest="rescale", type=float, nargs=2, help="Rescale the values in position rescale[0] by multiplying the values by rescale[1]. The same rescaling factor is shown in the parameter name.")

  p.add_option("-o", "--output", action="store", dest="outfile", default="latextable.tex", help="Output file name. [Default: %default]")

  p.add_option("--comment", action="store_true", dest="comment", default=False, help="Comment rows of the matrix when all the elements have null variance")

  p.add_option("-c", "--chain-ext", action="store", dest="cext", default=".0.5_total.txt", help="'file_rooti+cext': file containing the chains. [Default: %default]")
  p.add_option("-p", "--parname-ext", action="store", dest="pext", default=".paramnames", help="'file_rooti+cext': file containing the name of the parameters. [Default: %default]")

  return p.parse_args()

def check_optarg((opt, args)):
  """
  Check if the input options and parameters are ok.
  
  Parameters
  ----------
  options: dictionary
    dictionary of options returned by optparse
  args: list
    list of arguments returned by optparse

  output:
    opt: dictionary of the options
    file_roots: list of file roots
    tags: list of tags
  """

  if(len(args) < 2 or len(args)%2!=0):
    raise SystemExit("The program requires file name roots and followed by the associated tags: file_root1 tag1 [file_root2 tag2 ... file_rootn tagn]") 
  else:
    file_roots = args[::2]
    tags = args[1::2]

  #check 'opt.what'
  if not (opt.what[0] in ("s", "p")):
    raise KeyError("The option '-w','--what' accepts only 's' and 'p' as firts elements")

  return opt, file_roots, tags

def print_table(fname, mean, std_perc, tags, paramnames, std=True, fmt="%7.6f", rescale=None, comment=False):
  """
  Print on 'fname' the mean and standard deviation or percentile 
  for the parameters in paramnames and uses tags as the first row ofthe matrix

  Parameters
  ----------
  fname: string
    output file name
  mean: 1D array
    mean of the parameters
  std_perc: 1D or 2D array
    standard deviation of the parameters and scores of the percentiles. 
    The first column are the percentiles to consider. For each percentile, 
    the two scores at (percentile/2.) and (100 - percentile/2.) must be listed
  tags: list
    list of tags to be used in the first row of the matrix
  paramnames: list
    list of strings containing the parameter names with the format "name \t latex code"
  std: bool (optional)
    if 'True', 'std_perc' expected as 1D array containing the standard deviation;
    otherwise 'std_perc' expected as 2D array containing the percentile.
  fmt: string (optional)
    formatting of the output
  rescale: int or list (optional)
    multply the mean and standard deviation by rescale. 
    If None no rescaling apply, if a single float given, it is used for all the values, 
    if len(rescale)=len(paramnames) to each parameter is multiplied according to the corresponding value.
    If rescale[i] is different from one, the same factor is written in front of the parameter name
  comment: bool (optional)
    comment rows where all the columns are '--'
  """
  #check rescale
  if(rescale==None):
    rescale = [1 for i in range(len(paramnames))] #initialise the list for the rescaling of the given parameters
  elif(type(rescale).__name__ == 'float' or type(rescale).__name__ == 'int'):
    rescale = [float(rescale) for i in range(len(paramnames))] #initialise the list for the rescaling of the given parameters
  elif(type(rescale).__name__ == 'list'):
    if(len(rescale) != len(paramnames)):
      raise SystemExit("If 'rescale' is a list it must be the same length of paramnames")
  else:
    raise TypeError("Accepted types for 'rescale': None, float, integer, list.")


  out = open(fname, 'w')    #open file
  out.write("\\documentclass[10pt]{article}\n \\begin{document}\n\n")
  out.write("\\begin{table*}\n  \\centering\n  \\begin{minipage}{160mm}\n")   #LaTex table
  out.write("    \\caption{1D marginalised constraints for the cases: %s} \\label{tab:1Dmarg}\n" % (", ".join(tags)))  #caption and label
  out.write("    \\begin{tabular}{ l%s }\n    \\hline\n" %(" c"*len(tags) ))   # table type
  out.write("      & %s \\\\\n      \\hline\n" %(" & ".join(tags)))   #write the header 
  for i,p,r in it.izip(it.count(), paramnames[0], rescale):   #for 1: rows of the table
    outrow='      $'  #row of the matrix
    if(np.absolute(r-1) > 1e-3):   #if a rescale required
      if(r.is_integer()==True):
	strr = str(int(r))
      else:
	strr = str(r)
      outrow += '{0}'.format(strr)
    outrow += '{0[2]}$'.format(p)
    skipped = 0  #skipped columns
    for m,sp in it.izip(mean,std_perc):  #for 2: mean +- std (or percentile)
      if(std==True):
	if(np.absolute(sp[i]/m[i]) > 1e-07):
	  outrow += str(" & $"+fmt+"\pm"+fmt+"$") %(r*m[i], r*sp[i])   #add param mean +- std
	else:
	  outrow += " & -- "    #skip constant parameter
	  skipped += 1   #has been skipped
      else:
	if(np.absolute( (sp[0,i]+sp[1,i])/(2.*m[i]) ) > 1e-07):
	  outrow += str(" & $"+fmt+"_{"+fmt+"}^{+"+fmt+"}$") %((r*m[i]), (r*sp[0,i]), (r*sp[1,i]))   #write mean +- 68perc
	else:
	  outrow += " & -- "
	  skipped += 1   #has been skipped
    #end for 2: mean +- std (or percentile)
    if( comment==True and skipped == len(tags) ):  #if the empty lines are to be commented
      outrow = "%"+outrow
    out.write(outrow+"\\\\[1.5mm]\n")
  #end for 1: rows of the table
  out.write("      \\hline\n    \\end{tabular}\n  \\end{minipage}\n\\end{table*}")   #LaTex table
  out.write("\n\n\\end{document}\n")
  out.close()  #close file

  pass

if __name__ == "__main__":   # if is the main

  """read a list of chains and the corresponding parameter name files and 
  creates a latex table"""

  opt, froots, tags = check_optarg(options(op.OptionParser(version="%prog version 1")))   #create the object optparse

  paramnames = cf.get_paramnames(froots, ext=opt.pext, verbose=opt.verbose)

  rescale = [1 for i in range(len(paramnames))] #initialise the list for the rescaling of the given parameters
  if(opt.rescale != None):   #if a tuple of values to draw the horizontal line is given order it
    for h in opt.rescale:
      rescale[int(h[0])] = h[1]

  mean, std_perc = [],[]  #initialise the lists that will contain the mean and stddev or percentiles
  if(opt.what[0] == "s"):  #check if standard deviation of percentile required
    std, perc = True, None
  else:
    std, perc = False, float(opt.what[1])
  for fn in froots:   #read the chains and computed what is required
    if(opt.verbose == True):
      print "Computed mean and standard deviation or percentile on file '%s'" %(fn+opt.cext)
    tm, tsp = uc.marg(np.loadtxt(fn+opt.cext), sd=std, percentile=perc, verbose=opt.verbose)
    if(std == False):
      tsp = tsp[:,1:]    #do not consider the first column that contains the values of the percentile
      tsp = tsp - np.vstack([tm,]*2)
    mean.append(tm)
    std_perc.append(tsp)

  print_table(opt.outfile, mean, std_perc, tags, paramnames, std=std, fmt=opt.fmt, rescale=rescale, comment=opt.comment)

  sys.exit(0)

