#!/usr/bin/python
"""Container of misc functions"""

import numpy as np
import re
import sys

def convert2array(a):
  """check if 'a' is an iterable or not and convert it to a numpy array and flatten it.
  Useful to convert a variable to iterable quantity
  Parameters
  ----------
  a: 
    float, list, tuble, numpy array
  output: 
    numpy array with the same number of elements of a
  """
  try:
    iter(a)
  except TypeError:
    a = np.array([a])
  else:
    a = np.array(a)
  return a.flatten()

def check_file_exist(fname):
  """Check if file 'fname' already exists
  
  Parameters
  ----------
  fname: string
    file name to check

  output
  ------
  exists: bool
    true if it exists, false otherwise
  """
  try:
    fh = open(fname)
  except IOError as e:
    return False
  else:
    return True

def insert_src(src, insert, overwrite=False, skip=False):
  """Create a new string from 'src' inserting insert[0] before insert[1]

  Parameters
  ----------
  src: string
    string to modify
  insert: list
    insert insert[0] before insert[1] in scr
  overwrite: bool (optional)
    if false consider the string as an output file name and check if exist
  skip: bool (optional)
    if false the code is aborted

  output
  ------
  ostr: string
    output string
  toskip: bool
    true if 'skip' is true and ostr already exists as file
  """
  ostr = src.replace(insert[1], insert[0]+insert[1])
  if(overwrite == True):  
    toskip = False
  else:  #check if the output file already exists
    toskip = check_file_exist(ostr)   #check if the ostr exists as file
    if(toskip == True):
      print("The file '{}' already exists. ".format(ostr))
      if(skip != True):
	print("Aborting")
	sys.exit(10)
  return ostr, toskip

def replace_src(src, replace, overwrite=False, skip=False):
  """Create a new string from 'src' replacing replace[0] with replace[1]

  Parameters
  ----------
  src: string
    string to modify
  replace: list
    replace replace[0] with replace[1] in scr
  overwrite: bool (optional)
    if false consider the string as an output file name and check if exist
  skip: bool (optional)
    if false the code is aborted

  output
  ------
  ostr: string
    output string
  toskip: bool
    true if 'skip' is true and ostr already exists as file
  """
  ostr = src.replace(replace[0], replace[1])
  if(overwrite == True):  
    toskip = False
  else:  #check if the output file already exists
    toskip = check_file_exist(ostr)   #check if the ostr exists as file
    if(toskip == True):
      print("The file '{}' already exists. ".format(ostr))
      if(skip != True):
	print("Aborting")
	sys.exit(10)
  return ostr, toskip

def create_ofile_name(f, **kwargs):
    """
    Create a file name out of the name (of) 'f' substituting or inserting.
    Is uses "insert_src" and "replace_src"
    Parameters
    ----------
    f: file object or string
        file containing the catalogue
    output
    ------
    ofile: string
        output file name or *None*, if *skip* is *True* and ofile alread exists 

    accepted kwargs that affects the function
    +verbose: verbose mode [True|False] 
    +replace: replace string *replace[0]* with *replace[1]* in f.name
    +insert: insert string *insert[0]* before *insert[1]* in f.name
    +skip: existing file names skipped [True|False]
    +overwrite: existing file names overwritten [True|False]
    """

    try:
        fname = f.name
    except AttributeError:
        fname = f

    if(kwargs['verbose'] == True):
        print("Process catalogue '{0}'.".format(fname))

    #create the output file name and check it
    if(kwargs['replace'] == None):
        ofile, skip = insert_src(fname, kwargs['insert'],
            overwrite=kwargs['overwrite'], skip=kwargs['skip'])
    else:
        ofile, skip = replace_src(fname, kwargs['replace'],
            overwrite=kwargs['overwrite'], skip=kwargs['skip'])
    if(skip == True):
        print("Skipping file '{0}'".format(fname))
        return None
    return ofile

def n_lines_comments(f, comment='#'):
    """
    Count the number of lines starting with a comment at the beginning of the file.
    The function returns at the first non commented line
    Parameters
    ----------
    f: file object or string
        file to check
    comment: string
        comment character
    output
    ------
    n_lines: int
        number of lines starting with the comment
    """
    def search_commented_lines(fo, c):
        "fo: file object; c: character"
        n_lines = 0
        pattern = re.compile("^\s*{0}".format(c))
        for l in fo:
            if pattern.search(l) is None:
                break
            else:
                n_lines += 1
        return n_lines

    try:   # if it is a file name
        with open(f, 'r') as fo:  # open it and read up to the end of the comments
            return search_commented_lines(fo, comment)
    except TypeError: # if it is a file object
        initial_position = f.tell()  # get the current file position
        f.seek(0)  # go back to the beginning of the file
        n_lines = search_commented_lines(f, comment)
        f.seek(initial_position) # put back the file in the previous position
        return n_lines

