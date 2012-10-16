#!/usr/bin/python
"""Container of misc functions"""

import numpy as np
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
