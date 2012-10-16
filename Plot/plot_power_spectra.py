#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Plot rows of the window matrix"""

import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import sys


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

  p.add_argument( "" )

  p = apc.version_verbose( p, '1' )

  return p.parse_args(args=argv)

if __name__ == "__main__":   # if is the main

  args = parse(sys.argv[1:])



  exit()

