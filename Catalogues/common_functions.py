"""Collections of functions specific to the analysis of catalogues
used in more than one code"""

import numpy as np
import os
import subprocess as sp

def dec2rad(deg):
  """Convert degrees to radians

  Parameters
  ----------
  deg: Nd-arrays
    angle(s) in degrees 

  output
  ------
  rad: Nd-array
    angle(s) in radians
  """
  return np.array(deg)*np.pi/180.
def rad2dec(rad):
  """Convert degrees to radians

  Parameters
  ----------
  rad: Nd-arrays
    angle(s) in radians 

  output
  ------
  dec: Nd-array
    angle(s) in degrees
  """
  return np.array(rad)*180./np.pi
def dec22rad2(deg2):
  """Convert degrees^2 to radians^2

  Parameters
  ----------
  deg: Nd-arrays
    solid angle(s) in degrees^2 

  output
  ------
  arg: Nd-array
    solid angle(s) in radians^2
  """
  return np.array(deg2)*np.pi*np.pi/32400.

def radec2rad(ra, dec):
  """Convert ra and dec to angle in radians as expected by healpix

  Parameters
  ----------
  ra, dec: Nd-arrays
    latitude and longitude in degrees 

  output
  ------
  theta, phi: Nd-array
    longitude and co-latitude in radians
  """

  theta = dec2rad( 90.0-np.array(dec) )
  phi = dec2rad( np.array(ra) )
  return theta, phi

def arcsec2deg(arcsec):
  """Convert arcsec to degrees 

  Parameters
  ----------
  arcsec: Nd-arrays
    angle(s) in arcseconds 

  output
  ------
  deg: Nd-array
    angle(s) in degrees
  """
  return np.array(arcsec)/3600.

def pow2(string):
  """
  This function check if the input is integer and if it a power of 2
  
  Parameters
  ----------
  string: string
    string to parse

  output: int
    return the integer
  """
  try:
    string = int(string)
  except:
    msg = "%r is not an integer" % string
    raise ap.ArgumentTypeError(msg)
  if( ((string & (string - 1)) == 0) and string > 0 ):
    return string   #return the integer
  else:
    msg = "%r is a power of 2" % string
    raise ap.ArgumentTypeError(msg)

def read_polygons(polyfile, catalogue, tempfile):
  """ Run the Mangle command 'polyid' on a mangle polygon file and a catalogue 
  and returns the list of poligons to which the objects in catalogue belong
  Parameters
  ----------
  polyfile: string
    Mangle polygon file name
  catalogue: string
    catalogue file name containing ra and dec in the first two columns
  tempfile: string
    temporary file where to store the output of 'polyid'
    is deleted before exiting the function

  output
  ------
  polyid: array of int
    returns the poligon numbers to which the objects in 'catalogue' belong
  """

  cmd = ['polyid', '-q', polyfile, catalogue, tempfile]  #create the command
  exe = sp.Popen(cmd, stdout=sp.PIPE).wait()   #execute the command and wait for its end
  polyid = np.loadtxt(cmd[-1], usecols=(2,), skiprows=1, dtype=int)  #read the file with the polygon ID
  os.remove(cmd[-1])   #remove the temporary file
  return polyid  #return the polygons id

def read_poly2sect(fname, verbose=False):
  """Read the file containing the polygon to sector pointers and return its content

  Parameters
  ----------
  fname: string
    file name containing in the first column the sector number to which the polygon,
    indicated by the row index, belong and, in the second column, the number of tiles covering the sector
  verbose: bool (optional)
    verbose mode

  output
  ------
  poly2sect, n_tiles: 1D arrays
    arrays with the sector number for each polygon (indicated by the index) and 
    the number of tiles per sector
  """

  poly2sect = np.loadtxt(fname, dtype=int)  #read the file
  if(verbose == True):
    print("Polygon and the files read.")
  return poly2sect[:,0], poly2sect[:,1]  #return the sector id and the number of tiles

def set_cosmology(om, ok, wde, h0=None):
  """Initialise a cosmology object from 'cosmologydir'
  Parameters
  ----------
  om: float
    omega_matter
  ok: float
    omega_curvature
  wde: float
    dark energy equation of state parameter
  h0: float (optional)
    reduced hubble constant
  output
  ------
  cosmo: cosmologydir.setcosmology.Setcosmology object
  """

  import cosmologydir as c
  if(h0 == None):
    cosmo = c.setcosmology.Setcosmology(om=om, ok=ok, wde=wde)
  else:
    cosmo = c.setcosmology.Setcosmology(om=om, ok=ok, wde=wde, h=h0)

  return(cosmo)
