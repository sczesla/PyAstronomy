from __future__ import print_function, division
from PyAstronomy.pyaC import pyaErrors as PE
import gzip
import re
import os
import numpy as np
import six.moves as smo

def readUnit7(fn, minwl=None, maxwl=None):
  """
    Read PHOENIX Unit 7 files.
    
    The documentation of the PHOENIX output
    files may be found in the PHOENIX manual (Sect. PHOENIX files):
    `http://www.hs.uni-hamburg.de/EN/For/ThA/phoenix/manual.html`.
    
    Parameters
    ----------
    fn : string
        Filename. Name ending in .gz will be treated
        as gzipped files.
    minwl : float, optional
        Minimum wavelength to be considered [A].
    maxwl : float, optional
        Maximum wavelength to be considered [A].
    
    Returns
    -------
    Model spectrum : array
        Array with 4 column. 1) The wavelength [A],
        2) log10 of the flux [erg/s/cm^2/cm],
        3) log10 of the Planck flux [erg/s/cm^2/cm]
        function, 4) log10(flux(tau_std(2))) [erg/s/cm^2/cm].
  """
  if not os.path.isfile(fn):
    raise(PE.PyAValError("There is no file named: " + str(fn)))
  # Check whether it is a gzipped file
  r = re.match(".*\.gz", fn)
  if r is not None:
    f = gzip.open(fn, 'r')
  else:
    f = open(fn, 'r')
  lines = f.readlines()
  f.close()
  result = np.zeros((len(lines), 4))
  for i in smo.range(len(lines)):
    result[i,::] = np.array(lines[i].replace('D','E').split()[0:4], dtype=np.float)
  # Remove entries using min and max wavelength
  if minwl is not None:
    indi = np.where(result[::,0] >= minwl)[0]
    result = result[indi,::]
  if maxwl is not None:
    indi = np.where(result[::,0] <= maxwl)[0]
    result = result[indi,::]
  # Sort with respect to wavelength
  result = result[np.argsort(result[::,0]),::]
  return result


def readDTable(fn, sort=None):
  """
    Read a table with numbers in d-format.
    
    In d-format, a number, e.g., looks 1.3d-3.
    
    Parameters
    ----------
    fn : string
        Filename.
    sort : int, optional
        If given, the table will be sorted with respect
        to that column. The first columns is 0.
    
    Returns
    -------
    table : array
        An array with as many columns and rows as the table
        read from file.
  """
  if not os.path.isfile(fn):
    raise(PE.PyAValError("There is no file named: " + str(fn)))
  # Check whether it is a gzipped file
  r = re.match(".*\.gz", fn)
  if r is not None:
    f = gzip.open(fn, 'r')
  else:
    f = open(fn, 'r')
  lines = f.readlines()
  f.close()
  # Find number of columns
  ncol = len(lines[0].split())
  result = np.zeros((len(lines), ncol))
  for i in smo.range(len(lines)):
    result[i,::] = np.array(lines[i].replace('D','E').split(), dtype=np.float)
  if sort is not None:
    # Sort with respect to column
    result = result[np.argsort(result[::,sort]),::]
  return result


def decomposeFilename(fn):
  """
    Decompose PHOENIX filename.
    
    Parameters
    ----------
    fn : string
        The filename.
    
    Returns
    -------
    Parameters : dictionary
        A dictionary with the following keys:
         - teff: The effective temperature in K
         - logg: Log(g [cm/s**2])
         - met: Metallicity (M/H)
         - fn: Complete filename
         - notParsed: Part of the filename not parsed for teff, logg, and metallicity.
        Note that `None` is returned if the filename could not be
        parsed.
  """
  r = re.match("lte(\d+)([+-])(\d+\.\d+)([-+]\d+\.\d+)(.*)", fn)
  if r is None:
    PE.warn(PE.PyAValError("Cannot decompose PHOENIX filename: " + str(fn)))
    return None
  result = {}
  result["fn"] = fn
  teff = int(r.group(1))
  if teff < 1000:
    result["teff"] = teff * 100
  else:
    result["teff"] = teff
  result["logg"] = float(r.group(3))
  if r.group(2) == "+":
    result["logg"] *= -1.0
  result["met"] = float(r.group(4))
  result["notParsed"] = r.group(5)
  return result