from __future__ import print_function
from __future__ import division
import six.moves as smo
import numpy as np
import re
import os
from PyAstronomy.pyaC import pyaErrors as PE
import scipy.interpolate as sci


class TransmissionCurves:
  """
    Photometric transmission curves for various bands.
    
    Parameters
    ----------
    fn : string, optional
        The filename of the data file. Default is 'default',
        which means that the shipped data file will be used.
  """
  
  def _createSet(self, bn, lines):
    """
      Create data set from band name and set of lines.
      
      Parameters
      ----------
      bn : string
          Name of the band
      lines : list of string
          The data (table) specifying transmission as
          a function of wavelength.
    """
    if bn in self.bands:
      PE.warn(PE.PyAValError("Band named '" + str(bn) + "' already exists.", \
                             where="TransmissionCurves (reading data from file)", \
                             solution="Use unique band names."))
    self.bands[bn] = np.zeros( (len(lines), 2) )
    for i in smo.range(len(lines)):
      self.bands[bn][i,::] = np.array(lines[i].split(), dtype=np.float)
  
  def _readData(self, fn):
    """
      Read data file.
      
      The file must introduce each band via a row like
      "# BAND: 'band-name'". Otherwise, rows preceded by
      a hash will be ignored. The data must be specified as
      a table with columns wavelength [A] and transmission (0-1)
      below the line introducing the band; an arbitrary number
      of comment lines is allowed at all positions.
      
      Parameters
      ----------
      fn : string
          The filename of the data file.
    """
    if not os.path.isfile(fn):
      raise(PE.PyAValError("No such file: " + fn, \
                           where="TransmissionCurves", \
                           solution="Use valid file or use 'default'."))
    self.bands = {}
    cb = None
    lcol = None
    for l in open(fn):
      r = re.match("^\s*#\s*BAND:\s+'([^']+)'.*", l)
      if not r is None:
        # New current band
        # Add the previous data set
        if not cb is None:
          self._createSet(cb, lcol)
        # Prepare collecting data for new band
        cb = r.group(1)
        lcol = []
      else:
        if l.find('#') != -1:
          # It is a comment not specifying a band
          # Ignore this...
          continue
        else:
          # It is not a comment but data
          # Collect data
          lcol.append(l.rstrip("\n"))
    if len(lcol) > 0:
      # Add last data set
      self._createSet(cb, lcol)
  
  def _checkBand(self, bn):
    """
      Check whether band name is defined.
    """
    if not bn in self.bands:
      raise(PE.PyAValError("No such band: " + str(bn), \
                           where="TransmissionCurves", \
                           solution="Use one of: " + ', '.join(self.availableBands()) ))
  
  def availableBands(self):
    """
      Band names of available transmission curves.
      
      Returns
      -------
      Availabale names : list of strings
          All bands for which data are available.
    """
    return sorted(list(self.bands))
  
  def getTransCurve(self, bn, ik='linear'):
    """
      Get a tranmission curve.
      
      Parameters
      ----------
      bn : string
          Name of the band.
      ik : string, optional
          The type of interpolation. Accepts all values also
          accepted by the `kind` keyword of scipy's `interp1d` rountine.
          Default is 'linear'.
      
      Returns
      -------
      Tranmission curve : callable
          An object (scipy.interpolate.interp1d) that can be called
          with wavelength (float or array in [A]) as argument and returns
          the transmission. 
    """
    self._checkBand(bn)
    fi = sci.interp1d(self.bands[bn][::,0], self.bands[bn][::,1], kind=ik, bounds_error=False, fill_value=0.0)
    return fi
  
  def getTransCurveData(self, bn):
    """
      Get data specifying the transmission curve.
      
      Returns
      -------
      Transmission table : 2d array
          A table (array) with wavelength [A] in first column and
          tranmission (0-1) in the second column. 
    """
    self._checkBand(bn)
    return self.bands[bn]
    
  def __init__(self, fn="default"):
    self._readData(os.path.join(os.path.dirname(__file__), "transCurves.dat"))
    