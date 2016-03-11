from __future__ import print_function
from __future__ import division
import six.moves as smo
import numpy as np
import re
import os
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyaC import pyaPermanent as PP
import scipy.interpolate as sci


class TransmissionCurves:
  """
    Photometric transmission curves for various photometric bands.
    
    A list of supported bands can be obtained by calling `availableBands`.
    
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
        elif len(l.rstrip("\n")) == 0:
          # Ignore empty lines
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
  
  def convolveWith(self, wvl, spec, bn, ik="linear"):
    """
      Convolve spectrum with transmission curve.
      
      Parameters
      ----------
      wvl, spec : arrays
          Wavelength axis [A] and spectral data.
      bn : string
          Name of the band.
      ik : string, optional
          The type of interpolation. Accepts all values also
          accepted by the `kind` keyword of scipy's `interp1d` rountine.
          Default is 'linear'.
      
      Returns
      -------
      Convolved spectrum : array
          Input spectrum multiplied with transmission curve
          of the specified band.
    """
    self._checkBand(bn)
    tc = self.getTransCurve(bn, ik)
    return spec*tc(wvl)
  
  def addPassband(self, name, wvl, trans, snc=False):
    """
      Add a new passband to the inventory.
      
      Parameters
      ----------
      name : string
          The name of the passband.
      wvl : array
          Wavelength in A.
      trans : array
          Transmission of the passband.
      snc : boolean, optional
          A `Skip Name Check` flag. If False (default), an
          exception is raised if the passband name is already
          present in the inventory of passbands. Otherwise
          the old passband is replaced by the new specification.
    """
    if not snc:
      if name in self.bands:
        raise(PE.PyANameClash("A passband with name '" + str(name) + "' is already present.", \
                              solution=["Change the name.", "Use `snc=True` to ignore and overwrite old passband."]))
    self.bands[name] = np.vstack((wvl, trans)).transpose()

  
  def addSpitzerIRACPassbands(self, forceDownload=False, verbose=True):
    """
      Adds Spitzer IRAC passbands.
      
      On first call, the passband files are downloaded.
      The files are downloaded from:
      
      http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/spectralresponse/
      
      Parameters
      ----------
      forceDownload : boolean, optional
          If True, a re-download of the passband files is triggered.
          Default is False.
      verbose : boolean, optional
          If True (default), download process will print information
          on progress.
    """
    
    fns = ["080924ch1trans_full.txt", "080924ch2trans_full.txt", "080924ch3trans_full.txt",
           "080924ch4trans_full.txt", "080924ch1trans_sub.txt", "080924ch2trans_sub.txt",
           "080924ch3trans_sub.txt", "080924ch4trans_sub.txt"]
    
    path = "pyasl/resBased/"
    
    for fn in fns:
      fno = path + fn
      if (not self._fs.fileExists(fn)) or forceDownload:
        self._fs.downloadToFile("http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/spectralresponse/" + fn, fno, forceDownload, verbose)
      
      dat = np.loadtxt(self._fs.requestFile(fno))
      dat[::,0] *= 1e4

      wl = re.match(".*IRAC\s+(\d+\.\d+)\s+.*", self._fs.requestFile(fno).readline()).group(1)
      if fn.find("_full") != -1:
        win = "_full"
      else:
        win = "_sub"
      
      bn = "IRAC" + wl + win
      self.addPassband(bn, dat[::,0], dat[::,1], snc=True)
    
  def __init__(self, fn="default"):
    self._fs = PP.PyAFS()
    self._readData(os.path.join(os.path.dirname(__file__), "transCurves.dat"))
    