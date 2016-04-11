from __future__ import print_function, division
import numpy as np
import os
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo

class Ramirez2005:
  """
    Relation between effective temperature and color given by Ramirez and Melendez 2005.
    
    Ramirez and Melendez 2015, ApJ 626, 465-485 (please not that Ramirez has a non-ASCII
    accent on the i and Melendez an accent on the second e) give a relation between the
    stellar effective temperature and the color concerning various bands. This class
    allows to carry out the conversion in both directions. 
  """
  
  def _extractTableData(self, lines, tableno):
    """
      Extract lines pertaining to specified table.
      
      Parameters
      ----------
      lines : list of strings
          Content of the data file.
      tableno : int
          Number of the table to be extracted.
      
      Returns
      -------
      Table data : list of strings
          Part of the file belonging to the
          specified table.
    """
    result = []
    inTable = False
    tid = "Table " + str(tableno) + ":"
    for l in lines:
      if l.find(tid) == -1:
        if not inTable:
          continue
      else:
        inTable = True
        continue
      if inTable:
        if len(l) > 0:
          result.append(l)
          continue
      break
    return result
  
  def _convertBandName(self, bn):
    """
      Convert band name used in tables to internal representation.
      
      Parameters
      ----------
      bn : string
          Band name used in table.
      
      Returns
      -------
      Band ID : string
          Identifier used in the class. 
    """
    name = bn.replace("{", '').replace("}", '').replace(" ", '')
    if name.startswith("("):
      name = name[1:-1]
    return name
        
  def _readTab23(self, lines, tableno):
    """
      Read tables 2 and 3.
      
      Parameters
      ----------
      lines : list of strings
          Content of the data file.
      tableno : int
          Number of the table to be extracted.
      
      Returns
      -------
      bands : list of strings
          IDs of all bands in the table.
      result : array
          Table data as array.
    """
    dat = self._extractTableData(lines, tableno)
    result = np.zeros((len(dat), 8))
    bands = [None]*len(dat)
    for i, l in enumerate(dat):
      s = l.split()
      bands[i] = self._convertBandName(''.join(s[0:len(s)-8]))
      result[i,::] = np.array(s[-8:], dtype=np.float)
    return bands, result
  
  def _readTab45(self, lines, tableno):
    """
      Read tables 4 and 5.
      
      Parameters
      ----------
      lines : list of strings
          Content of the data file.
      tableno : int
          Number of the table to be extracted.
      
      Returns
      -------
      result : array
          Table data as array.
    """
    dat = self._extractTableData(lines, tableno)
    result = np.zeros((len(dat), 10))
    for i, l in enumerate(dat):
      # zero can be used instead of ldots
      l = l.replace("\ldots", "0.0")
      s = l.split()
      if len(s) > 10:
        # A new band identifier is given
        result[i,::] = np.array(s[-10:], dtype=np.float)
      else:
        result[i,::] = np.array(s[:], dtype=np.float)
    return result
    
  def _readData(self, fn):
    """
      Read the table data.
      
      Parameters
      ----------
      fn : string
          Filename.
    """
    # Read data and remove newlines
    lines = [s.rstrip("\n") for s in open(fn).readlines()]

    # Tables 2 and 3 have the same structure
    self._bands, self._tab2 = self._readTab23(lines, 2)
    _, self._tab3 = self._readTab23(lines, 3)
    # Tables 4 and 5 have the same structure
    self._tab4 = self._readTab45(lines, 4)
    self._tab5 = self._readTab45(lines, 5)
  
  def _checkBand(self, band):
    """
      Check whether band identifier is valid.
    """
    if not band in self._bands:
      raise(PE.PyAValError("No such band: " + str(band), \
                           solution="Use either of: " + ", ".join(self._bands)))
  
  def _checkST(self, st):
    """
      Check whether stellar type (main-sequence/giant) is valid.
    """
    if not st in ["ms", "g"]:
      raise(PE.PyAValError("No such stellar type: " + str(st), \
                           solution="Either use 'ms' for main sequence or 'g' for giant."))
  
  def _resolveMetallicityIndex(self, feh):
    """
      Determine where to find coefficients for given metallicity in Tables 4 and 5.
      
      Parameters
      ----------
      feh : float
          Metallicity
    """
    fbs = [[-0.5,0.5], [-1.5,-0.5], [-2.5,-1.5], [-4.0,-2.5]]
    for i, fb in enumerate(fbs):
      if (feh >= fb[0]) and (feh < fb[1]):
        return i
    raise(PE.PyAValError("No data for given metallicity of " + str(feh), \
                         solution="See RM05, Tables 4 and 5 for ranges of validity."))
    
  def colorToTeff_nop(self, band, X, feH, stype="ms"):
    """
      Converts color into effective temperature according to Eq. 1.
      
      The conversion using to Eq. 1 neglects a polynomial correction
      for metallicity. According to RM05, this causes a systematic
      error on the order of '30 or 40 K'.
      
      Parameters
      ----------
      band : string
          Band identifier.
      X : float
          The color index (e.g., value of B-V).
      feH : float
          Metallicity
      stype : string, {ms, g}
          Type of star (main sequence or giant).
      
      Returns
      -------
      Teff : float
          The effective temperature in K.
    """
    self._checkBand(band)
    self._checkST(stype)
    if stype == "ms":
      dat = self._tab2
    else:
      dat = self._tab3
    
    # Row of coefficients in the table.
    j = self._bands.index(band)
    
    # Calculate according to Eq. 1
    result = 0.0
    for i in smo.range(3):
      result += dat[j, i] * X**i
    result += dat[j, 3]*X*feH
    for i in smo.range(4,6):
      result += dat[j, i] * feH**i
    
    teff = 5040./result
    return teff
  
  def colorToTeff(self, band, X, feH, stype="ms", ignoreRange=False):
    """
      Converts color into effective temperature according to Eq. 2.
      
      This method takes the polynomial correction into account. Note
      that no interpolation is implemented between the polynomials
      defined in Tables 4 and 5, but the appropriate polynomial
      (according to footnote (a) on under the tables) is used.
      
      Parameters
      ----------
      band : string
          Band identifier.
      X : float
          The color index (e.g., value of B-V).
      feH : float
          Metallicity
      stype : string, {ms, g}
          Type of star (main sequence or giant).
      ignoreRange : boolean, optional
          If True, the validity range of the relations will
          be ignored. Otherwise (default) an exception will
          be raised when a value outside the range is
          encountered.
      
      Returns
      -------
      Teff : float
          The effective temperature in K.
    """
    self._checkBand(band)
    self._checkST(stype)
    if stype == "ms":
      dat = self._tab2
      datp = self._tab4
    else:
      dat = self._tab3
      datp = self._tab5
    
    # Row of coefficients in the table.
    j = self._bands.index(band)
    
    # Index where to find coefficients for polynomial
    # correction in tables 4 and 5. The factor 4 occurs,
    # because there are 4 distinct metallicity bands.
    mj = j*4 + self._resolveMetallicityIndex(feH)
    
    # Check range of validity
    if (X < datp[mj,1] or X > datp[mj,2]) and (not ignoreRange):
      raise(PE.PyAValError("The given color index of " + str(X) + " is outside the validity range of " + str(datp[mj,1]) + " - " + str(datp[mj,2]), \
                           where="colorToTeff", \
                           solution="Use index within valid range."))
    
    # Calculate according to Eq. 1
    result = 0.0
    for i in smo.range(3):
      result += dat[j, i] * X**i
    result += dat[j, 3]*X*feH
    for i in smo.range(4,6):
      result += dat[j, i] * feH**(i-3)
    
    teff = 5040./result
    
    # Polynomial correction
    for i in smo.range(7):
      teff += (datp[mj, 3+i] * X**i)
    
    return teff
  
  def teffToColor(self, band, teff, feH, stype="ms", dteff=0.01, maxiter=100):
    """
      Converts effective temperature into color according to Eq. 2.
      
      This method inverts Eq. 2 using an iterative scheme.
      
      Parameters
      ----------
      band : string
          Band identifier.
      teff : float
          Effective temperature in K.
      feH : float
          Metallicity
      stype : string, {ms, g}, optional
          Type of star (main sequence or giant).
      dteff : float, optional
          Temperature difference to be reached
          by the iteration [K]. Default is 0.01.
      maxiter : int, optional
          The maximum number of iterations to be
          carried out. Default is 100.
          
      Returns
      -------
      X : float
          Color in the specified band.
    """
    self._checkBand(band)
    self._checkST(stype)
    
    # Trial step-width in x (=color) to find iterative solution
    dx = 0.1
    
    # Obtain guess solution neglecting polynomial
    x0 = self.teffToColor_nop(band, teff, feH, stype)
    # Iteratively improve solution
    counter = 0
    while True:
      counter += 1
      tc0 = self.colorToTeff(band, x0, feH, stype, ignoreRange=True)
      tc1 = self.colorToTeff(band, x0+dx, feH, stype, ignoreRange=True)
      dtdx = (tc1-tc0) / dx
      dt = tc0 - teff
    
      if abs(dt) < dteff:
        break
      else:
        x0 -= dt/dtdx
    
      if counter > maxiter:
        raise(PE.PyAValError("Maxiter reached. Could not convert teff into color.", \
                             where="teffToColor", \
                             solutions="Try to lower precision via adjusting 'dteff'."))
      
    return x0

  def teffToColor_nop(self, band, teff, feH, stype="ms", noRaise=False):
    """
      Converts effective temperature into color according to Eq. 1.
      
      This method inverts Eq. 1. Note that the equation is parabolic
      in the color (i.e., X). Therefore, there are two solutions of
      which the one falling within the validity ranges specified in
      Tables 4 and 5 of RM05 is selected. If none or both of the
      solutions are valid, an exception is raised.  
      
      Parameters
      ----------
      band : string
          Band identifier.
      teff : float
          Effective temperature in K.
      feH : float
          Metallicity
      stype : string, {ms, g}
          Type of star (main sequence or giant).
      noRaise : boolean, optional
          If True, no exceptions will be raised, but warnings
          will be given Both candidate solutions will be
          returned in this case.
      
      Returns
      -------
      X : float
          Color in the specified band.
    """
    self._checkBand(band)
    self._checkST(stype)
    if stype == "ms":
      dat = self._tab2
      datp = self._tab4
    else:
      dat = self._tab3
      datp = self._tab5
    
    # Use warnings of 'noRaise' is True
    if noRaise:
      raiser = PE.warn
    else:
      def throw(x):
        raise(x)
      raiser = throw
    
    # Row of coefficients in the table.
    j = self._bands.index(band)

    c = dat[j,0] + dat[j, 4]*feH + dat[j,5]*feH**2 - 5040./teff
    a13p = dat[j,1] + dat[j,3]*feH
    a2 = dat[j,2]
    
    sq = np.sqrt(a13p**2 - 4.*a2*c)
    X = [(-a13p + sq)/(2.*a2), (-a13p - sq)/(2.*a2)]
    
    # Get validity ranges to decide, which solution is correct
    mj = j*4 + self._resolveMetallicityIndex(feH)
    
    xmin = datp[mj,1]
    xmax = datp[mj,2]
    xin = [False]*2
    for i in smo.range(2):
      if X[i] >= xmin and X[i] <= xmax:
        xin[i] = True
    if sum(xin) == 2:
      raiser(PE.PyAAlgorithmFailure("No unique solution in inversion of teff-color relation.", \
                                   where="teffToColor_nop",\
                                   solution="Consider carefully(!) using 'noRaise'."))
      # Only relevant if 'noRaise' is True
      return tuple(X)
    if sum(xin) == 0:
      raiser(PE.PyAAlgorithmFailure("No valid solution in inversion of teff-color relation.\n" + \
                                   "Band: " + str(band) + ", Range of validity: %.4e - %.4e" % (xmin, xmax) + \
                                   ", candidate solutions: %.4e and %.4e" % tuple(X) + ", Input Teff: %5.1f" % teff, \
                                   where="teffToColor_nop",\
                                   solution="Consider carefully(!) using 'noRaise'."))
      # Only relevant if 'noRaise' is True
      return tuple(X)
    
    if xin[0]:
      return X[0]
    else:
      return X[1]
  
  def availableBands(self):
    """
      Get a list of available band identifiers.
      
      Returns
      -------
      Band IDs : list of strings
          All strings used to identify bands. 
    """
    return self._bands[:]
  
  def __init__(self):
    self._readData(os.path.join(os.path.dirname(__file__), "ramirez2005tables.dat"))
