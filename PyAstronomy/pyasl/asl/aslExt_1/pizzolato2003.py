import numpy as np
import os
from PyAstronomy.pyaC import pyaErrors as PE

class Pizzolato2003:
  
  def _readData(self, fn):
    """
      Read data from file.
    """
    data = np.loadtxt(fn)
    self._tab3 = data[0:7,::].copy()
    self._tab4 = data[7:,::].copy()
  
  def _calcVal(self, tab, row, isat, iprsat, pr):
    """
      Evaluates log10(Lx) or log10(Lx/Lbol) and error
      
      Parameters
      ----------
      tab : array
          The data table, i.e., Table 3 or 4 from
          Pizzolato et al. 2003.
      isat : int
          Index of value representing saturated level.
      iprsat : int
          Index of limiting rotation period for reaching
          saturation.
      
      Returns
      -------
      log val : float
          The logarithm (base 10) of Lx or Lx/Lbol
      error log val : float
          The error on that value derived using error
          propagation.
    """
    if pr <= tab[row, iprsat]:
      # Saturated regime
      return tab[row, isat], tab[row, isat+1]
    lnx = tab[row, isat] + 2.0*np.log10(tab[row, iprsat]) - 2.0*np.log10(pr)
    errlnx = np.sqrt( tab[row, isat+1]**2 + \
                      (2.0/(tab[row, iprsat]*np.log(10.0)))**2*tab[row, iprsat+1]**2  )
    return lnx, errlnx
  
  def _findRow(self, tab, val, valname):
    """
      Find corrected row in table.
      
      Parameters
      ----------
      tab : array
          The data table, i.e., Table 3 or 4 from
          Pizzolato et al. 2003.
      val : float
          Value to be searched for.
      valname : string
          The name of the parameter. Needed in
          case of error.
      
      Returns
      -------
      row : int
          The index of the correct row.
    """
    row = np.where((val >= tab[::,0]) & (val < tab[::,1]))[0]
    if len(row) == 0:
      conc = np.concatenate((tab[::,0], tab[::,1]))
      raise(PE.PyAValError(valname + " value of " + str(val) + " is out of range. " +
                           "Allowed range: " + str(np.min(conc)) + " - " + str(np.max(conc)))
                           )
    return row[0]
  
  def _checkPr(self, pr):
    """
      Check whether rotation period is valid.
    """
    if pr <= 0.0:
      raise(PE.PyAValError("Rotation period must be >= 0. Given value: " + str(pr)))
  
  def log10lxbv(self, bv, pr):
    """
      Estimate log10(Lx)
      
      Parameters
      ----------
      bv : float
          B-V color [mag]
      pr : float
          Rotation period [d]
      
      Returns
      -------
      log10(Lx) : float
          X-ray luminosity estimate
      Error log10(Lx) : float
          Uncertainty
    """
    self._checkPr(pr)
    row = self._findRow(self._tab4, bv, "B-V")
    return self._calcVal(self._tab4, row, 3, 7, pr)
  
  def log10lxlbolbv(self, bv, pr):
    """
      Estimate log10(Lx/Lbol)
      
      Parameters
      ----------
      bv : float
          B-V color [mag]
      pr : float
          Rotation period [d]
      
      Returns
      -------
      log10(Lx/Lbol) : float
          X-ray luminosity estimate
      Error log10(Lx/Lbol) : float
          Uncertainty
    """
    self._checkPr(pr)
    row = self._findRow(self._tab4, bv, "B-V")
    return self._calcVal(self._tab4, row, 5, 9, pr)
  
  def log10lxmass(self, m, pr):
    """
      Estimate log10(Lx)
      
      Parameters
      ----------
      m : float
          Stellar mass [MSun]
      pr : float
          Rotation period [d]
      
      Returns
      -------
      log10(Lx) : float
          X-ray luminosity estimate
      Error log10(Lx) : float
          Uncertainty
    """
    self._checkPr(pr)
    row = self._findRow(self._tab3, m, "Mass")
    return self._calcVal(self._tab3, row, 3, 7, pr)
  
  def log10lxlbolmass(self, m, pr):
    """
      Estimate log10(Lx/Lbol)
      
      Parameters
      ----------
      m : float
          Stellar mass [MSun]
      pr : float
          Rotation period [d]
      
      Returns
      -------
      log10(Lx/Lbol) : float
          X-ray luminosity estimate
      Error log10(Lx/Lbol) : float
          Uncertainty
    """
    self._checkPr(pr)
    row = self._findRow(self._tab3, m, "Mass")
    return self._calcVal(self._tab3, row, 5, 9, pr)
  
  def __init__(self):
    """
      X-ray luminosity/rotation period relations by Pizzolato et al. 2003.
      
      Implementation of the relations between X-ray luminosity and
      rotation period given by Pizzolato et al. 2003, A&A 397, 147P.
      
      To find the X-ray luminosity, the following expression is evaluated
      
      .. math::
          x = \log_{10} (x_{sat}) + 2\log_{10}(P_{r,sat}) - 2\log_{10}(P_r)
      
      where x denotes either Lx or Lx/Lbol and 'sat' indicates the value
      at saturation. The coefficients can be found in Tables 3 and 4
      of Pizzolato et al. 2003.
    """
    self._readData(os.path.join(os.path.dirname(__file__), "pizzolato2003tables.dat"))
  
  