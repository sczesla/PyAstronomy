from __future__ import print_function, division
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo
import re


class SpecTypeDeJager(object):
  """
    Spectral type calibration from de Jager and Nieuwenhuijzen.
    
    This class implements the spectral type calibration presented by
    de Jager and Nieuwenhuijzen 1987, A&A 177, 217-227 (DJ87). Specifically,
    the authors calibrate the relation between spectral type and stellar
    luminosity and effective temperature for a wide range of spectral types
    and luminosity classes.
    
    DJ87 give their results in the form of a polynomial relation with 20
    coefficients (Eqs. 2a and 2b) as well as in the form of their
    Tables 5 and 6.
    
    .. Note:: Tables 5 and 6 in DJ87 are calculated on the basis of a
              40 parameter solution, which is not given in the paper.
              The numbers based on the 20 parameter polynomial solution
              deviate from the tabulated numbers by typically a few percent
              for the effective temperature and 10 percent for the
              luminosities (not their logarithm).
    
  """
  
  def __init__(self):
    self._typeLetter = ["O", "B", "A", "F", "G", "K", "M"]
    # Generate all valid types (without LK)
    self._types = [l+str(n) for l in self._typeLetter for n in list(smo.range(10))]
    # Adjust to dJ87
    self._types.pop(0)
    self._types.append("M10")
    
    # Table to resolve variable s (Table 1)
    self._stable = [("O1", 0.1, 0.1),
                    ("O9", 0.9, 0.3),
                    ("B2", 1.8, 0.15),
                    ("A0", 3.0, 0.1),
                    ("F0", 4.0, 0.1),
                    ("G0", 5.0, 0.05),
                    ("K0", 5.5, 0.1),
                    ("M0", 6.5, 0.2)]
    # Table to resolve b (Table 2)
    self._btable = {"V":5., "IV":4., "III":3., "II":2., "Ib":1.4, "Iab":1.0, "Ia":0.6, "Ia+":0.0}
    # Table to resolve Chebychev polynomials (use self._cij[i,j])
    self._cij = np.array([[ 3.82573, -2.13868, -0.46357,  0.02076, -0.11937],
                          [-1.55607, -1.89216, -0.96916, -0.08869, -0.20423],
                          [ 1.05165,  0.42330, -0.94379, -0.07438, np.NaN],
                          [-0.01663, -0.20024, -0.18552, np.NaN,   np.NaN],
                          [-0.07576, -0.10934, np.NaN  , np.NaN,   np.NaN],
                          [ 0.11008, np.NaN  , np.NaN  , np.NaN,   np.NaN]])
    # Table to resolve Chebychev polynomials (use self._dij[i,j])
    self._dij = np.array([[ 3.96105,  0.03165, -0.02963,  0.01307, -0.01172],
                          [-0.62945,  0.02596, -0.06009,  0.01881, -0.01121],
                          [ 0.14370, -0.00977, -0.03265,  0.01649, np.NaN],
                          [ 0.00791,  0.00076, -0.03006, np.NaN,   np.NaN],
                          [ 0.00723, -0.02621, np.NaN  , np.NaN,   np.NaN],
                          [ 0.02755, np.NaN  , np.NaN  , np.NaN,   np.NaN]])
  
  def _indexSpecType(self, spt):
    """
      Get index of spectral type in list and check validity.
      
      Returns
      -------
      index : int
          Position in self._types list
    """
    if not spt in self._types:
      raise(PE.PyAValError("'" + str(spt) + "' is not a valid spectral type.", \
                           where="SpecTypeDeJager", \
                           solution="Use either of: " + ', '.join(self._types)))
    return self._types.index(spt)
  
  def _resolve_s(self, bt, st):
    """
      Resolve numerical value of variable s.
      
      Parameters
      ----------
      bt : string
          Spectral type (e.g., F3), no LK included
      st : float
          Subtype (e.g., 3.5 from F3.5)
      
      Returns
      -------
      s : float
          Value of s
    """
    spi = self._indexSpecType(bt)
    for i in smo.range(len(self._stable)-1, -1, -1):
      # Range start index
      rsi = self._indexSpecType(self._stable[i][0])
      if spi >= rsi:
        break
    s = self._stable[i][1] + float(spi+(st-int(st))-rsi)*self._stable[i][2]
    return s
      
  def _resolve_b(self, lk):
    """
      Resolve numerical value of variable b.
      
      Parameters
      ----------
      lk : string
          Luminosity class (e.g., IV)
      
      Returns
      -------
      b : float
          The value of b
    """
    if not lk in self._btable:
      raise(PE.PyAValError("'" + str(lk) + "' is not a valid luminosity class.", \
                           where="SpecTypeDeJager", \
                           solution="Use either of: " + ', '.join(list(self._btable))))
    return self._btable[lk]
  
  def _cheby(self, i, x):
    """
      Evaluate the Chebychev polynomial.
    """
    return np.cos(i*np.arccos(x))
  
  def _decomposeSPT(self, spt):
    """
      Decomposes spectral type into 'basic type' and subtype.
      
      Parameters
      ----------
      spt : string
          Spectral type (e.g, F8.5)
      
      Returns
      -------
      Basic type : string
          Type with integer subtype (e.g., F8)
      Subtype : float
          Subtype as a float (e.g., 8.5)
    """
    r = re.match("^\s*([OBAFGKM])(([\d])(\.\d+)?)\s*$", spt)
    if r is None:
      raise(PE.PyAValError("'" + str(spt) + "' is not a valid spectral type.", \
                           where="SpecTypeDeJager", \
                           solution="Use either of (also with float subtype): " + ', '.join(list(self._btable))))
    return r.group(1)+r.group(3), float(r.group(2))
  
  def lumAndTeff(self, spt, lk):
    """
      Get effective temperature and luminosity.
      
      Parameters
      ----------
      spt : string
          The spectral type (may include float subtype).
          For instance, F3, F3.5, or K2
      lk : string, {V, IV, III, II, Ib, Iab, Ia, Ia+}
          The luminosity class.
      
      Returns
      -------
      log10(L) : float
          The base-10 logarithm of the luminosity in units
          of the solar luminosity.
      log10(teff) : float
          The base-10 logarithm of the effective temperature.
    """
    bt, st = self._decomposeSPT(spt)
    s = self._resolve_s(bt, st)
    b = self._resolve_b(lk)    
    # Normalize to -1,+1 (Eqs. 2ab in dJ87)
    snorm = (s-4.25)/4.25
    bnorm = (b-2.5)/2.5
    
    # Prepare result
    llum = 0.0
    lteff = 0.0
    
    for n in smo.range(0,6):
      for i in smo.range(0, n+1):
        j = n - i
        t = self._cheby(i, snorm)*self._cheby(j, bnorm)
        try:
          llum += (self._cij[i,j] * t)
          lteff += (self._dij[i,j] * t)
        except:
          pass
    return llum, lteff