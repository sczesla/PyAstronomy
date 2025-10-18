from __future__ import division
from . import pal
import PyAstronomy.funcFit as fuf
from PyAstronomy.pyaC import pyaErrors as PE
import six

class SecPalLC(fuf.OneDFit):
  """
    Calculate and fit analytical secondary transit light-curves using the formulae
    provided by Pal 2008.
    
    .. note :: The **evaluation of elliptical integrals** is essential in
               calculating the transit model. While both
               the *mpmath* module and the *Boost* libraries implement those
               integrals, it is the Boost library, which
               evaluates them by far more quickly. Yet, the support
               for Boost has to be added manually. 
    
    *Model parameters*:
      - `p` - Radius ratio between planet and star.
      - `a` - Semi-major axis of planetary orbit [stellar radii].
      - `i` - Inclination of orbit in degrees (90 deg is *edge on* view).
      - `T0` - Time offset of secondary transit center.
      - `per` - Period of planetary orbit.
      - `b` - Describes the flux ratio between a stellar companion and the main star (default is 0).
      - `brat` - Brightness ratio between planet and star (one mean that planet has the same brightness as the star).

    @TODO - subtractOne docu

    This class inherits the functionality of funcFit's OneDFit object.
    You can, therefore, set and get the parameter using the brackets:
      e.g., pallc["p"] = 0.12345

    .. warning::
        Time units have to be consistent.
  """
  
  def __init__(self, subtractOne=True):
    fuf.OneDFit.__init__(self, ["p", "a", "i", "T0", "per", "brat"  ,"b"])
    self.freeze(["p", "a", "i", "T0", "per", "b", "brat"])
    self.pal = pal.PalLC()
    self.pal["quadLimb"] = 0.0
    self.pal["linLimb"] = 0.0
    self.pal["b"] = 0.0
    self.subtractOne = subtractOne
    self.setRootName("SecPal")

  def evaluate(self, time):
    """ 
     Calculate a light curve for secondary transit.
        
     Parameters
     ----------
     time : array
         An array of time points at which the light curve
         shall be calculated.
        
     .. note:: time = 0 -> Planet is exactly in the line of sight (phase = 0).

     Returns
     -------
     Model : array
         The analytical light curve.
    """
    # Assign all variables
    for k in six.iterkeys(self.propMap):
      if k == "brat" or k == "b":
        continue
      self.pal[k] = self[self.propMap[k]]
    # Calculate LC
    df = ((1.0 - self.pal.evaluate(time)) / self["p"]**2 * self["brat"])
    # Account for possible stellar companion
    y = 1.0 - ( df / (1.0 + self["b"]) )
    if self.subtractOne:
      y -= 1.0
    return y
