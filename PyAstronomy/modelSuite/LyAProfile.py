from .voigtAstro import VoigtAstroP
from PyAstronomy import funcFit as fuf
from PyAstronomy.pyasl import convertDampingConstant
import numpy as np

class LyaTransmission(fuf.OneDFit):
  """
    Lyman alpha transmission profile including Deuterium absorption.
    
    The transmission is given by
      
    .. math::
       
        e^{-\\sigma N}\\, ,
        
    where N is the column density and :math:`\\sigma` is the wavelength-dependent
    cross-section.
    
    *Fit Parameters*:
    
    ===== ============================= =====
    N     Hydrogen column density       /cm^2
    b     Doppler parameter             km/s
    Dfrac Deuterium fraction            --
    ===== ============================= =====
    
    Parameters
    ----------
    N : float, optional
        Hydrogen column density [/cm^2]. The default is 0.
    b  : float, optional
        The Doppler parameter (corresponds to sqrt(2) times
        the velocity dispersion) to model thermal line width [km/s].
        The default is 10 km/s.
    D_fraction : float, optional
        Fractional abundance of Deuterium with respect to Hydrogen.
        The default is 1.5e-5.
  """
  
  def __init__(self, N=0.0, b=10.0, D_fraction=1.5e-5):
    fuf.OneDFit.__init__(self, ["N","b", "Dfrac"], rootName="LyaTransmission")
    self["b"] = b
    self["Dfrac"] = D_fraction
    # (Only) Einstein coefficient relevant for LyA
    elya = 6.258085e8
    # Hydrogen
    self._absH = VoigtAstroP()
    self._absH["w0"] = 1215.67
    self._absH["b"] = self["b"]
    self._absH["gamma"] = convertDampingConstant(elya, 1215.67)
    self._absH["f"] = 0.416
    # Deuterium, its width is sqrt(2) times smaller
    self._absD = VoigtAstroP()
    self._absD["w0"] = 1215.34
    self._absD["b"] = self["b"]/np.sqrt(2.0)
    self._absD["gamma"] = convertDampingConstant(elya, 1215.67)
    self._absD["f"] = 0.416
  
  def evaluate(self, x):   
    """
      Evaluate the transmission profile.
      
      Parameters
      ----------
      x : array of floats
          Contains the wavelength values in Angstrom.
      
      Returns
      -------
      Model : array of floats
          The line profile.
    """
    self._absH["b"] = self["b"]
    self._absD["b"] = self["b"]/np.sqrt(2.0)
    y = np.exp(-self["N"] * self._absH.evaluate(x)) * \
        np.exp(-self["N"] * self["Dfrac"] * self._absD.evaluate(x))
    return y
