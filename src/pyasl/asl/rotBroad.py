from __future__ import print_function, division
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo

class _Gdl:
  
  def __init__(self, vsini, epsilon):
    """
      Calculate the broadening profile.
      
      Parameters
      ----------
      vsini : float
          Projected rotation speed of the star [km/s]
      epsilon : float
          Linear limb-darkening coefficient
    """
    self.vc = vsini / 299792.458
    self.eps = epsilon
  
  def gdl(self, dl, refwvl, dwl):
    """
      Calculates the broadening profile.
      
      Parameters
      ----------
      dl : array
          'Delta wavelength': The distance to the reference point in
          wavelength space [A].
      refwvl : array
          The reference wavelength [A].
      dwl : float
          The wavelength bin size [A].
      
      Returns
      -------
      Broadening profile : array
          The broadening profile according to Gray. 
    """
    self.dlmax = self.vc * refwvl
    self.c1 = 2.*(1.- self.eps) / (np.pi * self.dlmax * (1. - self.eps/3.))
    self.c2 = self.eps / (2.* self.dlmax * (1. - self.eps/3.))
    result = np.zeros(len(dl))
    x = dl/self.dlmax
    indi = np.where(np.abs(x) < 1.0)[0]
    result[indi] = self.c1*np.sqrt(1. - x[indi]**2) + self.c2*(1. - x[indi]**2)
    # Correct the normalization for numeric accuracy
    # The integral of the function is normalized, however, especially in the case
    # of mild broadening (compared to the wavelength resolution), the discrete
    # broadening profile may no longer be normalized, which leads to a shift of
    # the output spectrum, if not accounted for.
    result /= (np.sum(result) * dwl)
    return result
  

def rotBroad(wvl, flux, epsilon, vsini, edgeHandling="firstlast"):
  """
    Apply rotational broadening to a spectrum.
    
    This function applies rotational broadening to a given
    spectrum using the formulae given in Gray's "The Observation
    and Analysis of Stellar Photospheres". It allows for
    limb darkening parameterized by the linear limb-darkening law.
    
    The `edgeHandling` parameter determines how the effects at
    the edges of the input spectrum are handled. If the default
    option, "firstlast", is used, the input spectrum is internally
    extended on both sides; on the blue edge of the spectrum, the
    first flux value is used and on the red edge, the last value
    is used to extend the flux array. The extension is neglected
    in the return array. If "None" is specified, no special care
    will be taken to handle edge effects.
    
    .. note:: Currently, the wavelength array as to be regularly
              spaced.
    
    Parameters
    ----------
    wvl : array
        The wavelength array [A]. Note that a
        regularly spaced array is required.
    flux : array
        The flux array.
    vsini : float
        Projected rotational velocity [km/s].
    epsilon : float
        Linear limb-darkening coefficient (0-1).
    edgeHandling : string, {"firstlast", "None"}
        The method used to handle edge effects.
    
    Returns
    -------
    Broadened spectrum : array
        An array of the same size as the input flux array,
        which contains the broadened spectrum.
  """
  # Check whether wavelength array is evenly spaced
  sp = wvl[1::] - wvl[0:-1]
  if abs(max(sp) - min(sp)) > 1e-6:
    raise(PE.PyAValError("Input wavelength array is not evenly spaced.",
                         where="pyasl.rotBroad",
                         solution="Use evenly spaced input array."))
  if vsini <= 0.0:
    raise(PE.PyAValError("vsini must be positive.", where="pyasl.rotBroad"))
  if (epsilon < 0) or (epsilon > 1.0):
    raise(PE.PyAValError("Linear limb-darkening coefficient, epsilon, should be '0 < epsilon < 1'.",
                         where="pyasl.rotBroad",
                         solution="Adapt epsilon."))
  
  # Wavelength binsize
  dwl = wvl[1] - wvl[0]
  
  # Indices of the flux array to be returned
  validIndices = None
  
  if edgeHandling == "firstlast":
    # Number of bins additionally needed at the edges 
    binnu = int(np.floor(((vsini / 299792.458) * max(wvl)) / dwl)) + 1
    # Defined 'valid' indices to be returned
    validIndices = np.arange(len(flux)) + binnu
    # Adapt flux array
    front = np.ones(binnu) * flux[0]
    end = np.ones(binnu) * flux[-1]
    flux = np.concatenate( (front, flux, end) )
    # Adapt wavelength array
    front = (wvl[0] - (np.arange(binnu) + 1) * dwl)[::-1]
    end = wvl[-1] + (np.arange(binnu) + 1) * dwl
    wvl = np.concatenate( (front, wvl, end) )
  elif edgeHandling == "None":
    validIndices = np.arange(len(flux))
  else:
    raise(PE.PyAValError("Edge handling method '" + str(edgeHandling) + "' currently not supported.",
                         where="pyasl.rotBroad",
                         solution="Choose ones of the valid edge handling methods"))
    
  
  result = np.zeros(len(flux))
  gdl = _Gdl(vsini, epsilon)
  
  for i in smo.range(len(flux)):
    dl = wvl[i] - wvl
    g = gdl.gdl(dl, wvl[i], dwl)
    result[i] = np.sum(flux * g)
  result *= dwl
  
  return result[validIndices]


def fastRotBroad(wvl, flux, epsilon, vsini, effWvl=None):
  """
    Apply rotational broadening using a single broadening kernel.
    
    The effect of rotational broadening on the spectrum is
    wavelength dependent, because the Doppler shift depends
    on wavelength. This function neglects this dependence, which
    is weak if the wavelength range is not too large.
    
    .. note:: numpy.convolve is used to carry out the convolution
              and "mode = same" is used. Therefore, the output
              will be of the same size as the input, but it
              will show edge effects.
    
    Parameters
    ----------
    wvl : array
        The wavelength
    flux : array
        The flux
    epsilon : float
        Linear limb-darkening coefficient
    vsini : float
        Projected rotational velocity in km/s.
    effWvl : float, optional
        The wavelength at which the broadening
        kernel is evaluated. If not specified,
        the mean wavelength of the input will be
        used.
    
    Returns
    -------
    Broadened spectrum : array
        The rotationally broadened output spectrum.
  """
  # Check whether wavelength array is evenly spaced
  sp = wvl[1::] - wvl[0:-1]
  if abs(max(sp) - min(sp)) > 1e-6:
    raise(PE.PyAValError("Input wavelength array is not evenly spaced.",
                         where="pyasl.rotBroad",
                         solution="Use evenly spaced input array."))
  if vsini <= 0.0:
    raise(PE.PyAValError("vsini must be positive.", where="pyasl.rotBroad"))
  if (epsilon < 0) or (epsilon > 1.0):
    raise(PE.PyAValError("Linear limb-darkening coefficient, epsilon, should be '0 < epsilon < 1'.",
                         where="pyasl.rotBroad",
                         solution="Adapt epsilon."))
  
  # Wavelength binsize
  dwl = wvl[1] - wvl[0]
  
  if effWvl is None:
    effWvl = np.mean(wvl)
  
  gdl = _Gdl(vsini, epsilon)
  
  # The number of bins needed to create the broadening kernel
  binnHalf = int(np.floor(((vsini / 299792.458) * effWvl / dwl))) + 1
  gwvl = (np.arange(4*binnHalf) - 2*binnHalf) * dwl + effWvl
  # Create the broadening kernel
  dl = gwvl - effWvl
  g = gdl.gdl(dl, effWvl, dwl)
  # Remove the zero entries
  indi = np.where(g > 0.0)[0]
  g = g[indi]
  
  result = np.convolve(flux, g, mode="same") * dwl
  return result