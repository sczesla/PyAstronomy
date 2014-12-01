import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyasl import _ic

def crosscorrRV(w, f, tw, tf, rvmin, rvmax, drv, mode="doppler", skipedge=0):
  """
    Cross-correlate a spectrum with a template.
    
    The algorithm implemented here works as follows: For
    each RV shift to be considered, the wavelength axis
    of the template is shifted, either linearly or using
    a proper Doppler shift depending on the `mode`. The
    shifted template is then linearly interpolated at
    the wavelength points of the observation
    (spectrum) to calculate the cross-correlation function.
    
    Parameters
    ----------
    w : array
        The wavelength axis of the observation.
    f : array
        The flux axis of the observation.
    tw : array
        The wavelength axis of the template.
    tf : array
        The flux axis of the template.
    rvmin : float
        Minimum radial velocity for which to calculate
        the cross-correlation function [km/s].
    rvmax : float
        Maximum radial velocity for which to calculate
        the cross-correlation function [km/s].
    drv : float
        The width of the radial-velocity steps to be applied
        in the calculation of the cross-correlation
        function [km/s].
    mode : string, {lin, doppler}, optional
        The mode determines how the wavelength axis will be
        modified to represent a RV shift. If "lin" is specified,
        a mean wavelength shift will be calculated based on the
        mean wavelength of the observation. The wavelength axis
        will then be shifted by that amount. If "doppler" is
        specified (the default), the wavelength axis will
        properly be Doppler-shifted.
    skipedge : int, optional
        If larger zero, the specified number of bins will be
        skipped from the begin and end of the observation. This
        may be useful if the template does not provide sufficient
        coverage of the observation.
    
    Returns
    -------
    dRV : array
        The RV axis of the cross-correlation function. The radial
        velocity refer to a shift of the template, i.e., positive
        values indicate that the template has been red-shifted and
        negative numbers indicate a blue-shift of the template.
        The numbers are given in km/s.
    CC : array
        The cross-correlation function.
  """
  if not _ic.check["scipy"]:
    raise(PE.PyARequiredImport("This routine needs scipy (.interpolate.interp1d).", \
                               where="crosscorrRV", \
                               solution="Install scipy"))
  import scipy.interpolate as sci
  # Copy and cut wavelength and flux arrays
  w, f = w.copy(), f.copy()
  if skipedge > 0:
    w, f = w[skipedge:-skipedge], f[skipedge:-skipedge]
  # Speed of light in km/s
  c = 299792.458
  # Check order of rvmin and rvmax
  if rvmax <= rvmin:
    raise(PE.PyAValError("rvmin needs to be smaller than rvmax.",
                         where="crosscorrRV", \
                         solution="Change the order of the parameters."))
  # Check whether template is large enough
  if mode == "lin":
    meanWl = np.mean(w)
    dwlmax = meanWl * (rvmax/c)
    dwlmin = meanWl * (rvmin/c)
    if (w.min() + dwlmin) < tw[0]:
      raise(PE.PyAValError("The minimum wavelength is not covered by the template.", \
                           where="crosscorrRV", \
                           solution=["Provide a larger template", "Try to use skipedge"]))
    if (w.max() + dwlmax) > tw[-1]:
      raise(PE.PyAValError("The maximum wavelength is not covered by the template.", \
                           where="crosscorrRV", \
                           solution=["Provide a larger template", "Try to use skipedge"]))
  elif mode == "doppler":
    maxwl = w[-1] * (1.0+rvmax/c)
    minwl = w[0] * (1.0+rvmin/c)
    if minwl < tw[0]:
      raise(PE.PyAValError("The minimum wavelength is not covered by the template.", \
                           where="crosscorrRV", \
                           solution=["Provide a larger template", "Try to use skipedge"]))
    if maxwl > tw[-1]:
      raise(PE.PyAValError("The maximum wavelength is not covered by the template.", \
                           where="crosscorrRV", \
                           solution=["Provide a larger template", "Try to use skipedge"]))
  else:
    raise(PE.PyAValError("Unknown mode: " + str(mode), \
                         where="crosscorrRV", \
                         solution="See documentation for available modes."))
  # Calculate the cross correlation
  drvs = np.arange(rvmin, rvmax, drv)
  cc = np.zeros(len(drvs))
  for i, rv in enumerate(drvs):
    if mode == "lin":
      # Shift the template linearly
      fi = sci.interp1d(tw+meanWl*(rv/c), tf)
    elif mode == "doppler":
      # Apply the Doppler shift
      fi = sci.interp1d(tw*(1.0 + rv/c), tf)
    # Shifted template evaluated at location of spectrum
    cc[i] = np.sum(f * fi(w))
  return drvs, cc