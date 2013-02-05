# -*- coding: utf-8 -*-
import numpy as np

def flux2photons(wave, flux):
  """
    Convert flux (erg/s) to photons/s.
    
    Parameters
    ----------
    wave : float, array
        The wavelength in Angstrom.
    flux : float, array
        The associated flux array in erg/s.
   
    Returns
    -------
    photons : array
        Number of photons/s for the given input flux

  """
  # h*c = 1.98644568327e-08 erg*A; ephoton = energy(erg) per photon
  ePhoton = 1.98644568327e-08 / np.array(wave)
  nPhotons = np.array(flux) / ePhoton
  return nPhotons
  
def photons2flux(wave, photons):
  """
    Convert photons/s to flux (erg/s).
    
    Parameters
    ----------
    wave : float, array
        The wavelength in Angstrom.
    photons : float, array
        The associated photons array in photons/s.
   
    Returns
    -------
    flux : array
        Flux in erg/s for the given photon flux
  """
  # h*c = 1.98644568327e-08 erg*A; ephoton = energy(erg) per photon
  ePhoton = 1.98644568327e-08 / np.array(wave)
  flux = np.array(photons) * ePhoton
  return flux  