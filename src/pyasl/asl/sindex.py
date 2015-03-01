from __future__ import print_function 
from PyAstronomy.pyaC import pyaErrors as PE
import numpy as np

class SMW_RHK:

  def __init__(self, ccfs="rutten", afc="middelkoop", rphot="noyes"):
    """
      Converting Mount-Wilson S-index into RHK index.
      
      The Mount-Wilson S-index is a measure of the emission-line cores
      of the Ca II H and K lines at about 3933 A and 3968 A in two
      narrow bands normalized by two adjacent continuum bands.
      
      The activity index RHK is closely related to the S-index. In
      particular, it gives the emission in the narrow bands normalized
      by the bolometric brightness of the star
      
      .. math::
      
          R_{HK} = \\frac{4\\pi R_s^2 (F_H + F_K)}{4\\pi R_s^2\\sigma T_{eff}^4} = \\frac{F_H+F_K}{\\sigma T_{eff}^4} \\; .
      
      The stellar surface flux in "arbitrary units" in the narrow
      H and K bands (fH + fK) is related to the Mount-Wilson S-index
      through the relation
      
      .. math::
      
          f_H + f_K = S C_{cf} T_{eff}^4 10^{-14} \\; ,
      
      where Ccf is a conversion factor, which can be parameterized in terms of
      the B-V color and the luminosity class. The conversion between
      arbitrary units and physical units, needed to derive the true surface
      flux and the RHK index, has been derived by several authors
      starting with Middelkoop 1982.
      Their factor was also used by Noyes et al. 1984---in particular in their
      appendix a, where is appears implicitly. Later, the value of the conversion
      factor has been revised by several authors, e.g., Oranje 1983 and Rutten 1984,
      who estimated a value about 70% larger than previously proposed. Hall et al.
      2007 derive a value 40% larger than that of Middelkoop 1982 and provide
      a thorough discussion on the differences between the individual approaches
      and results given in the literature.
      
      Finally, the RHK index thus derived still covers a photospheric
      contribution, which is always present and not related to the
      chromosphere. To obtain the purely chromospheric, primed RHK index,
      an estimate of the photospheric surface flux in the H and K pass-bands
      has to be subtracted. For active stars, the photospheric correction
      is usually quite irrelevant. For inactive, quiet stars, it can,
      however, be important.
      
      The issue of the Mount-Wilson S-index conversion has been revisited
      by Mittag et al. 2013, who provide an alternative conversion procedure
      and revised photospheric corrections for various luminosity classes. 
      
      .. note:: In the default configuration, the conversion of the S-index
                into RHK is identical to the relation stated by Noyes et al. 1984
                in their appendix (a)
                
                .. math::
                
                    R_{HK} = 1.340 \\times 10^{-4} C_{cf} S
                
                where the factor 1.34e-4 is a combination of the conversion from
                arbitrary to physical units, 1e-14, and the Stefan-Boltzmann
                constant, in particular 1.34e-4 = 7.6e5*1e-14/5.67e-5. The Ccf factor
                is, however, calculated according to Rutten 1984.
      
      The relations and coefficients used here are taken from the
      following publications (and references therein):
        - Middelkoop 1982, A&A 107, 31
        - Oranje 1983, A&A 124, 43
        - Noyes et al. 1984, A&A 279, 763
        - Rutten 1984, A&A 130, 353
        - Hall et al. 2007, AJ 133, 862
        - Mittag et al. 2013, A&A 549, 117
      
      Parameters
      ----------
      ccfs : string, {rutten, noyes}, optional
          Source of the conversion factor between S-index and RHK.
      afc : string, {rutten, oranje, middelkoop, hall}, optional
          Source of conversion factor between "arbitrary units"
          and physical units of surface flux.
      rphot : string, {noyes}
          The source for the photospheric correction for the RHK
          index.
    """
    if ccfs == "rutten":
      self._log10ccf = self.log10ccfRutten
    elif ccfs == "noyes":
      self._log10ccf = self.log10ccfNoyes
    else:
      raise(PE.PyAValError("No such Ccf source: " + str(ccfs), \
                           solution="Use 'rutten' or 'noyes'."))
    self._ccfs = ccfs
    
    if rphot == "noyes":
      self._logrphot = self.logRphotNoyes
    else:
      raise(PE.PyAValError("No such source for the photospheric correction: " + str(rphot), \
                           solution="Use 'noyes'."))
    self._rphots = rphot
    
    self._afc = afc
    # The conversion factor from "arbitrary units" to physical units.
    self._absCal = {"rutten":1.29e6, "oranje":1.21e6, "middelkoop": 7.6e5, "hall":1.07e6}
    if not self._afc in self._absCal:
      raise(PE.PyAValError("No such source for the conversion from arbitrary to physical units: " + str(self._afc), \
                           solution="Use either of: " + ', '.join(self._absCal.keys()) ))
    
    from PyAstronomy.pyasl import _ic
    if not _ic.check["quantities"]:
      raise(PE.PyARequiredImport("The 'quantities' package is not installed, which is required to use 'SMW_RHK'.", \
                                 where="SMW_RHK", \
                                 solution="Install quantities (https://pypi.python.org/pypi/quantities/0.10.1)."))
    from PyAstronomy import constants as PC
    self._pc = PC.PyAConstants()
    self._pc.setSystem("cgs")
  
  def logRphotNoyes(self, bv, lc="ms"):
    """
      Photospheric contribution to surface flux in the H and K pass-bands.
      
      Relation given by Noyes et al. 1984.
      
      Parameters
      ----------
      bv : float
          B-V color [mag]
      lc : string, {ms, g}, optional
          Luminosity class.
      
      Returns
      -------
      log10(Rphot) : float
          Logarithm of the photospheric contribution.
    """
    if (bv < 0.44) or (bv > 0.82):
      PE.warn(PE.PyAValError("Noyes et al. 1984 give a validity range of 0.44 < B-V < 0.82 for the " + \
                             "photospheric correction. However, the authors use it for B-V > 0.82, " + \
                             "where it quickly decreases."))
    if lc != "ms":
      PE.warn(PE.PyAValError("Noyes et al. 1984 specify the photospheric correction only for main-sequence stars."))
    rp = -4.898 + 1.918*bv**2 - 2.893*bv**3
    return rp
      
  def log10ccfNoyes(self, bv, **kwargs):
    """
      Ccf conversion factor according to Noyes et al. 1984.
      
      Parameters
      ----------
      bv : float
          The B-V color [mag].
      
      Returns
      -------
      log10(Ccf) : float
          The logarithm of the conversion factor.
    """
    if ("lc" in kwargs) and (kwargs["lc"] != "ms"):
      PE.warn(PE.PyAValError("The Ccf conversion factor by Noyes et al. 1984 is only valid for main-sequence stars.", \
              solution="Use the conversion factor by Rutten 1984 for giants."))
    logccf = 1.13*bv**3 - 3.91*bv**2 + 2.84*bv - 0.47
    if bv <= 0.63:
      x = 0.63 - bv
      dlogccf = 0.135*x - 0.814*x**2 + 6.03*x**3
      logccf += dlogccf
    return logccf
    
  def log10ccfRutten(self, bv, lc="ms"):
    """
      Ccf conversion factor from Rutten 1984 (Eqs. 10a and 10b).
      
      Parameters
      ----------
      bv : float
          B - V color [mag].
      lc : string, {ms, g}, optional
          Specifies whether the relation for
          main-sequence (ms) or giant (g) stars
          shall be evaluated.
      
      Returns
      -------
      log10(Ccf) : float
          The logarithm of the conversion factor.
    """
    if lc == "ms":
      if (bv < 0.3) or (bv > 1.6):
        PE.warn(PE.PyAValError("B-V color out of range. Rutten 1984 states a validity range of 0.3 <= b-v <= 1.6 " +
                               "for main-sequence stars. You specified: " + str(bv) + "."))
      logccf = 0.25*bv**3 - 1.33*bv**2 + 0.43*bv + 0.24
    elif lc == "g":
      if (bv < 0.3) or (bv > 1.7):
        PE.warn(PE.PyAValError("B-V color out of range. Rutten 1984 states a validity range of 0.3 <= b-v <= 1.7 " +
                               "for giant stars. You specified: " + str(bv) + "."))
      logccf = -0.066*bv**3 - 0.25*bv**2 - 0.49*bv + 0.45
    else:
      raise(PE.PyAValError("No such luminosity class: " + str(lc), \
                           solution="Specify either 'ms' or 'g'."))
    return logccf
  
  def FHFK(self, S, Teff, log10ccf):
    """
      Calculate the FH+FK flux in arbitrary units.
      
      Parameters
      ----------
      S : float
          Mount-Wilson S-index.
      Teff : float
          The effective temperature [K].
      log10ccf : float
          The logarithm of the Ccf conversion factor.
      
      Returns
      -------
      FH + FK : float
          The stellar surface flux in the H and K pass-bands
          in arbitrary units (not erg/cm**2/s). 
    """
    ccf = 10.0**log10ccf
    fhfk = S * ccf * Teff**4 * 1e-14
    return fhfk
  
  def SMWtoRHK(self, S, Teff, bv, lc="ms", verbose=False):
    """
      Convert Mount-Wilson S-index into R_HK.
      
      Parameters
      ----------
      S : float
          Mount-Wilson S-index.
      Teff : float
          Effective temperature [K].
      bv : float
          B-V color [mag]
      lc : String, {ms, g}, optional
          Luminosity class; Main-sequence (ms) or giants (g)
      verbose : boolean, optional
          If True, the details of the calculation are printed
          to stdout.
      
      Returns
      -------
      RHK prime : float
          RHK parameter corrected for photospheric contribution. The primed
          number measures the purely chromospheric emission.
      RHK : float
          RHK parameter without correction for photospheric contribution.
      ccf : float
          The Ccf conversion factor used.
      fhfk : float
          The FH+FK surface flux in arbitrary units.
      fhfk (physical) : float
          The FH+FK surface flux in physical units [erg/cm^2/s].
      R_phot : float
          Photospheric flux contribution used in translating RHK into
          RHK prime.
    """
    # Get Ccf conversion factor
    log10ccf = self._log10ccf(bv, lc=lc)
    ccf = 10.0**log10ccf
    # Get FH+FK
    fhfk = self.FHFK(S, Teff, log10ccf)
    # Convert arbitrary units to physical units
    surfaceFlux = fhfk * self._absCal[self._afc]
    # Get RHK (includes photospheric contribution)
    rhk = surfaceFlux/(self._pc.sigma * Teff**4)
    # Get the photospheric correction
    logrphot = self._logrphot(bv, lc=lc)
    # Correct RHK for photospheric contribution
    rhkprime = rhk - 10.0**logrphot
    if verbose:
      print("Converting Mount-Wilson S-index to RHK")
      print("  Source of photospheric correction: " + self._rphots)
      print("  Source of Ccf conversion factor: " + self._ccfs)
      print("  log10ccf = %6.3e, ccf = %6.3e" % (log10ccf, ccf))
      print("  Surface flux in H and K pass-bands in arbitrary units: %6.3e" % (fhfk))
      print("  Arbitrary unit to flux conversion factor: %6.3e" % (self._absCal[self._afc]) + \
            " from source: " + self._afc)
      print("  Surface flux in physical units [erg/cm^2/s]: %6.3e" % (surfaceFlux))
      print("  R_HK (including photosphere): %6.3e" % (rhk))
      print("  log10(R_HK) (including photosphere): %6.3e" % (np.log10(rhk)))
      print("  Photospheric contribution (log10(R_phot)): %6.3e" % logrphot)
      print("  R_HK prime (corrected for photospheric correction): %6.3e" % (rhkprime))
      print("  log10(R_HK prime) (corrected for photospheric correction): %6.3e" % ((np.log10(rhkprime))))
    return rhkprime, rhk, ccf, fhfk, surfaceFlux, 10.0**logrphot
    