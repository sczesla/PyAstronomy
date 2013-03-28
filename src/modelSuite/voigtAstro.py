from PyAstronomy import funcFit as fuf
import numpy as np

class VoigtAstroP(fuf.OneDFit):
    """
      Astronomical parameterization of Voigt profile.
      
      This class provides a convenient parameterization of the Voigt
      profile, as it is frequently used in astronomy. In particular,
      the line is parameterized in terms of wavelength,
      Doppler parameter, damping width, and oscillator strength.
      
      *Fit parameters*:
      
      ======= ================================================ =====
      
      w0      Wavelength of the transition                     A
      b       Doppler parameter (corresponds                   km/s
              to sqrt(2) times the velocity dispersion).
      gamma   Damping width (full width at half maximum of
              the Lorentzian)                                  cm
      f       Oscillator strength (unitless)                   --
      ======= ================================================ =====
      
    """
    
    def __init__(self):
      fuf.OneDFit.__init__(self, ["w0", "b", "gamma", "f", "lin", "off"], rootName="VoigtAstroP")
      self._profile = fuf.Voigt1d()
      # Set default parameter, which may not be changed
      self._profile["ad"] = 1./np.sqrt(2.0)
      # Set to zero, because evaluation is done in velocity space
      self._profile["mu"] = 0.0
    
    def evaluate(self, x):
      """
        Evaluate the absorption-line profile.

        Parameters
        ----------
        x : array of floats
            Contains the wavelength values in Angstrom.
        
        Returns
        -------
        Model : array of floats
            Return the cross-section in cm^2.
      """
      # Wavelength in cm
      w0cm = self["w0"] / 1e8
      # Doppler width in cm
      bl = w0cm * (self["b"] * 100000.0) / 29979245800.0
      # The constant equals (pi e^2)/(m_e c^2)
      self._profile["A"] =  8.85282064473e-13 * self["f"]  * w0cm**2 / bl
      # A factor of 2.0 because `al` defines the half FWHM in Voigt profile
      self._profile["al"] = (self["gamma"] / bl) / 2.0
      self._profile["lin"] = self["lin"]
      self._profile["off"] = self["off"]
      # Convert the given wavelength axis into velocity units
      u = (x/1e8 - w0cm) / bl
      return self._profile.evaluate(u)