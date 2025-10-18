from PyAstronomy import funcFit as fuf
import numpy as np


class VoigtAstroP(fuf.OneDFit):
    """
    Astronomical parameterization of Voigt profile.

    This class provides a convenient parameterization of the Voigt
    profile, frequently used in astronomy. In particular,
    the line is parameterized in terms of wavelength,
    Doppler parameter, damping width, and oscillator strength.
    
    If the velocity components of the constituents follow
    independent Gaussian distributions with standard deviation
    :math:`\\sigma` in all three dimensions of space
    
    .. math::
        
        v_{x,y,z} \\sim N(0,\sigma^2) \; ,
        
    the mode of the distribution of absolute velocity,
    :math:`v_a = \\sqrt{v_x^2 + v_y^2 + v_z^2}`
    is given by :math:`v_{mode} = \sqrt{2}\\sigma`. This is
    the velocity referred to as the Doppler parameter, b.
    Consequently, the velocity dispersion in any individual dimension,
    i.e., the standard deviation of the Gaussian velocity distribution,
    is related to the Doppler parameter by :math:`\\sigma = b/\\sqrt{2}`.
    For zero damping width (gamma),
    the resulting model line is, therefore, a Gaussian with a standard
    deviation of :math:`\\sigma = b/\\sqrt{2}` in velocity units.
    
    Instrumental resolution can be applied via the parameter R.
    The instrumental profile is assumed to be a Gaussian with
    FWHM of w0/R. The additional broadening is implemented by
    using an internal, effective Doppler parameter. The case
    R=0 corresponds to infinite instrumental resolution (i.e.,
    no additional broadening).

    *Fit parameters*:

    ======= ================================================ =====

    w0      Wavelength of the transition                     A
    b       Doppler parameter (corresponds                   km/s
            to sqrt(2) times the velocity dispersion).
    gamma   Damping width (full width at half maximum of
            the Lorentzian)                                  cm
    f       Oscillator strength (unitless)                   --
    R       Instrumental resolution                          --
    ======= ================================================ =====

    """

    def __init__(self):
        fuf.OneDFit.__init__(
            self, ["w0", "b", "gamma", "f", "lin", "off", "R"], rootName="VoigtAstroP")
        self._profile = fuf.Voigt1d()
        # Set default parameter, which may not be changed
        self._profile["ad"] = 1. / np.sqrt(2.0)
        # Set to zero, because evaluation is done in velocity space
        self._profile["mu"] = 0.0
        # Conversion factor between Gaussian FWHM and STD
        self._fwhmstd = 2*np.sqrt(2*np.log(2.))

    def bl(self):
        """
        Doppler width in cm
        """
        bl = self["w0"] / 1e8 * (self["b"] * 100000.0) / 29979245800.0
        if self["R"] > 0:
            # Account for instrumental resolution
            bl = np.sqrt(bl**2 + 2*(self["w0"]/(1e8*self["R"]*self._fwhmstd))**2)
        return bl

    def FWHM(self):
        """
        Estimate FWHM
        
        Applies same approximation is Voigt1d
        
        Returns
        -------
        FWHM : float
            FWHM of line profile in wavelength units [A]
        """
        fwhm = self._profile.FWHM()
        c = 1e8 * self.bl()
        return fwhm * c
            
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
        bl = self.bl()
        # The following formulation is, e.g., according to D.F. Gray "Stellar photospheres"
        # (third edition, Eq. 11.46, p. 257). For the relation between H(u, a) and the Voigt
        # profile see Pagnini and Saxena "A note on the Voigt profile" (Eq. 6,
        # https://arxiv.org/abs/0805.2274); the relation quoted by Gray is probably not entirely
        # accurate. Note that the Voigt profile is here evaluated at u = dw/bl.
        # 
        # The constant equals (pi e^2)/(m_e c^2) (pi times the classical electron radius)
        self._profile["A"] = 8.85282064473e-13 * self["f"] * w0cm**2 / bl
        # A factor of 2.0 because `al` defines the half FWHM in Voigt profile (division by bl
        # accounts for evaluation at dw/bl).
        self._profile["al"] = (self["gamma"] / bl) / 2.0
        self._profile["lin"] = self["lin"]
        self._profile["off"] = self["off"]
        # Convert the given wavelength axis into velocity units
        u = (x / 1e8 - w0cm) / bl
        return self._profile.evaluate(u)
