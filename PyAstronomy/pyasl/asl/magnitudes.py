from PyAstronomy.pyaC import pyaErrors as PE
import numpy as np


def magToFluxDensity_bessel98(band, mag, mode="nu"):
    """
    Convert magnitude into flux density according to Bessel et al. 1998

    The conversion implemented here is based on the data given in Table A2 of
    Bessel et al. 1998, A&A 333, 231-250, which gives "Effective wavelengths (for an A0 star),
    absolute fluxes (corresponding to zero magnitude) and zeropoint magnitudes for the UBVRI-
    JHKL Cousins-Glass-Johnson system". Note that zp(f_nu) and zp(f_lam) are exchanged in
    the original table.

    Parameters
    ----------
    band : string
        Any of U, B, V, R, I, J, H, K, Kp, L, and L*
    mag : float, array
        The magnitude value to be converted
    mode : string, {nu, mod}
        Determines whether f_nu or f_lam will be calculated.

    Returns
    -------
    f_nu/lam : float
        The corresponding flux density in units if erg/cm**2/s/Hz in the
        case of mode 'nu' and erg/cm**2/s/A in the case of 'lam'.
    lam_eff : float
        Effective filter wavelength in Angstrom
    """
    if not mode in ["nu", "lam"]:
        raise (
            PE.PyAValError(
                "Unknown mode",
                where="magToFluxDensity_bessel98",
                solution="Use 'nu' or 'lam'",
            )
        )

    # Data from Bessel 1998 Table A2
    d = """bands U B V R I J H K Kp L L*
    lam_eff 0.366 0.438 0.545 0.641 0.798 1.22 1.63 2.19 2.12 3.45 3.80
    f_nu 1.790 4.063 3.636 3.064 2.416 1.589 1.021 0.640 0.676 0.285 0.238
    f_lam 417.5 632 363.1 217.7 112.6 31.47 11.38 3.961 4.479 0.708 0.489
    zp(f_nu) 0.770 -0.120 0.000 0.186 0.444 0.899 1.379 1.886 1.826 2.765 2.961
    zp(f_lam) -0.152 -0.602 0.000 0.555 1.271 2.655 3.760 4.906 4.780 6.775 7.177"""
    # Digest table
    s = d.split()
    t = {}
    for i in range(6):
        t[s[i * 12]] = s[i * 12 + 1 : (i + 1) * 12]
        if s[i * 12] != "bands":
            t[s[i * 12]] = np.array(t[s[i * 12]], dtype=float)

    if not band in t["bands"]:
        raise (
            PE.PyAValError(
                "Unknown passband: " + str(band),
                where="magToFluxDensity_bessel98",
                solution="Use any of: " + ", ".join(t["bands"]),
            )
        )

    # Band index
    bi = t["bands"].index(band)

    c = {"nu": 48.598, "lam": 21.1}

    f = 10.0 ** ((mag + c[mode] + t["zp(f_" + mode + ")"][bi]) / -2.5)
    return f, t["lam_eff"][bi] * 1e4


def absMagToPower(am, absMagSun=4.75, absLumSun=3.846e33):
    """
    Convert absolute magnitude to power scale

    The default values for the absolute magnitude and luminosity
    of the Sun are adopted from Harmanec and Prsa 2011
    (2011PASP..123..976H).

    Parameters
    ----------
    am : float
        Absolute magnitude.
    absMagSun : float, optional
        Absolute magnitude of the Sun.
    absLumSun : float, optional
        Absolute luminosity of the Sun.
        The default is given in units of erg/s.

    Returns
    -------
    Power : float
        Total emitted power. Same units as `absLumSun`;
        the default corresponds to erg/s.
    """
    power = 10.0 ** ((am - absMagSun) / (-2.5)) * absLumSun
    return power


def absModuleToDist(magApp, magAbs):
    """
    Convert apparent and absolute magnitude into distance.

    Parameters
    ----------
    magApp : float
        Apparent magnitude of object.
    magAbs : float
        Absolute magnitude of object.

    Returns
    -------
    Distance : float
        The distance resulting from the difference in
        apparent and absolute magnitude [pc].
    """
    d = 10.0 ** (-(magAbs - magApp) / 5.0 + 1.0)
    return d
