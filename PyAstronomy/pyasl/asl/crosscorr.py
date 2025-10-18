import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyasl import _ic


def crosscorrRV(
    w,
    f,
    tw,
    tf,
    rvmin,
    rvmax,
    drv,
    mode="doppler",
    skipedge=0,
    edgeTapering=None,
    weights=None,
    meanwvl=None,
):
    """
    Cross-correlate a spectrum with a template.

    The algorithm implemented here works as follows: For
    each RV shift to be considered, the wavelength axis
    of the template is shifted, either linearly or using
    a proper Doppler shift depending on the `mode`. The
    shifted template is then linearly interpolated at
    the wavelength points of the observation
    (spectrum) to calculate the cross-correlation function.

    If N data points :math:`f_i` at wavelengths :math:`w_i` are given, the cross
    correlation function :math:`CC(v_j)`, depending on the velocity :math:`v_j` and
    optional weights :math:`\\alpha_i` is calculated as:

    .. math::

        CC(v_j) = \\sum_{i=1}^N \\alpha_i \\times \left(f_i \\times t(w_i-\Delta_{i,j}) \\right)

    If the mode is `lin`, the shift is implemented as :math:`\Delta_{i,j} = \\bar{w} (v_j/c)`, where
    :math:`\\bar{w}` is the mean wavelength (i.e., the mean of the input w unless meanwvl is specified).
    For any velocity, the shift is a constant across the entire wavelength
    axis. If the mode is `doppler`, the shift is implemented as
    :math:`\Delta_{i,j} = w_i (v_j/c)`, which depends both on the wavelength and the velocity.
    The evaluation of the shifted template :math:`t(w_i-\Delta_{i,j})` is carried out using linear interpolation.
    
    If no weights are specified, :math:`\\alpha_i = 1` is adopted. If uncertainties, :math:`\\sigma_i`, for the
    :math:`f_i` are known, a useful choice for the weights is :math:`\\alpha_i = \sigma_i^{-2}`, because this makes
    maximizing the CCF equivalent to maximizing the likelihood 
    (see, e.g., `Lookwood et al. 2014, ApJ 783 <https://ui.adsabs.harvard.edu/abs/2014ApJ...783L..29L/abstract>`_).
    
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
    edgeTapering : float or tuple of two floats
        If not None, the method will "taper off" the edges of the
        observed spectrum by multiplying with a sine function. If a float number
        is specified, this will define the width (in wavelength units)
        to be used for tapering on both sides. If different tapering
        widths shall be used, a tuple with two (positive) numbers
        must be given, specifying the width to be used on the low- and
        high wavelength end. If a nonzero 'skipedge' is given, it
        will be applied first. Edge tapering can help to avoid
        edge effects (see, e.g., Gullberg and Lindegren 2002, A&A 390).
    weights : array of floats, optional
        If given, it specifies weights for individual data points (given by f).
        A possible choice for weights is 1/sigma_i**2 where sigma_i is the
        uncertainty of the i-th data point.
    meanwvl : float, optional
        If given, this value will be adopted for the mean wavelength
        :math:`\\bar{w}` to calculate shifts in 'lin' mode instead of
        the mean of the input array w.

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
        raise (
            PE.PyARequiredImport(
                "This routine needs scipy (.interpolate.interp1d).",
                where="crosscorrRV",
                solution="Install scipy",
            )
        )
    import scipy.interpolate as sci

    if weights is None:
        weights = np.ones_like(w)
    else:
        weights = weights.copy()

    if (w.shape != f.shape) or (w.shape != weights.shape):
        raise (
            PE.PyAValError(
                f"w, f, and weights must have the same shape. Given shapes are (w,f,weights): {w.shape}, {f.shape}, {weights.shape}",
                where="crosscorrRV",
                solution="Adjust arrays (w,f) and/or weights if specified.",
            )
        )

    # Copy and cut wavelength and flux arrays
    w, f = w.copy(), f.copy()
    if skipedge > 0:
        w, f = w[skipedge:-skipedge], f[skipedge:-skipedge]
        weights = weights[skipedge:-skipedge]

    if edgeTapering is not None:
        # Smooth the edges using a sine
        if isinstance(edgeTapering, (float, int)):
            edgeTapering = [float(edgeTapering), float(edgeTapering)]
        if len(edgeTapering) != 2:
            raise (
                PE.PyAValError(
                    "'edgeTapering' must be a float or a list of two floats.",
                    where="crosscorrRV",
                )
            )
        if edgeTapering[0] < 0.0 or edgeTapering[1] < 0.0:
            raise (
                PE.PyAValError(
                    "'edgeTapering' must be (a) number(s) >= 0.0.", where="crosscorrRV"
                )
            )
        # Carry out edge tapering (left edge)
        indi = np.where(w < w[0] + edgeTapering[0])[0]
        f[indi] *= np.sin((w[indi] - w[0]) / edgeTapering[0] * np.pi / 2.0)
        # Carry out edge tapering (right edge)
        indi = np.where(w > (w[-1] - edgeTapering[1]))[0]
        f[indi] *= np.sin(
            (w[indi] - w[indi[0]]) / edgeTapering[1] * np.pi / 2.0 + np.pi / 2.0
        )

    # Speed of light in km/s
    c = 299792.458
    # Check order of rvmin and rvmax
    if rvmax <= rvmin:
        raise (
            PE.PyAValError(
                "rvmin needs to be smaller than rvmax.",
                where="crosscorrRV",
                solution="Change the order of the parameters.",
            )
        )
    # Check whether template is large enough
    if mode == "lin":
        if meanwvl is None:
            meanWl = np.mean(w)
        else:
            meanWl = meanwvl
        dwlmax = meanWl * (rvmax / c)
        dwlmin = meanWl * (rvmin / c)
        if (tw[0] + dwlmax) > w[0]:
            raise (
                PE.PyAValError(
                    "The minimum wavelength is not covered by the template for all indicated RV shifts.",
                    where="crosscorrRV",
                    solution=["Provide a larger template", "Try to use skipedge"],
                )
            )
        if (tw[-1] + dwlmin) < w[-1]:
            raise (
                PE.PyAValError(
                    "The maximum wavelength is not covered by the template for all indicated RV shifts.",
                    where="crosscorrRV",
                    solution=["Provide a larger template", "Try to use skipedge"],
                )
            )
    elif mode == "doppler":
        # Ensure that the template covers the entire observation for all shifts
        maxwl = tw[-1] * (1.0 + rvmin / c)
        minwl = tw[0] * (1.0 + rvmax / c)
        if minwl > w[0]:
            raise (
                PE.PyAValError(
                    "The minimum wavelength is not covered by the template for all indicated RV shifts.",
                    where="crosscorrRV",
                    solution=["Provide a larger template", "Try to use skipedge"],
                )
            )
        if maxwl < w[-1]:
            raise (
                PE.PyAValError(
                    "The maximum wavelength is not covered by the template for all indicated RV shifts.",
                    where="crosscorrRV",
                    solution=["Provide a larger template", "Try to use skipedge"],
                )
            )
    else:
        raise (
            PE.PyAValError(
                "Unknown mode: " + str(mode),
                where="crosscorrRV",
                solution="See documentation for available modes.",
            )
        )
    # Calculate the cross correlation
    drvs = np.arange(rvmin, rvmax, drv)
    cc = np.zeros(len(drvs))
    for i, rv in enumerate(drvs):
        if mode == "lin":
            # Shift the template linearly
            fi = sci.interp1d(tw + meanWl * (rv / c), tf)
        elif mode == "doppler":
            # Apply the Doppler shift
            fi = sci.interp1d(tw * (1.0 + rv / c), tf)
        # Shifted template evaluated at location of spectrum
        cc[i] = np.sum(f * fi(w) * weights)
    return drvs, cc
