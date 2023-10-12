from __future__ import print_function, division
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import funcFit as fuf
import scipy.interpolate as sci
from quantities.dimensionality import assert_isinstance


def equidistantInterpolation(x, y, dxmode, ifct=sci.interp1d):
    """
    Construct evenly sampled data set by interpolation
    
    Parameters
    ----------
    x : array
        Input x-axis
    y : array or list/tuple of arrays
        The y axis or axes
    dxmode : string {2x, mean} or float
        If a float is given, it specifies the spacing of the new x axis. If '2x' is
        given, an equidistant x axis with twice as many points as the input one will
        be used. If 'mean' is given, the spacing of the new x axis corresponds to the
        mean of that of the input axis. If dxmode is an array, it is used as the new x-axis
        (without check for equidistant sampling).
    ifct : interpolation callable
        A callable taking x and y as arguments and returning a callable, which takes
        the new equidistant wavelength axis as argument. The default is interp1d from
        scipy called by ny = ifct(x, y)(nx).
    
    Returns
    -------
    x, y : arrays
        New equidistant x axis and linearly interpolated y axis. If input y is a list or tuple of arrays,
        a list of interpolated arrays is returned.
    """
    if type(dxmode) == float:
        nx = np.arange(x[0], x[-1], dxmode)
    elif isinstance(dxmode, np.ndarray):
        nx = dxmode
    elif dxmode == "2x":
        nx = np.linspace(x[0], x[-1], 2*len(x))
    elif dxmode == "mean":
        dx = np.mean(np.diff(x))
        nx = np.arange(x[0], x[-1], dx)
    else:
        raise(PE.PyAValError("Cannot interpret 'dxmode'", \
                             where="makeEquidistant", \
                             solution=["Use a float to specify spacing.", \
                                       "Use '2x' or 'mean' to give rule to select spacing"]))
    
    if isinstance(y, (list,tuple)):
        ny = [ifct(x, yy)(nx) for yy in y]
    else:
        ny = ifct(x, y)(nx)
    
    return nx, ny


def broadGaussFast(x, y, sigma, edgeHandling=None, maxsig=None):
    """
    Apply Gaussian broadening. 

    This function broadens the given data using a Gaussian
    kernel.

    Parameters
    ----------
    x, y : arrays
        The abscissa and ordinate of the data.
    sigma : float
        The width (i.e., standard deviation) of the Gaussian
        profile used in the convolution.
    edgeHandling : string, {None, "firstlast"}, optional
        Determines the way edges will be handled. If None,
        nothing will be done about it. If set to "firstlast",
        the spectrum will be extended by using the first and
        last value at the start or end. Note that this is
        not necessarily appropriate. The default is None.
    maxsig : float, optional
        The extent of the broadening kernel in terms of
        standard deviations. By default, the Gaussian broadening
        kernel will be extended over the entire given spectrum,
        which can cause slow evaluation in the case of large spectra.
        A reasonable choice could, e.g., be five.  

    Returns
    -------
    Broadened data : array
        The input data convolved with the Gaussian
        kernel.
    """
    # Check whether x-axis is linear
    dxs = x[1:] - x[0:-1]
    if abs(max(dxs) - min(dxs)) > np.mean(dxs) * 1e-6:
        raise(PE.PyAValError("The x-axis is not equidistant, which is required.",
                             where="broadGaussFast",
                             solution=["Use pyasl.equidistantInterpolation to produce equidistantly sampled data."]))

    if maxsig is None:
        lx = len(x)
    else:
        lx = int(((sigma * maxsig) / dxs[0]) * 2.0) + 1
    # To preserve the position of spectral lines, the broadening function
    # must be centered at N//2 - (1-N%2) = N//2 + N%2 - 1
    nx = (np.arange(lx, dtype=int) - sum(divmod(lx, 2)) + 1) * dxs[0]
    gf = fuf.GaussFit1d()
    gf["A"] = 1.0
    gf["sig"] = sigma
    e = gf.evaluate(nx)
    # This step ensured that the
    e /= np.sum(e)

    if edgeHandling == "firstlast":
        nf = len(y)
        y = np.concatenate((np.ones(nf) * y[0], y, np.ones(nf) * y[-1]))
        result = np.convolve(y, e, mode="same")[nf:-nf]
    elif edgeHandling is None:
        result = np.convolve(y, e, mode="same")
    else:
        raise(PE.PyAValError("Invalid value for `edgeHandling`: " + str(edgeHandling),
                             where="broadGaussFast",
                             solution="Choose either 'firstlast' or None"))
    return result


def instrBroadGaussFast(wvl, flux, resolution, edgeHandling=None, fullout=False, maxsig=None, equid=False):
    """
    Apply Gaussian instrumental broadening. 

    This function broadens a spectrum assuming a Gaussian
    kernel. The width of the kernel is determined by the
    resolution. In particular, the function will determine
    the mean wavelength and set the Full Width at Half
    Maximum (FWHM) of the Gaussian to
    (mean wavelength)/resolution. 

    Parameters
    ----------
    wvl : array
        The wavelength
    flux : array
        The spectrum
    resolution : int
        The spectral resolution.
    edgeHandling : string, {None, "firstlast"}, optional
        Determines the way edges will be handled. If None,
        nothing will be done about it. If set to "firstlast",
        the spectrum will be extended by using the first and
        last value at the start or end. Note that this is
        not necessarily appropriate. The default is None.
    fullout : boolean, optional
        If True, also the FWHM of the Gaussian will be returned.
    maxsig : float, optional
        The extent of the broadening kernel in terms of
        standard deviations. By default, the Gaussian broadening
        kernel will be extended over the entire given spectrum,
        which can cause slow evaluation in the case of large spectra.
        A reasonable choice could, e.g., be five.
    equid : boolean, optional
        If True, linear interpolation will be used to obtain an
        intermediate spectrum with uniformly sampled wavelength grid
        to carry out the broadening.

    Returns
    -------
    Broadened spectrum : array
        The input spectrum convolved with a Gaussian
        kernel.
    FWHM : float, optional
        The Full Width at Half Maximum (FWHM) of the
        used Gaussian kernel.
    """
    # Check whether wvl axis is linear
    dwls = wvl[1:] - wvl[0:-1]
    # True if wavelength axis is not equidistant
    dwlned = abs(max(dwls) - min(dwls)) > np.mean(dwls) * 1e-6
    if dwlned and (not equid):
        raise(PE.PyAValError("The wavelength axis is not equidistant, which is required.",
                             where="instrBroadGaussFast", \
                             solution=["Define equidistant wavelength axis", "Set 'equid' to True",
                                       "Use pyasl.equidistantInterpolation to produce equidistantly sampled data."]))
    elif dwlned and equid:
        # Define an equidistant wavelength grid
        wvlorig = wvl
        wvl = np.linspace(wvl[0], wvl[-1], 2*len(wvl))
        flux = sci.interp1d(wvlorig, flux)(wvl)
    elif (not dwlned) and equid:
        # Wavelength axis is equidistant and equid is also True (not required) 
        equid = False
    elif (not dwlned) and (not equid):
        # Wavelength axis is equidistant and no intermediate wvl axis requested
        pass
        
    meanWvl = np.mean(wvl)
    fwhm = 1.0 / float(resolution) * meanWvl
    sigma = fwhm / (2.0 * np.sqrt(2. * np.log(2.)))

    result = broadGaussFast(
        wvl, flux, sigma, edgeHandling=edgeHandling, maxsig=maxsig)

    if equid:
        # Recover wvl grid
        result = sci.interp1d(wvl, result)(wvlorig)

    if not fullout:
        return result
    else:
        return (result, fwhm)


def thermalBroadeningWidth(lam0, T, m=None, fwhm=True):
    """
    Calculate the width for thermal broadening.

    Thermal motion of particles causes a Doppler
    broadening of the line profile. The resulting
    line profile is Gaussian with FWHM given by

    .. math::

        fwhm = \\lambda_0\\sqrt{\\frac{8k_B T\\ln(2)}{m c^2}}

    See, e.g.,
    http://hyperphysics.phy-astr.gsu.edu/hbase/atomic/broaden.html

    Parameters
    ----------
    lam0 : float
        Wavelength at which to calculate the width.
    T : float
        Temperature [K].
    m : float, optional
        Mass of the particles [kg]. If not specified,
        the proton mass is assumed.
    fwhm : boolean, optional
        If True (default), the FWHM of the Gaussian
        broadening kernel will be returned. Otherwise,
        the standard deviation is returned.

    Returns
    -------
    Width : float
        The width of the Gaussian broadening kernel.
        By default, the FWHM is returned. To obtain
        the standard deviation, set the `fwhm` flag to
        False.
    """
    from PyAstronomy import constants as PC
    pc = PC.PyAConstants()
    pc.setSystem("SI")
    if m is None:
        # Use proton mass if not specified otherwise
        m = pc.mp
    result = lam0 * np.sqrt(8.0 * pc.k * T * np.log(2.0) / (m * (pc.c)**2))
    if not fwhm:
        return result / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    else:
        return result
    
    
    
def tempFromthermalBroadeningWidth(lam0, width, m=None, awidth=0.0):
    """
    Calculate the temperature required to obtain thermal broadening width.

    Thermal motion of particles causes Doppler
    broadening of the line profile. The resulting
    line profile is Gaussian with FWHM given by

    .. math::

        fwhm = \\lambda_0\\sqrt{\\frac{8k_B T\\ln(2)}{m c^2}}

    See, e.g.,
    http://hyperphysics.phy-astr.gsu.edu/hbase/atomic/broaden.html

    Here, the relation is reversed, so that the temperature, T, can be
    obtained from the width. The relation is only strictly valid in the case
    of Gaussian lines.

    .. note:: Units of lam0, width, and awidth have to be consistent.

    Parameters
    ----------
    lam0 : float
        Wavelength at which to width (FWHM) was measured.
    width : float
        Measured width of the feature (FWHM) 
    m : float, optional
        Mass of the particles [kg]. If not specified,
        the proton mass is assumed.
    awidth : float
        Additional width contributed by other effects such as
        instrumental resolution. The thermal width, ft, is obtained by
        quadratic subtraction, i.e., ft**2 = width**2 - awidth**2.
        Default is zero.

    Returns
    -------
    Temperature : float
        The temperature required to achieve the specified broadening
        width.
    Thermal width : float
        The thermal width used in the calculations (may be modified
        if `awidth` has been specified.
    """
    from PyAstronomy import constants as PC
    pc = PC.PyAConstants()
    pc.setSystem("SI")
    if m is None:
        # Use proton mass if not specified otherwise
        m = pc.mp
        
    if awidth > width:
        raise(PE.PyAValError("The additional width (awidth) is larger than the specified width of the feature.", \
                             where="tempFromthermalBroadeningWidth", \
                             solution="Check width specification.")) 
        
    # Take into account additional width
    width = np.sqrt(width**2 - awidth**2)
        
    result = width**2 / lam0**2 * (m * (pc.c)**2) / (8.0 * pc.k * np.log(2.0))
    return result
    

