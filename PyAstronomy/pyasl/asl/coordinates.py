from PyAstronomy.pyaC import pyaErrors as PE
import re


def hmsToDeg(h, m, s):
    """
    Convert hour-minute-second specification into degrees.
    
    Parameters
    ----------
    h : float
        Hours (0-24)
    m : float
        Minutes (time, 0-60)
    s : float
        Seconds (time, 0-60)
    
    Returns
    -------
    Angle : float
        The corresponding angle in degrees.
    """
    if (h < 0.0) or (h >= 24.0):
        raise(PE.PyAValError("Hour (" + str(h) + ") out of range (0 <= h < 24).",
                             solution="Specify a value between 0 and 24."))
    if (m < 0.0) or (m >= 60.0):
        raise(PE.PyAValError("Minute (" + str(m) + ") out of range (0 <= m < 60)",
                             solution="Specify a value between 0 and 60."))
    if (s < 0.0) or (s >= 60.0):
        raise(PE.PyAValError("Second (" + str(s) + ") out of range (0 <= m < 60)",
                             solution="Specify a value between 0 and 60."))
    result = h * 15.0 + m * 0.25 + s * (0.25 / 60.0)
    return result


def degToHMS(d):
    """
    Convert degrees into time units (hours-minutes-seconds)
    
    Parameters
    ----------
    d : float
        Degrees (0-360)
    
    Returns
    -------
    h, m, s: float
        Hours, minutes, and seconds
    """
    if (d < 0.0) or (d >= 360.0):
        raise(PE.PyAValError("The argument should be between 0 and <360 degrees."))
    h = int(d / 15.0)
    d = d - h * 15.0
    m = int(d / 0.25)
    d = d - m * 0.25
    s = d / (0.25 / 60.0)
    return h, m, s


def degToDMS(g):
    """
    Convert degrees into arc units (degrees, (arc)minutes, (arc)seconds)
    
    Parameters
    ----------
    g : float
        Value in degrees.
    
    Returns
    -------
    d, m, s, sign : float, int
        Degrees, (arc)minutes, and (arc)seconds. Note that only the
        degree number is signed. The sign (+1/-1) is also returned to
        yield a complete result if the value of degree (d) is zero. 
    """
    sign = 1
    if g < 0.0:
        sign = -1
    g = abs(g)
    d = int(g) * sign
    g = g - int(g)
    m = int(g * 60.0)
    g = g - m / 60.0
    s = g * 3600.0
    return d, m, s, sign


def dmsToDeg(d, m, s, esign=0):
    """
    Convert degree-arcminute-arcsecond specification into degrees.
    
    Parameters
    ----------
    d : float
        Degrees.
    m : float
        Arcminutes (0-60)
    s : float
        Arcseconds (0-60)
    esign : int, optional, {-1,0,1}
        Explicit sign with -1 representing negative sign, +1 representing
        positive sign, and 0 indicating no explicit sign specification.
        The explicit sign is necessary if negative southern coordinates are
        specified but d is 0 and, thus, cannot carry the sign. 
    
    Returns
    -------
    Angle : float
        The corresponding angle in degrees.
    """
    if (m < 0.0) or (m >= 60.0):
        raise(PE.PyAValError("Minute (" + str(m) + ") out of range (0 <= m < 60)",
                             solution="Specify a value between 0 and 60."))
    if (s < 0.0) or (s >= 60.0):
        raise(PE.PyAValError("Second (" + str(s) + ") out of range (0 <= m < 60)",
                             solution="Specify a value between 0 and 60."))
    if not esign in (-1, 0, +1):
        raise(PE.PyAValError("Invalid value for 'esign' (value given is " + str(esign) + ").",
                             solution="Use +1, 0, or -1.",
                             where="dmsToDeg"))
    if (d == 0.0) and (esign == 0):
        raise(PE.PyAValError("When a value of zero is specified for degrees (d), the sign has to be given explicitly via the `esign` flag.",
                             where="dmsToDeg",
                             solution="Specify esign=+/-1 according to sign (i.e., northern or southern coordinates).",
                             addInfo="Safety device added to prevent common conversion errors from coordinates like '-00 12 08'."))

    sign = 1.0
    if d < 0.0:
        sign = -1.0

    if (d != 0.0) and (esign != 0):
        # Check consistency of sign and esign
        if int(sign) != esign:
            raise(PE.PyAValError("The sign of 'd' and the specification of 'esign' must not contradict.",
                                 solution="Adjust 'esign' specification.",
                                 where="dmsToDeg"))
    elif d == 0.0:
        sign = esign

    result = d + sign * (m / 60.0) + sign * (s / 3600.0)
    return result


def coordsSexaToDeg(c, fullOut=False):
    """
    Convert sexagesimal coordinate string into degrees.
    
    Parameters
    ----------
    c : string
        The coordinate string. Valid formats are, e.g.,
        "00 05 08.83239 +67 50 24.0135" or "00:05:08.83239 -67:50:24.0135".
        Spaces or colons are allowed as separators for the individual
        components of the coordinates.
    fullOut : boolean, optional
        If True, two additional tuples holding the individual components of
        the right ascension and declination specification will be returned.
        The default is False.
    
    Returns
    -------
    ra, dec : float
        Right ascension and declination in degrees.
    hms, dms : tuples of three floats
        If `fullOut` is True, two tuples of three numbers each will be
        returned, which hold the individual constituents making up the
        right ascension (hms) and declination (dms) specifiers in sexagesimal
        representation.
    """
    r = re.match(
        "^\s*(\d+)([\s:])(\d+)([\s:])(\d+(\.\d*)?)\s+(([+-])(\d+))([\s:])(\d+)([\s:])(\d+(\.\d*)?)\s*$", c)
    if r is None:
        raise(PE.PyAValError("Could not decompose coordinate string: \"" + str(c) + "\"",
                             solution="Use, e.g., 00 05 08.83239 +67 50 24.0135 or 00:05:08.83239 +67:50:24.0135"))
    # Check separator consistency
    for n in [4, 10, 12]:
        if r.group(2) != r.group(n):
            raise(PE.PyAValError("Separators (space and colon) mixed in coordinate definition.",
                                 solution="Use consistent separator."))

    ra = (float(r.group(1)), float(r.group(3)), float(r.group(5)))
    esign = {'+': +1, '-': -1}[r.group(8)]
    de = (float(r.group(7)), float(r.group(11)), float(r.group(13)), esign)
    result = (hmsToDeg(*ra), dmsToDeg(*de))
    if not fullOut:
        return result
    return result[0], result[1], ra, de


def coordsDegToSexa(ra, dec, asString=True, fmt=("%02d %02d %06.3f  ", "%s%02d %02d %06.3f")):
    """
    Convert right ascension and declination from degrees into sexagesimal representation.
    
    Parameters
    ----------
    ra : float
        Right ascension in degrees.
    dec : float
        Declination in degrees.
    asString : boolean, optional
        If True (default), the result will be a string formatted
        according to the rules specified by `fmt`.
    fmt : tuple of strings, optional
        The output format used to create the output string (first ra, second dec). Only used
        if `asString` is True (default).
    
    Returns
    -------
    Coordinates : string or tuple
        If `asString` is True (default), a string holding the coordinates
        is returned, which is formatted according to the rules specified by
        the `fmt` parameter. If False, a tuple of two tuples is returned,
        of which the first holds three numbers representing the right ascension
        (hms) and the second three numbers representing the declination (dms, sign).    
    """
    if (ra < 0.0) or (ra >= 360.):
        raise(PE.PyAValError("Right ascension should be between 0 and 360 degrees."))
    if (dec < -90.0) or (dec > 90.0):
        raise(PE.PyAValError("Declination should be between -90 and +90 degrees."))
    rat = degToHMS(ra)
    dect = degToDMS(dec)

    if asString:
        result = (fmt[0] % rat) + (fmt[1] %
                                   ({-1: '-', +1: '+'}[dect[3]], abs(dect[0]), dect[1],  dect[2]))
        return result
    return rat, dect
