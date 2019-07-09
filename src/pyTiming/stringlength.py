import numpy as np
from PyAstronomy import pyasl
from PyAstronomy.pyaC import pyaErrors as PE 

def stringlength_norm(m, norm):
    """
    Normalize string length data set.
    
    Parameters
    ----------
    m : array
        Data array
    norm : string, {default, no}
        If 'default' (default), the data points (mi) will be renormalized according
        to Eq. 3 in Dworetsky 1983, i.e., using mnew = (mi - min(m)) / (2*(max(m) - min(m))).
        If 'no' is specified, the data will not be changed.
    
    Returns
    -------
    ms : array
        The normalized data
    """
    if norm == "default":
        mmin = np.min(m)
        mmax = np.max(m)
        # Dworetsky 1983, Eq. 3
        ms = (m-mmin) / (2 * (mmax-mmin)) - 0.25
    elif norm == "no":
        # Do nothing wrt normalization
        ms = m
    else:
        raise(PE.PyAValError("Unknown value for 'norm': " + str(norm), \
                             where="stringlength_norm", \
                             solution="Use 'default' or 'no'"))
    return ms


def stringlength_pm(p, m, norm="default", closed=True):
    """
    Compute the string length for phased data set.
    
    Parameters
    ----------
    p : array
        The phase array (0-1). Sorted in ascending order.
    m : array
        Data array
    norm : string, {default, no}
        If 'default' (default), the data points (mi) will be renormalized according
        to Eq. 3 in Dworetsky 1983, i.e., using mnew = (mi - min(m)) / (2*(max(m) - min(m))).
        If 'no' is specified, the data will not be changed.
    closed : boolean, optional
        If True (default), first and last point on the phase axis
        will be connected (close the loop).
    
    Returns
    -------
    sl : float
        The string length
    """
    
    if np.any(np.diff(p) < 0):
        raise(PE.PyAValError("Phase array (p) is not sorted in ascending order.", \
                             where="stringlength_pm"))
    
    ms = stringlength_norm(m, norm)

    # String length
    sl = np.sum(np.sqrt( np.diff(p)**2 + np.diff(ms)**2 ))
    if closed:
        # Close the loop ...
        sl += np.sqrt( (p[0] - p[-1] + 1)**2 + (ms[0] - ms[-1])**2 )
    
    return sl


def stringlength_dat(x, m, tps, norm='default', isFreq=False, closed=True):
    """
    Compute string length for data set.
    
    Parameters
    ----------
    x, m : arrays
        x and y coordinates of data points.
    tps : tuple or array
        The trial periods (or frequencies): Either a three-tuple specifying
        (pmin, pmax, nperiods) used by numpy's linspace or an array of trial periods.
        If isFreq is True, tps is interpreted as frequencies.
    isFreq : boolean, optional
        If True, the input tps will be assumed to refer to frequencies instead of periods.
    norm : string, {default, no}
        If 'default' (default), the data points (mi) will be renormalized according
        to Eq. 3 in Dworetsky 1983, i.e., using mnew = (mi - min(m)) / (2*(max(m) - min(m))).
        If 'no' is specified, the data will not be changed.
    closed : boolean, optional
        If True (default), first and last point on the phase axis
        will be connected (close the loop).
    
    Returns
    -------
    Trial periods/frequencies : array
        The tested periods (or frequencies if isFreq is True).
    String length : array
        Associated string lengths
    """
    
    if isinstance(tps, tuple) and len(tps) == 3:
        pers = np.linspace(*tps)
    elif isinstance(tps, np.ndarray):
        pers = tps
    else:
        raise(PE.PyAValError("Unrecognized specification of trial periods or frequencies (tps).", \
                             where="stringlength_dat"))
    
    if isFreq:
        # The input is frequencies. Still need the periods.
        pers = 1.0/pers
    
    ms = stringlength_norm(m, norm)
    
    sls = np.zeros(len(pers))
    for i, p in enumerate(pers):
        phase = pyasl.foldAt(x, p)
        indi = np.argsort(phase)
        sls[i] = stringlength_pm(phase[indi], ms[indi], norm='no', closed=closed)
    
    if isFreq:
        return 1.0/pers, sls
    
    return pers, sls

