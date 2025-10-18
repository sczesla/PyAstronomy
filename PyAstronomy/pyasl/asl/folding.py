import numpy as np


def foldAt(time, period, T0=0.0, sortphase=False, centralzero=False, getEpoch=False):
    """
    Fold time series with a particular period.

    Calculate the phase, P, from time, period, and reference point, T0,
    according to 

    .. math::
    
        P = (time - T0)/period - [(time-T0)/period].
    
    Here, square brackets indicate Gaussian brackets (i.e., the floor function), and
    the phase is a number between 0 and 1 by definition (and not between 0 and 2pi).
    
    Optionally, also the epoch, E, can be returned, which is an integer corresponding
    to the second term in the above equation. For any point of the series, therefore,
    the following relation applies
    
    .. math::
    
        time = T0 + (E+P) * period .

    Of course the series to be folded does not necessarily have to be a time
    series although this particular example guides the naming convention here.
    
    Parameters
    ----------
    time : array
        The time stamps.
    period : float
        The period to fold with (same units as time stamps).
    T0 : float
        Time reference point. The point T0 as well as all points T0+n*period with
        integer n are mapped to phase zero. Default is 0.0.
    sortphase : boolean, optional
        If True, return values will be sorted w.r.t. phase
    centralzero: boolean, optional
        If True, phase will be between -0.5 and +0.5 instead of 0 to 1. 
    getEpoch : boolean, optional
        If True, an array holding the epoch for every point
        in time will be returned; the default is False.
        Note that the first epoch, corresponding to times between T0 and T0+per, is 0.

    Returns
    -------
    Phases : array
        The (unsorted) phase array pertaining to the input time axis.
        Sorted if `sortphase` is set True.
    Epoch : array, optional
        An array holding the epoch for every given point in time.
        The counting starts at zero. Only returned if `getEpoch`
        is True.
    """
    epoch = np.floor((time - T0) / period)
    phase = (time - T0) / period - epoch
    
    if centralzero:
        indi = np.where(phase >= 0.5)[0]
        phase[indi] -= 1.0
    
    if sortphase:
        indi = np.argsort(phase)
        phase, epoch = phase[indi], epoch[indi]
    
    if getEpoch:
        return phase, epoch
    return phase
