import numpy

def foldAt(time, period, T0=0.0, getEpoch=False):
  """
    Fold time series at a particular period.
    
    In this function, the phase is a number between 0 and 1
    (and not between 0 and 2pi).
    
    Parameters
    ----------
    time : array
        The time stamps.
    period : float
        The period to fold at (same time units is period.)
    T0 : float
        Time reference point
    getEpoch : boolean, optional
        If True, an array holding the epoch for every point
        in time will be returned; the default is False.
        Note that the first epoch is 0.
  
    Returns
    -------
    Phases : array
        The (unsorted) phase array pertaining to the input time axis.
    Epoch : array, optional
        An array holding the epoch for every given point in time.
        The counting starts at zero. Only returned if `getEpoch`
        is True.
  """
  epoch = numpy.floor( (time - T0)/period )
  phase = (time - T0)/period - epoch
  if getEpoch:
    return phase, epoch
  return phase