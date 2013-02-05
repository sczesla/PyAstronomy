import numpy
from PyAstronomy.pyaC import pyaErrors as PE

def isInTransit(time, T0, period, halfDuration, boolOutput=False):
  """
    Check whether time is inclosed by transit interval.
    
    This function uses the given ephemerides (T0, period, and
    halfDuration) to check whether the time point(s)
    specified by `time` are within a transit window or not.
    The edges of the window are counted as transit times.
    
    Parameters
    ----------
    time : float or array like
        The time(s) which are to be checked.
    T0 : float
        The time reference point (center of transit).
    period : float
        The orbital period.
    halfDuration : float
        The half-duration of the event.
        Must have same units as `time`.
    boolOutput : boolean, optional
        If set True and `time` is an array, the function will
        return a bool array holding True for time points in-
        and False for time points out-of-transit.
    
    Returns
    -------
    inTransit : boolean or array of int
        If `time` was a float, the return value is a boolean,
        which is True if the give time falls into a transit
        interval and False otherwise.
        If `time` was given as an array, the return value is
        an array holding the indices of those time points,
        which fall into a transit window. The `boolOutput`
        option may be used to obtain a boolean array holding
        True for in-transit points. 
    
  """
  if halfDuration > period/2.:
    raise(PE.PyAValError("The half-duration is longer than half the period. This cannot be true.", \
                         where="isInTransit"))
  if (period <= 0.0) or (halfDuration <= 0.0):
    raise(PE.PyAValError("Both period and half-duration must be larger 0.", \
                         where="isInTransit"))
  absPhase = numpy.abs((numpy.array(time) - T0)/period)
  absPhase -= numpy.floor(absPhase)
  dPhase = halfDuration/period
  isIn = numpy.logical_or(absPhase <= dPhase, absPhase >= (1.-dPhase))
  indi = numpy.where(isIn)[0]
  if isinstance(time, float):
    return (len(indi) == 1)
  if boolOutput:
    return isIn
  return indi
  