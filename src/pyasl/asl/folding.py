import numpy

def foldAt(time, period, T0=0.0):
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
  
    Returns
    -------
    Phases : array
        The (unsorted) phase array pertaining to the input time axis.
  """
  phase = (time - T0)/period - numpy.floor( (time - T0)/period )
  return phase