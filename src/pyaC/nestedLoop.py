from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo

class Counter:
  """
    A cyclic counter.
    
    On each call to `increment`, the counter
    will be increased. The `counter` attribute
    counts from minimum to (maximum-1). The
    `status` attribute counts from 0 to
    (maximum-minimum-1).
    
    Parameters
    ----------
    minimax : tuple of two int
        The minimum and maximum of the counter.
  """

  def __init__(self, minmax):
    self.mini, self.maxi = minmax
    self.maxmimin = self.maxi - self.mini
    self.status = 0
    self.counter = self.status + self.mini

  def increment(self):
    """
      Increment the counter.
      
      Returns
      -------
      Overflow flag : int
          Returns 1 if the counter reached the maximum and
          jumped back to the minimum and zero otherwise.
    """
    self.status += 1
    q, self.status = divmod(self.status, self.maxmimin)
    self.counter = self.mini + self.status
    return q


class NestedLoop:
  """
    Implements an iteration over a nested loop.
    
    Iterates over nested loops. First, increases the
    first counter, then increments the second and so
    on.
    
    Parameters
    ----------
    limits : list of int
        The upper limits for the loops.
    lowerLimits : list of int, optional
        The lower limits of the loops. The default
        is zero.
  """

  def __init__(self, limits, lowerLimits=None):
    self.limits = limits[:]
    if lowerLimits is None:
      self.lowerLimits = [0]*len(limits)
    else:
      self.lowerLimits = lowerLimits[:]
      if len(self.lowerLimits) != len(limits):
        raise(PE.PyAValError("You need to specify as many lower as upper limits.", \
                             solution="Make `limits` and `lowerLimits` the same length."))
    for i in smo.range(len(self.limits)):
      if self.limits[i] <= self.lowerLimits[i]:
        raise(PE.PyAValError("Upper limit needs to be larger than lower limit", \
                             solution="Adapt ranges!"))
    # Number of loops
    self._nl = len(limits)
    # Total number of points 
    self._np = 1
    for i in smo.range(self._nl):
      self._np *= (self.limits[i] - self.lowerLimits[i])
    # List of counters
    self.counters = list(map(Counter, list(zip(self.lowerLimits, self.limits))))

  def __iter__(self):
    """
      Iterator for the loop.
      
      Returns
      -------
      Loop index : tuple
          A tuple with an int for every loop index.
    """
    for j in smo.range(self._np):
      yield tuple([x.counter for x in self.counters])
      for i in smo.range(self._nl):
        x = self.counters[i].increment()
        if x == 0: break
