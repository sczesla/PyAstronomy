from . import pyaErrors as PE
import numpy as np

def invertIndexSelection(a, indi):
  """
    Invert index selection in one-dimensional array.
    
    Say, e.g., a numpy.where operation produced an array
    of indices (`indi`), which you needed for one reason,
    but now, for another reason, you need all elements of
    the array, which were not selected by the operation.
    This is the situation handled by this function. 
    
    Parameters
    ----------
    a : int or array
        Either the length of the array to which `indi` refers
        to or the array (one-dimensional) itself.
    indi : array
        An array of indices selected from `a`.
    
    Returns
    -------
    Inverse selection : array
        An array containing the indices of all array elements
        not referenced by `indi`. 
  """
  if isinstance(a, int):
    n = a
  elif isinstance(a, np.ndarray):
    if len(a.shape) != 1:
      raise(PE.PyAValError("Found first parameter to be a numpy array, but it is not of dimension one!", \
            solution=["Specify a one-dimensional array", "Specify merely the length of the array as an int."], \
            where="invertIndexSelection"))
    n = len(a)
  else:
    raise(PE.PyAValError("First parameter is neither a one-dimensional numpy array or an integer.", \
          solution="Use one of these types to describe the array."))
  
  x = np.ones(n, dtype=np.bool)
  x[indi] = False
  return np.arange(n)[x]
