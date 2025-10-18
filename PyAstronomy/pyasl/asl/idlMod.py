import numpy as np

def idlMod(a, b):
  """
    Emulate 'modulo' behavior of IDL.
    
    Parameters
    ----------
    a : float or array
        Numerator
    b : float
        Denominator
    
    Returns
    -------
    IDL modulo : float or array
        The result of IDL modulo operation.
  """
  if isinstance(a, np.ndarray):
    s = np.sign(a)
    m = np.mod(a, b)
    m[(s < 0)] -= b
  else:
    m = a % b
    if a < 0: m -= b
  return m