from __future__ import print_function, division
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo

def expCorrRN(n, tau, mean=0.0, std=1.0, rnos=None, fullOut=False):
  """
    Generate exponentially correlated random numbers.
    
    This procedure implements the prescription given by
    Deserno 2002 ("How to generate exponentially correlated Gaussian
    random numbers"). The autocorrelation function of the resulting
    numbers decays with the predefined "decay time", tau.
    The correlation coefficient of the resulting numbers is given
    by exp(-1/tau).
    
    Parameters
    ----------
    n : int
        Number of numbers to be generated.
    tau : float
        Decay time
    mean : float, optional
        Mean of the numbers to be generated.
        Default is 0.0.
    std : float, optional
        Standard deviation of the generated numbers.
        Default is 1.0.
    rnos : array, optional
        Uncorrelated Gaussian random numbers with mean 0.0
        and standard deviation 1.0 to be used to generate
        correlated random numbers. If not given, Gaussian
        random numbers will be obtained using numpy.random.normal.
    fullOut : boolean, optional
        If False (default), only the correlated random numbers
        will be returned. 
    
    Returns
    -------
    Correlated RNs : array
        Correlated Gaussian random numbers.
    Uncorrelated RNs : array, optional
        The uncorrelated random numbers used to generate
        the correlated numbers (only of `fullOut` is True).
    Correlation coefficient : float
        The correlation coefficient (exp(-1/tau), only of
        `fullOut` is True).
  """
  if rnos is None:  
    # Get n uniformly distributed random numbers between 0 and 1
    g = np.random.normal(0.0, 1.0, n)
  else:
    g = rnos
  
  if len(g) != n:
    raise(PE.PyAValError("The length of `rnos` must be n.",
                         where="expCorrRN",
                         solution=["Adjust `rnos`.", \
                                   "Remove rnos argument to use internally generated Gaussian random numbers."]))
    
  # Prepare result
  result = np.zeros(n)
  result[0] = g[0]
  
  # Correlation coefficient
  f = np.exp(-1.0/tau)
  
  # Obtain correlated numbers
  somf = np.sqrt(1.0 - f**2)
  for i in smo.range(1,n):
    result[i] = f*result[i-1] + somf*g[i]
  result = mean + std*result
  
  if fullOut:
    return result, g, f
  return result