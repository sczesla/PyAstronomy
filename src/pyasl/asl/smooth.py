# -*- coding: utf-8 -*-
from PyAstronomy.pyaC import pyaErrors as PE
import numpy as np


def smooth(x, windowLen, window='flat'):
  """
    Smooth data using a window function.
    
    This method is based on the convolution of a window function with the signal.
    The window function is normalized so that the sum of its entries amounts to
    one. The signal is prepared by adding reflected copies of the signal 
    (with the window size) to both ends of the input array, so that the output
    array can have the same length as the input. Consequently the smoothing at
    the edges is actually based on extrapolation.
    
    .. note:: This algorithm was adopted from the scipy cookbook
              (http://www.scipy.org/Cookbook/SignalSmooth). The copyright
              of the original algorithm belongs to the authors of that
              cookbook algorithm.
    
    Parameters
    ----------
    x : array
        The input signal 
    windowLen : int
        The dimension of the smoothing window. It must be an odd integer.
    window : string, {'flat', 'hanning', 'hamming', 'bartlett', 'blackman'}
        The window function to be used. A flat window will
        produce a moving average smoothing.

    Returns
    -------
    Smoothed signal : array
        The smoothed signal. Same length as input array.
  """

  if x.ndim != 1:
    raise(PE.PyAValError("Only one dimensional arrays can be smoothed. Dimension of " + \
          "given array is " + str(x.ndim),
          solution="Change dimension of input."))

  if x.size < windowLen:
    raise(PE.PyAValError("Input vector needs to be bigger than window size.", \
          solution="Check the length of the input array and window size."))

  if windowLen < 3:
    PE.warn(PE.PyAValError("Length of window is smaller then 3. No smoothing is done.", \
        solution="Check window size."))
    return x
  
  if windowLen % 2 != 1:
      raise(PE.PyAValError("Parameter `windowLen` should be an odd integer"))
  
  # Extend input array at the edges to have the same
  # length for the output array. Insert a mirrored version
  # of the first part of the data array in front of the
  # first data point; apply the same scheme to the end of the
  # data array. 
  s = np.r_[x[windowLen-1:0:-1], x, x[-1:-windowLen:-1]]
  
  if window == 'flat':
    # This is a moving average
    w = np.ones(windowLen, 'd')
  elif window == "hanning":
    w = np.hanning(windowLen)
  elif window == "hamming":
    w = np.hamming(windowLen)
  elif window == "bartlett":
    w = np.bartlett(windowLen)
  elif window == "blackman":
    w = np.blackman(windowLen)
  else:
    raise(PE.PyAValError("Current `window` parameter (" + str(window) + ") is not supported. " + \
                         "Must be one of: 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'", \
          solution="Choose one of the above window types."))

  y = np.convolve(w/w.sum(), s, mode='valid')
  return y[(windowLen/2):-(windowLen/2)]
