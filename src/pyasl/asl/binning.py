# -*- coding: utf-8 -*-
from PyAstronomy.pyaC import pyaErrors as PE
import numpy as np
import copy

def binningx0dt(x, y, yerr=None, x0=None, dt=None, nbins=None, reduceBy=None, removeEmpty=True, \
                removeNoError=False, useBinCenter=True, useMeanX=False, nanHandling=None):
  """
    A simple binning algorithm.
    
    This algorithm uses a fixed bin-width to produce a binned
    data set. Either the bin-width, `dt`, or the number of bins,
    `nbins`, must be specified. The number of output bins may
    also depend on other flags such as, for example, `removeNoError`.
    
    If no errors are specified via `yerr`, the errors for the binned
    data are estimated as the standard deviation of the input data
    points divided by the square root of their number. If `yerr` has
    been specified, error propagation is used to determine the error.
    
    The behavior of the x-axis can be controlled via the
    `useBinCenter` flag.
    
    Values which cannot be determined will be indicated by NaN.
    Various flags can be used to remove such bins from the binned
    data set.
    
    Parameters
    ----------
    x, y : array
        The x and y data values.
    yerr : array, optional
        Errors on the data values.
    x0 : float, optional
        Starting time of first bin.
        Default is lowest given x value.
    dt : float, optional
        Width of a bin (either `dt`, `nbins` or `reduceBy` must
        be given).
    nbins : int, optional
        Number of bins to use (either `dt`, `nbins` or `reduceBy` must
        be given). Note that this specifies the number of bins into which
        the range from `x0` to the last data point is subdivided.
    reduceBy : int, optional
        Reduce the number of elements in the array by the given factor 
        (either `dt`, `nbins` or `reduceBy` must be given). Note that
        in this case, `x0` is set to the first (minimum x-value) and
        the number of bins, n, is calculated according to the
        prescription: :math:`n = int(round(len(x)/reduceBy))`
    removeEmpty : boolean, optional
        If True (default), bins with no data points will be
        removed from the result.
    removeNoError : boolean, optional
        If True, bins for which no error can be determined
        will be removed from the result. Default is False.
    useBinCenter : boolean, optional
        If True (default), the time axis will refer to the
        center of the bins. Otherwise the numbers refer to
        the start of the bins.
    useMeanX : boolean, optional
        If True, the binned x-values refer to the mean x-value
        of all points in that bin.
        Therefore, the new time axis does not have to be equidistant.
    nanHandling : None, "ignore", float, (optional)
        Controls how NaNs in the data are handled.
          - None: By default (None),
            nothing is done and NaNs are treated as if they were valid
            input data, so that they are carried over into the binned data.
            This means that output bins containing NaN(s) will also end up
            as NaN(s). If 'ignore'
          - 'ignore': In this case, NaNs contained in the input data are
            removed from the data prior binning. Note however, that `x0`,
            unless specified explicitly, will still refer to the first data
            point, whether or not this holds a NaN value. 
          - float: If a float is given, input data values containing NaNs
            are replaced by the given float before binning. Note that no error on
            the data (yerr) can be considered in this case, to avoid 
            erronous treatment of un- or misspecified error values.
    
    Returns
    -------
    Binned data set : array
        An array with four columns: 1) The new x-axis,
        2) The binned data (the mean value of the data points located
        in the individual bins), 3) Error of binned data, 4) The
        number of input data points used to create the bin.
        For instance, the new x-values can be accessed
        using result[::,0].
    dt : float
        The width of the bins.
  """
  if len(x) != len(y):
    raise(PE.PyAValError("x and y need to have the same length."))
  if ((not dt is None) + (not nbins is None) + (not reduceBy is None)) != 1:
    raise(PE.PyAParameterConflict("Specify one of `dt`, `nbins`, or `reduceBy`."))
  if ((not x0 is None) + (not reduceBy is None)) != 1:
    raise(PE.PyAParameterConflict("Specify either `x0` or `reduceBy`."))
  if x0 is None:
    # Use first time as starting point
    x0 = np.min(x)
  if x0 > np.max(x):
    raise(PE.PyAValError("The starting point, `x0`, is larger than the end time of the data.", \
                         solution="Use a smaller value."))

  # Use arrays in calculation. Only copy if conversion to numpy array
  # is required
  xl, yl, yerrl = False, False, False
  if not isinstance(x, np.ndarray):
    x = np.array(x)
    xl = True
  if not isinstance(y, np.ndarray):
    y = np.array(y)
    yl = True
  if not yerr is None:
    if not isinstance(yerr, np.ndarray):
      yerr = np.array(yerr)
      yerrl = True 

  # nanHandling
  if nanHandling is not None:
    # As some manipulation of the data may be required, generate
    # local copies, unless this has already happened.
    if not xl:
      x = copy.copy(x)
    if not yl:
      y = copy.copy(y)
    if not yerrl:
      yerr = copy.copy(yerr)
    
    if nanHandling == "ignore": # Remove bins containing NaN from the input data
      gi = np.isfinite(y)
      x, y = x[gi], y[gi]
      if not yerr is None:
        # Check if bins containing NaN from the input data
        yerr = yerr[gi]
    else: 
      try:
        nanValue = float(nanHandling)
      except ValueError as ve:
        raise(PE.PyAValError("Invalid value given for 'nanHandling'. Error was: " + str(ve), \
                             solution="Use None, 'ignore', or float number."))

      gi = np.isnan(y)
      y[gi] = nanValue
      if yerr is not None:
        raise(PE.PyAValError("yerr is not permitted for option 'nanHandling=float'.", \
                             solution="Remove yerr=... from the call."))
          
  # Calculate the new number of array elements.                         
  if reduceBy is not None:
    nbins = int(round(len(x)/float(reduceBy))) 
    if nbins == 0: nbins=1 # Prevent empty return arrays
  if nbins is not None:
    # Use a specified number of bins.
    # Calculate bin length
    dt = (np.max(x) - x0)/float(nbins)
  # Start calculation
  # In which bin do the individual data points belong?
  inWhichBin = np.floor(((x-x0)/dt)).astype(np.int)
  # Lonely last bin correction
  # Brings the last data point into the last valid bin
  # instead of creating a new bin with that data point\
  # at its very beginning
  if nbins is not None:
    inWhichBin[np.where(inWhichBin == nbins)[0]] -= 1
  # Get the number of bins (start at x0 even if the
  # first bins do not contain any data points)
  nbins = np.max(inWhichBin) + 1
  # Bins with data
  bwd = np.unique(inWhichBin)
  # Sort data into the bins
  # Create output array (time, flux, error, data-point-counter)
  result = np.empty( (nbins, 4) )
  result[:] = np.NAN
  # Assign time axis (beginning of bins)
  result[::,0] = x0 + np.arange(nbins) * dt
  if useBinCenter:
    # Use the center of the bin for timing
    result[::,0] += (0.5 * dt)
  # Set data point counter (points/bin) to zero
  result[::,3] = 0
  for b in bwd:
    indi = np.where(inWhichBin == b)[0]
    result[b, 3] = len(indi)
    result[b, 1] = np.mean(y[indi])
    if useMeanX:
      # Overwrite the time axis using the mean x-value
      result[b, 0] = np.mean(x[indi])
    if yerr is None:
      # No errors on data points are given
      if len(indi) > 1:
        result[b, 2] = np.std(y[indi]) / np.sqrt(result[b, 3])
      else:
        # No error if there is only a single point in the bin
        result[b, 2] = np.NAN
    else:
      # There are errors on the data points
      # Use error propagation
      result[b, 2] = np.sqrt(np.sum(yerr[indi]**2)) / result[b, 3]
  
  if removeEmpty:
    # Remove bins without data points in it
    indi = np.where(np.invert(np.isnan(result[::,1])))[0]
    result = result[indi,::]
  if removeNoError:
    # Remove bins for which no error can be given
    indi = np.where(np.invert(np.isnan(result[::,2])))[0]
    result = result[indi,::]
  return result, dt