import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE

def estimateSNR(x, y, xlen, deg=1, controlPlot=False, xlenMode="dataPoints"):
  """
    Estimate Signal to Noise Ratio (SNR) in a data set.
    
    This function provides a simple algorithm to estimate
    the SNR in a data set. The algorithm subdivides the
    data set into sections with a length determined by the
    `xlen` parameter. Each individual subsection is
    fitted using a polynomial of degree `deg`. The SNR
    is then computed by assuming that the resulting
    chi square value is properly distributed according
    to the chi-square distribution.
    
    Parameters
    ----------
    x : array
        The abscissa values.
    y : array
        The ordinate values.
    xlen : int or float
        Length of the data subsection considered
        in the computation of the SNR. Whether `xlen`
        refers to the number of data points or a fixed
        subsection of the abscissa is determined by the
        `xlenMode` flag.
    deg : int, optional
        Degree of the polynomial used to fit the
        subsection of the data set. The default is
        one.
    controlPlot : boolean, optional
        If True, a control plot will be shown to
        verify the validity of the estimate.
    xlenMode : string, {"dataPoints", "excerpt", "all"}, optional
        Determines whether `xlen` refers to data points
        or a fixed length scale on the abscissa. If 'all' is specified,
        all available data will be used in the estimation.
    
    Returns
    -------
    SNR : dictionary
        Contains the SNR estimate for the individual subsections (key 'snrs')
        and the final SNR estimate (mean of all 'sub-SNRs', key 'SNR-Estimate').
  """
  # Stores the SNRs computed from the individual sections.
  snrs = []
  i = 0
  lastdof = None
  
  while True:
    if xlenMode == "dataPoints":
      # xlen is given in units of data points
      indi = range(i*xlen, (i+1)*xlen)
      i += 1
      if i > len(x)/xlen: break
    elif xlenMode == "excerpt":
      # xlen is given in units of the abscissa values
      indi = np.where(np.logical_and(x >= i*float(xlen), x < (i+1)*float(xlen)))[0]
      i += 1
      if i*float(xlen) > max(x): break
    elif xlenMode == "all":
      if i > 1:
        break
      indi = range(len(x))
      i += 1
    else:
      raise(PE.PyAValError("Unknown xlenMode (" + str(xlenMode) + ").", \
                           solution = "Use either 'dataPoints' or 'excerpt'."))
    # Check whether indi is long enough
    if len(indi) <= (deg + 1):
      # Skip this subsection
      continue
    
    # Fit polynomial, calculate model and residuals
    poly = np.polyfit(x[indi], y[indi], deg)
    model = np.polyval(poly, x[indi])
    residuals = y[indi] - model
    
    # The number of "degrees of freedom)
    dof = len(indi) - (deg + 1)
    # Calculate reduced chi square ...
    stdEstimate = np.sqrt((residuals**2).sum() / dof)
    # A brute-force way to compute the expectation of sqrt(1/(X(d)/d)),
    # where X(d) is chi-square distributed with d degrees of freedom.
    if dof != lastdof:
      ccSample = np.random.chisquare(dof, 1000)
      corrFac = np.mean(np.sqrt(1.0/(ccSample / dof)))
      lastdof = dof
    stdEstimate *= corrFac 
    # ... and SNR
    snrs.append(np.mean(model) / stdEstimate)
    
    if controlPlot:
      try:
        import matplotlib.pylab as plt
      except ImportError:
        raise(PE.PyARequiredImport("Cannot import matplotlib.pylab", \
              where="estimateSNR", \
              solution=["Install matplotlib.", "Change `controlPlot` to False."]))
      if i == 1:
        # This is the first call
        ax1 = plt.subplot(3,1,1)
        plt.title("Data (blue) and model (red)")
        ax2 = plt.subplot(3,1,2, sharex=ax1)
        plt.title("Residuals")
        ax3 = plt.subplot(3,1,3, sharex=ax1)
        plt.title("SNR")
      # Create plot
      plt.subplot(3,1,1)
      plt.plot(x[indi], y[indi], 'b.')
      plt.plot(x[indi], model, 'r-')
      plt.subplot(3,1,2, sharex=ax1)
      plt.plot(x[indi], residuals, 'g.')
      plt.subplot(3,1,3, sharex=ax1)
      plt.plot(np.mean(x[indi]), snrs[-1], 'bp')
  
  if controlPlot:
    plt.show()
  
  return {"SNR-Estimate":np.mean(snrs), "snrs":snrs}