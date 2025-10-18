from __future__ import print_function, division
import numpy
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo

class Scanner:
  
  def __init__(self, minVal=0.1, maxVal=10.0, dVal=0.05, mode="period"):
    """
      The Scanner class is used to define the period/frequency range and steps used in the \
      PDM analysis. It is iteratable.
      
      Parameters:
        - `minVal` - float, Minimum value,
        - `maxVal` - float, Maximum value,
        - `dVal`   - float, Delta value,
        - `mode`   - string, optional, Either "period" or "frequency" (default = "period").
      
      .. Note: Whether min, max, and dval refer to period or frequency depends on the ``mode'' parameter.
      .. Note: Time units have match the time axis given to the PDM analysis class.
    """
    if (mode != "period") and (mode != "frequency"):
      raise(PE.PyAValError("`mode` must either be 'period' or 'frequency'. Current value is: "+str(mode), \
                           where="Scanner::__init__", solution="Choose either 'period' or 'frequency' for mode parameter."))
    self.minVal = minVal
    self.maxVal = maxVal
    self.dVal = dVal
    self.mode = mode
  
  def __iter__(self):
    self.curVal = self.minVal
    while self.curVal <= self.maxVal:
      yield self.curVal
      self.curVal += self.dVal



class PyPDM:
  
  def __init__(self, time, mag):
    """
      This class allows to carry out a ``Phase Dispersion Minimization'' (PDM) analysis as described in
      Stellingwerf 1978.
      
      Parameters:
       - time - array, The time array (arbitrary time units),
       - mag - array, The corresponding flux in magnitudes. 
      
      Properties:
       - minBinPoints - The minimum number of data points,
                        which may be contained in an individual bin (default = 3).
    """
    # The minimum number of data points per bin.
    self.minBinPoints = 3
    # Save the data
    self.x = time.copy()
    self.y = mag.copy()
  
  def __setUpEquiBlocks(self, nbins, phases):
    """
      Set up a sequence of equidistant bins. Then check whether all bins contain enough data points \
      (defined by minBinPoints property). \
      If this is not the case, the border to the adjacent bin with fewer photons is dropped. Therefore,
      the resulting bins are not necessarily equidistant anymore.
      
      Parameters:
       - nbins - int, The number of bins.
       - phases - array, The phases of the data points. 
      
      Returns: A list with the beginnings and one list with the ends of the bins (in phase). 
    """
    blockBegin = (numpy.arange(nbins) / float(nbins)).tolist()
    blockEnd   = ((numpy.arange(nbins) + 1.0) / float(nbins)).tolist()
    Ns = []
    for i in smo.range(nbins):
      indi = numpy.where(numpy.logical_and(phases >= blockBegin[i], phases < blockEnd[i]))[0]
      Ns.append(len(indi))
    
    i = 0
    while(i < len(Ns)):
      if Ns[i] < self.minBinPoints:
        # The bin with number i does not contain enough data points.
        # Remove the border to make the bins larger...
        iPlus = i + 1
        iMinu = i - 1
        NPlus = 0; NMinu = 0
        if iPlus == len(blockEnd):
          NPlus = 1e100
        else:
          NPlus = Ns[iPlus]
        if iMinu == -1:
          NMinu = 1e100
        else:
          NMinu = Ns[iMinu]
        if NMinu <= NPlus:
          # Eliminate the border to bin below
          blockBegin.pop(i); blockEnd.pop(i-1)
          Ns[i] = Ns[i-1] + Ns[i]; Ns.pop(i-1)
          # Reiterate
          i = 0
          continue
        else:
          # Eliminate the border to bin behind
          blockBegin.pop(i+1); blockEnd.pop(i)
          Ns[i+1] = Ns[i] + Ns[i+1]; Ns.pop(i)
          # Reiterate
          i = 0
          continue
      i += 1
    return blockBegin, blockEnd

  def __setUpEquiBlocksCover(self, nbins, phases, covers):
    """
      Set up a sequence of equidistant bins. Then reproduce the same sequence but offset by \
      1/(nbins*covers). In this vein, data points are covered by multiple bins. The phase \
      axis is cyclic. This function pays no attention to the number of data points within \
      individual bins.
      
      Parameters:
       - nbins - int, The number of bins.
       - phases - array, The phases of the data points. 
       - covers - int, The number of covers (as explained above).
      
      Returns: A list with the beginnings and one list with the ends of the bins (in phase). 
    """
    blockBegin = numpy.arange(nbins) / float(nbins)
    blockEnd   = (numpy.arange(nbins) + 1.0) / float(nbins)
    
    blockBeginOrig = blockBegin.copy()
    blockEndOrig = blockEnd.copy()
    
    offset = 0
    for i in smo.range(1,covers):
      offset = float(i)/float(nbins * covers)
      blockBegin = numpy.append(blockBegin, blockBeginOrig + offset)
      blockEnd   = numpy.append(blockEnd, blockEndOrig + offset)
      
    return blockBegin, blockEnd

  def phase(self, time, period):
    """
      Calculate phases and ordering.
    
      Parameter:
        - `time` - array, Time stamps,
        - `period` - float, The period used to calculate phases.
    
      Returns:
        An unsorted phase array and an array with indices specifying the order. Thus,
        phase[indi], mag[indi] is the phase-sorted phased light curve.
    """
    phase = time/period - numpy.floor(time/period)
    indi = numpy.argsort(phase)
    return phase, indi
  
  def __getThet(self, phase, mag, bbegin, bend):
    """
      Calculate the Theta statistics defined by Stellinger '78.
      
      Parameters:
       - phase, mag - arrays, Phase-sorted(!) light curve.
       - bbegin, bend -lists, begin and end of individual bins (``blocks'').
      
      Returns the value of Theta.
    """
    meanMag = mag.mean()
    N = len(mag)
    M = len(bbegin)
    sigmaSqr = ((mag - meanMag)**2).sum() / float(N - 1)

    sSqrUp = 0.0; sSqrDown = 0.0
    for i in smo.range(len(bbegin)):
      # Points belonging to a chunk
      indi = numpy.where(numpy.logical_and(phase >= bbegin[i], phase < bend[i]))[0]
      nj = float(len(indi))
      if nj > 1:
        # Variance of individual block
        sigmaj = ((mag[indi] - mag[indi].mean())**2).sum() / float(nj - 1)
        sSqrUp += (nj - 1)*sigmaj
        sSqrDown += nj
    sSqrDown -= M
    
    sSqr = sSqrUp / sSqrDown
    theta = sSqr / sigmaSqr
    return theta
  
  def pdmEquiBin(self, nbins, scanner):
    """
      Carry out the PDM analysis using equidistant bins.
      
      This method uses equidistant bins, yet, it pays attention to the number of \
      data points contained in individual bins. If this number is insufficient \
      (as defined by the `minBinPoints` property), adjacent bins are combined
      (see documentation of *__setUpEquiBlocks*).
      
      Parameters:
        - `nbins` - int, Number of bins to use.
        - `scanner` - An instance of the *Scanner* class, defining, which periods/frequencies
                      are tested.
      
      Returns two arrays: the periods and the associated values of the Theta statistic. These
      are also saved to the class properties `periods`, `frequencies`, and `theta`.

      .. Note:: Whether the first return value is frequency or period
                depends on the mode of the scanner.
    """
    # Go through periods
    periods = []; theta = []
    for val in scanner:
      if scanner.mode == "period":
        period = val
      else:
        period = 1.0/val
      # Note the p[indi], self.y[indi] is the ordered, phased light curve
      p, indi = self.phase(self.x, period)
      bbegin, bend = self.__setUpEquiBlocks(nbins, p[indi])
      periods.append(period)
      theta.append(self.__getThet(p[indi], self.y[indi], bbegin, bend))
    
    # Save result
    self.periods = numpy.array(periods)
    self.frequencies = 1.0 / numpy.array(periods)
    self.theta = numpy.array(theta)
    
    if scanner.mode == "period":
      return self.periods.copy(), self.theta.copy()
    else:
      return self.frequencies.copy(), self.theta.copy()

  def pdmEquiBinCover(self, nbins, covers, scanner):
    """
      Carry out the PDM analysis using multiple sequences of equidistant bins.
      
      The bins used by this method are equidistant, but the phase axis and, thus, \
      the data points are covered by multiple bins. The first sequence of bins is \
      the usual devision of the 0-1 interval into *nbins* bins; the following sequences \
      are the same but offset of 1/(nbins*covers) in phase (the phase axis is cyclic).
      
      This method does not check whether bins contain "enough" data points, i.e., it neglects \
      the `minBinPoints` property taken into account by *pdmEquiBin*. If, however, less than \
      two data points are contained within a bin, it is neglected.
      
      Parameters:
        - `nbins` - int, Number of bins to use.
        - `scanner` - An instance of the *Scanner* class, defining, which periods/frequencies are tested.
        - `covers` - int, The number of covers, i.e., phase-shifted bin sets.
      
      Returns two arrays: the periods/frequencies and the associated values of the Theta statistic. These
      are also saved to the class properties ``periods'', ``frequencies'', and ``theta''.
      
      .. Note:: Whether the first return value is frequency or period depends on the mode of the scanner.
    """
    # Go through periods
    periods = []; theta = []
    for val in scanner:
      if scanner.mode == "period":
        period = val
      else:
        period = 1.0/val
      # Note the p[indi], self.y[indi] is the ordered, phased light curve
      p, indi = self.phase(self.x, period)
      # Extend phase and mag arrays to allow for phases > 1.0 (cyclic phase array...)
      phase = numpy.append(p[indi], p[indi]+1.0)
      mag = numpy.append(self.y[indi], self.y[indi])
      bbegin, bend = self.__setUpEquiBlocksCover(nbins, phase, covers)
      periods.append(period)
      theta.append(self.__getThet(phase, mag, bbegin, bend))
    
    # Save result
    self.periods = numpy.array(periods)
    self.frequencies = 1.0 / numpy.array(periods)
    self.theta = numpy.array(theta)
    
    if scanner.mode == "period":
      return self.periods.copy(), self.theta.copy()
    else:
      return self.frequencies.copy(), self.theta.copy()
    
    
    
    