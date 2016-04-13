# -*- coding: utf-8 -*-
from __future__ import print_function, division
from . import periodBase
from numpy import sum, cos, sin, arange, min, max, pi, \
                  exp, log, zeros, argmax, sqrt, arctan2
from PyAstronomy.pyaC import pyaErrors as PE 


class Gls(periodBase.PeriodBase):
  """
    Compute the Generalized Lomb-Scargle (GLS) periodogram.
  
    The *Gls* class computes the error-weighted Lomb-Scargle periodogram as
    developed by [ZK09]_ using various possible
    normalizations.

    The constructor of *Gls* takes a *TimeSeries* instance (i.e., a light curve)
    as first argument. As the algorithm is pretty slow, the constructor either
    allows to pass `ofac` and/or `hifac` keywords or an optional `freq` array,
    which will be used to calculate the periodogram.
    
    The main result of the calculation, i.e., the powers, are stored
    in the class property `power`.

    Parameters
    ----------
    lc : TimeSeries
        The light curve or any object providing the attributes
        time, flux, and error, which define the data.
    ofac : int
        Oversampling factor (default=10).
    hifac : float
        Maximum frequency `freq` = `hifac` * (average Nyquist frequency)
        (default=1).
    freq : array, optional
        Contains the frequencies at which to calculate the periodogram.
        If not given, a frequency array will be automatically generated.
    norm : string, optional
        The normalization; either of "Scargle", "HorneBaliunas", "Cumming".
        The default is "HorneBaliunas".
    stats : boolean, optional
        Set True to obtain some statistical output (default is False).
    
    Attributes
    ----------
    power : array
        The normalized power of the GLS.
    freq : array
        The frequency array.
    ofac : int
        The oversampling factor.
    hifac : float
        The maximum frequency.
    t : array
        The abscissa data values.
    y : array
        The ordinate data values.
    error : array
        The errors of the data values.
    norm : string
        The used normalization.
  """

  def _plot(self):
    """
      Create a plot.
    """
    try:
      import matplotlib.pylab as plt
    except ImportError:
      raise(PE.PyARequiredImport("Could not import matplotlib.pylab."))
    self.fig = plt.figure()
    self.fig.subplots_adjust(hspace=0.35)
    self.ax = self.fig.add_subplot(2,1,1)
    self.ax.set_title("Normalized periodogram")
    self.ax.set_xlabel("Frequency")
    self.ax.set_ylabel("Power")
    self.ax.plot(self.freq, self.power, 'b.--')

    self.ax = self.fig.add_subplot(2,1,2)
    self.ax.set_title("Data")
    self.ax.set_xlabel("Time")
    self.ax.set_ylabel("Data")
    self.ax.plot(self.t, self.y, 'r.--')
    
    plt.show()

  def _output(self):
    """
      Some statistical output.
    """
    # Index with maximum power
    bbin = argmax(self.power)
    # Maximum power
    pmax = self._upow[bbin]
    
    rms = sqrt(self._YY * (1.-pmax))
    
    # Get the curvature in the power peak by fitting a parabola y=aa*x^2
    if (bbin > 1) and (bbin < len(self.freq)-2):
      # Shift the parabola origin to power peak
      xh = (self.freq[bbin-1:bbin+2] - self.freq[bbin])**2
      yh = self._upow[bbin-1:bbin+2] - self._upow[bbin]
      # Calculate the curvature (final equation from least square)
      aa = sum(yh*xh)/sum(xh*xh)
      nt = float(self.N)
      f_err = sqrt(-2./nt * pmax/aa*(1.-pmax)/pmax)
      Psin_err = sqrt(-2./nt* pmax/aa*(1.-pmax)/pmax) / self.freq[bbin]**2
    else:
      f_err = None
      Psin_err = None
      PE.warn(PE.PyAValError("WARNING: Highest peak is at the edge of the frequency range.\nNo output of frequency error.\nIncrease frequency range to sample the peak maximum."))
    
    fbest = self.freq[bbin]
    amp = sqrt(self._a[bbin]**2 + self._b[bbin]**2)
    ph  = arctan2(self._a[bbin], self._b[bbin]) / (2.*pi)
    T0  = min(self.th) - ph/fbest
    # Re-add the mean
    offset = self._off[bbin] + self._Y
    
    # Statistics
    print("Generalized LS - statistical output")
    print("-----------------------------------")
    print("Number of input points:     %6d" % (nt))
    print("Weighted mean of dataset:   % e" % (self._Y))
    print("Weighted rms of dataset:    % e" % (sqrt(self._YY)))
    print("Time base:                  % e" % (max(self.th) - min(self.th)))
    print("Number of frequency points: %6d" % (len(self.freq)))
    print()
    print("Maximum power, p :    % e " % (self.power[bbin]))
    print("Maximum power (without normalization):   %e" % (pmax))
    print("Normalization    : ", self.norm)
    print("RMS of residuals :    % e " % (rms))
    if self.error is not None:
      print("  Mean weighted internal error:  % e" %(sqrt(nt/sum(1./self.error**2))))
    print("Best sine frequency : % e +/- % e" % (fbest, f_err))
    print("Best sine period    : % e +/- % e" % (1./fbest, Psin_err))
    print("Amplitude:          : % e +/- % e" % (amp, sqrt(2./nt)*rms))
    print("Phase (ph) : % e +/- % e" % (ph, sqrt(2./nt)*rms/amp/(2.*pi)))
    print("Phase (T0) : % e +/- % e" % (T0, sqrt(2./nt)*rms/amp/(2.*pi)/fbest))
    print("Offset     : % e +/- % e" % (offset, sqrt(1./nt)*rms))
    print("-----------------------------------")

  def __calcPeriodogram(self):

    # Build frequency array if not present
    if self.freq is None:
      self.__buildFreq()
    # Circular frequencies
    omegas = 2.*pi * self.freq

    err2 = self.error*self.error
    w = (1. / sum(1./err2)) / err2

    self._Y = sum(w*self.y)        # Eq. (7)
    yh = self.y - self._Y          # Subtract weighted mean
    self._YY = sum(w * yh**2)
    
    # Unnormalized power
    self._upow = zeros(len(omegas))
    self._a = zeros(len(omegas))
    self._b = zeros(len(omegas))
    self._off = zeros(len(omegas))
    for i, omega in enumerate(omegas):
      x = omega*self.th
      cosx = cos(x)
      sinx = sin(x)
      wcosx = w*cosx         # attach weights
      wsinx = w*sinx         # attach weights
      
      C = sum(wcosx)         # Eq. (8)
      S = sum(wsinx)         # Eq. (9)

      YC = sum(yh*wcosx)     # Eq. (11)
      YS = sum(yh*wsinx)     # Eq. (12)
      CCh = sum(wcosx*cosx)  # Eq. (13)
      CSh = sum(wcosx*sinx)  # Eq. (15)
      SSh = 1.-CCh
      CC = CCh-C*C           # Eq. (13)
      SS = SSh-S*S           # Eq. (14)
      CS = CSh-C*S           # Eq. (15)
      D = CC*SS-CS*CS        # Eq. (6)
      
      self._a[i] = (YC*SS-YS*CS) / D
      self._b[i] = (YS*CC-YC*CS) / D
      self._off[i] = -self._a[i]*C - self._b[i]*S

      self._upow[i] = (SS*YC*YC + CC*YS*YS - 2.*CS*YC*YS) / (self._YY*D) # Eq. (5) in ZK09

    self.N = len(self.y)
    # An ad-hoc estimate of the number of independent frequencies (ZK_09 Eq. 24)
    self.M = (max(self.freq)-min(self.freq)) * (self.th.max() - self.th.min())

    # Normalization:
    if self.norm == "Scargle":
      popvar=input('pyTiming::gls - Input a priori known population variance:')
      self.power = self._upow/float(popvar)
    if self.norm == "HorneBaliunas":
      self.power = (self.N-1.)/2.*self._upow
    if self.norm == "Cumming":
      self.power = (self.N-3.)/2. * self._upow/(1.-max(self._upow))
    
    # Output statistics
    if self._stats:
      self._output()
    if self._showPlot:
      self._plot()

  def __buildFreq(self):
    """
      Build frequency array (`freq` attribute).
    """
    nout = self.ofac * self.hifac * len(self.th)/2
    xdif = max(self.th)-min(self.th)
    self.freq = 1./(xdif*self.ofac) + arange(nout)/(self.ofac*xdif)

  def prob(self, Pn):
    """
      Probability of obtaining the given power.
    
      Calculate the probability to obtain a power higher than
      `Pn` from the noise, which is assumed to be Gaussian.
      
      .. note:: Normalization
        (see [ZK09]_ for further details).

        - `Scargle`: 
        .. math::
          exp(-Pn)

        - `HorneBaliunas`: 
        .. math::
          \\left(1 - 2 \\times \\frac{Pn}{N-1} \\right)^{(N-3)/2}
        
        - `Cumming`: 
        .. math:: 
          \\left(1+2\\times \\frac{Pn}{N-3}\\right)^{-(N-3)/2}
      
      Parameters
      ----------
      Pn : float
          Power threshold.
      
      Returns
      -------
      Probability : float
          The probability to obtain a power equal or
          higher than the threshold from the noise.
    """
    if self.norm=="Scargle": return exp(-Pn)
    if self.norm=="HorneBaliunas": return (1.-2.*Pn/(self.N-1.))**((self.N-3.)/2.)
    if self.norm=="Cumming": return (1.+2.*Pn/(self.N-3.))**(-(self.N-3.)/2.)

  def Prob(self, Pn):
    """
      Outdated---use `prob` instead. 
    """
    return self.prob(Pn)

  def probInv(self, Prob):
    """
      Calculate minimum power for given probability.
    
      This function is the inverse of `Prob(Pn)`.
      Returns the minimum power
      for a given probability threshold `Prob`.
      
      Parameters
      ----------
      Prob : float
          Probability threshold.
      
      Returns
      -------
      Power threshold : float
          The minimum power for the given
          false-alarm probability threshold.
    """
    if self.norm=="Scargle": return -log(Prob)
    if self.norm=="HorneBaliunas": return (self.N-1.)/2.*(1.-Prob**(2./(self.N-3.)))
    if self.norm=="Cumming": return (self.N-3.)/2.*(Prob**(-2./(self.N-3.))-1.)

  def ProbInv(self, Prob):
    """
      Outdated---use `probInv` instead. 
    """
    return self.probInv(Prob)
  
  def __init__(self, lc, ofac=10, hifac=1, freq=None, norm="HorneBaliunas", stats=False, plot=False):

    if lc.error is None:
      raise(PE.PyAValError("No errors given for light curve", where="Gls::__init__", \
            solution="Specify 'error' property in light curve object."))
      
    self.power = None
    self.freq = freq
    self.ofac = ofac
    self.hifac = hifac
    self.t = lc.time
    self.th = lc.time - min(self.t)
    self.y = lc.flux
    self.error = lc.error
    self.norm = norm
    self._showPlot = plot
    # Check and assign normalization
    self.label = {'title': 'Generalized Lomb Periodogram',\
                  'xlabel': 'Frequency',\
                  'ylabel': "Normalization not implemented!"}
    if self.norm == "Scargle":
      self.label["ylabel"] = "Normalized Power (Scargle 1982)"
    elif self.norm == "HorneBaliunas":
      self.label["ylabel"] = "Normalized Power (Horne & Baliunas 1986)"
    elif self.norm == "Cumming":
      self.label["ylabel"] = "Normalized Power (Cumming 1999)"
    else:
      raise(PE.PyAValError("Unknown norm: "+str(self.norm), \
            solution="Use either of Scargle, HorneBaliunas, or Cumming", \
            where="GLS"))
    
    self._stats = stats
    self.__calcPeriodogram()