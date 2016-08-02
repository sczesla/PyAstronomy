# -*- coding: utf-8 -*-
from __future__ import print_function, division
from . import periodBase
from numpy import var,mean,min,max,pi,sin,cos,sum,arange,zeros,arctan
from PyAstronomy.pyaC import pyaErrors as PE


class LombScargle(periodBase.PeriodBase):
  """
    The constructor of *LombScargle* takes a *TimeSeries* instance, i.e., a 
    light curve object, as first argument. It then computes the usual 
    Lomb-Scargle periodogram using a fast algorithm. The frequency array is
    constructed on the fly based on the oversampling keywords, which are
    mandatory in this case. The power is normalized according to the
    prescription of [HB86]_.
  """

  def __calcPeriodogram(self):

    ofac=int(self.ofac); hifac=float(self.hifac)

    if self.freq is None:
      self.__buildFreq(ofac,hifac)
    nout = int(ofac*hifac*len(self.t)/2.)
    vari = var(self.y)
    ave = mean(self.y)

    xdif = max(self.t)-min(self.t)
    xave = 0.5*(max(self.t)+min(self.t))

    # Get starting frequency
    pnow = 1./(xdif*ofac)

    arg = 2.*pi*(self.t-xave)*pnow
    wpr = -2.*sin(0.5*arg)**2
    wpi = sin(arg)
    wr = cos(arg)
    wi = wpi

    py = zeros(nout)
    yy = self.y-ave

    for i in range(nout):
      sumsh = sum(wi*wr)
      sumc = sum((wr-wi)*(wr+wi))
      wtau = 0.5 * arctan(2.*sumsh/sumc)
      swtau = sin(wtau)
      cwtau = cos(wtau)
      ss = wi*cwtau-wr*swtau
      cc = wr*cwtau+wi*swtau
      sums = sum(ss**2)
      sumc = sum(cc**2)
      sumsy = sum(yy*ss)
      sumcy = sum(yy*cc)
      wtemp = wr
      wr   =(wtemp*wpr-wi*wpi)+wr
      wi   =wi*wpr+wtemp*wpi+wi
      py[i]=0.5*(sumcy**2/sumc+sumsy**2/sums)/vari

    self.power = py.copy()
    self.M = 2.*nout/ofac
    self.N = len(self.y)

  def __buildFreq(self,ofac,hifac):
    nout = int(ofac*hifac*len(self.t)/2.)
    xdif = max(self.t)-min(self.t)
    pnow = 1./(xdif*ofac)
    px = pnow + arange(nout)/(ofac*xdif)
    self.freq = px
    return 0

  def __init__(self, lc, ofac, hifac):
    """
      Parameters:
       - `lc` - TimesSeries instance, The light curve to be analyzed.
       - `ofac` - int, Oversampling factor.
       - `hifac` - float, Maximum frequency `freq` = `hifac` * (average Nyquist frequency).
      
      .. note::
        Adapted from routine of the same routine in [NR]_ , 
        based on period.pro by Han Wen, August 1996.
      
      The result, i.e., the power, is stored in the class property `power`.
    """

    self.freq=None
    self.power=None
    self.ofac=ofac
    self.hifac=hifac
    self.t=lc.time
    self.y=lc.flux
    self.label = {'title': 'Lomb-Scargle periodogram', \
                  'xlabel': 'Frequency', \
                  'ylabel': 'Scargle Power'}
    self.__calcPeriodogram()

  def prob(self, Pn):
    """
      Returns the probability to obtain a power *Pn* or larger from the noise,
      which is assumes to be Gaussian.
      
      Parameters:
        - `Pn` - float, Power threshold.

      .. note::
        *LombScargle* calculates the quantity (N-1)/2.*p=p' (in the formalism of 
        [ZK09]_), which is de facto the normalization
        prescription of [HB86]_. In this
        scheme the probability P(p'>Pn) is given by the following statement:
        
        .. math::
          P(p'>Pn) = \\left(1 - 2 \\frac{Pn}{N-1} \\right)^{(N-3)/2}

        If properly normalized to the population variance of the time series, which
        must be known a priori (usually not the case), the
        power :math:`p/p_n=p"` is a direct measure of the SNR as proposed by [Scargle82]_:
        
        .. math::
          P(p">Pn) = exp(-Pn) \\; .
        
        This formula is often used erroneously in this context.
    """
    return (1.-2.*Pn/(self.N-1.))**((self.N-3.)/2.)

  def Prob(self, Pn):
    """
      Outdated -- use "prob" instaed.
    """
    return self.prob(Pn)

  def probInv(self, Prob):
    """
      Inverse of `Prob(Pn)`. Returns the minimum power
      for a given probability level `Prob`.
      
      Parameters:
        - `Prob` - float, probability
    """
    return (self.N-1.)/2.*(1.-Prob**(2./(self.N-3.)))

  def ProbInv(self, Prob):
    """
      Outdated -- use "probInv" instaed
    """
    return self.probInv(Prob)