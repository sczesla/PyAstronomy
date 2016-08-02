# -*- coding: utf-8 -*-
from __future__ import print_function, division
from . import periodBase
import numpy
import scipy.special


class Fourier(periodBase.PeriodBase):
  """
    This class computes the Fast Fourier Transform (FFT) of the input
    data using *numpy*'s FFT routines. The constructor takes the light
    curve, `lc` (*TimeSeries* instance), as input. The optional
    argument specifies the normalization of the Fourier Power.
    Currently, only the normalization according to [Leahy83]_
    is supported, which in the case of purely Poissonian
    noise results in a mean power of 2.
  """

  def __calcPeriodogram(self):

    dt=self.t[1:]-self.t[0:-1]
    if numpy.max(dt) > 1.1 * numpy.min(dt):
      print("timingAna::fourier() - Error: data must be equally spaced!")
      exit()

    power=numpy.abs(numpy.fft.fft(self.y))**2.
    if self.norm=="Leahy":
      power=power*2./numpy.sum(self.y)
    else:
      print("Normalization not yet implemented!")
      power=power*2./numpy.sum(self.y)
    freq = numpy.fft.fftfreq(numpy.size(self.y),numpy.mean(dt))
    indi = numpy.where(freq>0)
    self.freq = freq[indi]
    self.power = power[indi]
    self.M = len(self.freq)

  def __init__(self, lc, norm="Leahy"):
    """
      Parameters:
        - `lc` - *TimeSeries* instance,
                 The light curve to be analyzed.
        - `norm` - optional, string,
                   Normalization method; currently, only default ("Leahy") is supported.
    """

    self.freq = None
    self.power = None
    self.t = lc.time
    self.y = lc.flux
    self.norm = norm
    self.label = {'title': 'Fourier Transformation',\
                  'xlabel': 'Frequency',\
                  'ylabel': 'Fourier Power'}
    self.__calcPeriodogram()

  def Prob(self, Pn):
    """
      Parameters:
        - `Pn` - float, Power threshold.
    
      Returns the probability to obtain a power larger than the threshold, `Pn`, from
      the noise, which is assumed to Poisson-distributed.

      .. note::
        According to [vdK]_ the probability \
        to obtain a power larger than a given threshold from the noise is given by
        
        .. math::
          Prob(p>Pn) = Q(M \\times W \\times Pn, 2 \\times M \\times W)
        
        where :math:`Q(\chi^2, \\nu)` is the cumulative :math:`\\chi^2` distribution with :math:`\\nu` d.o.f.

        
    """
    nu = 2
    Q = scipy.special.chdtrc(nu, Pn) #-- Integral from x to infinity of Chi-square pdf.
    return Q



