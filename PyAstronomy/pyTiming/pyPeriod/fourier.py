# -*- coding: utf-8 -*-
from __future__ import print_function, division
from . import periodBase
import numpy
import scipy.special
from PyAstronomy.pyaC import pyaErrors as PE


class Fourier(periodBase.PeriodBase):
    r"""
    This class computes the Fast Fourier Transform (FFT) of the input data.

    This is actually a wrapper around numpy's FFT routines.
    The constructor takes the light
    curve, `lc` (*TimeSeries* instance), as input. The optional
    argument specifies the normalization of the Fourier Power.
    Currently, only the normalization according to [Leahy83]_
    is supported, which in the case of purely Poissonian
    noise results in a mean power of 2.

    Parameters
    ----------
    lc : TimeSeries instance
        The light curve to be analyzed.
    norm : optional, string
        Normalization method; currently, only default ("Leahy") is supported.
    """

    def __calcPeriodogram(self):

        dt = self.t[1:]-self.t[0:-1]
        if numpy.max(dt) > 1.1 * numpy.min(dt):
            raise(PE.PyAValError("Data must be evenly spaced in time.",
                                 where="Fourier",
                                 solution="Use, e.g., the Gls or LombScargle class for unevenly sampled data."))

        # Leahy normalization (the only one available, checked in __init__)
        power = numpy.abs(numpy.fft.fft(self.y))**2.
        power = power*2./numpy.sum(self.y)
        freq = numpy.fft.fftfreq(numpy.size(self.y), numpy.mean(dt))
        indi = numpy.where(freq > 0)
        self.freq = freq[indi]
        self.power = power[indi]
        self.M = len(self.freq)

    def __init__(self, lc, norm="Leahy"):
        self.freq = None
        self.power = None
        self.t = lc.time
        self.y = lc.flux
        if norm != "Leahy":
            raise(PE.PyAValError("Unknown normalization: " + str(norm),
                                 where="Fourier",
                                 solution="Use norm='Leahy' (the only available normalization)."))
        self.norm = norm
        self.label = {'title': 'Fourier Transformation',
                      'xlabel': 'Frequency',
                      'ylabel': 'Fourier Power'}
        self.__calcPeriodogram()

    def prob(self, Pn):
        r"""
        Returns the probability to obtain a power larger than the threshold, `Pn`.

        In the calculations, the noise is assumed to be Poisson-distributed.

        .. note::
            According to [vdK]_, the Leahy-normalized power of pure Poisson
            noise in a single frequency bin of a single spectrum follows a
            :math:`\chi^2` distribution with two degrees of freedom, so that

            .. math::
              Prob(p>Pn) = Q(Pn, 2) = \exp(-Pn/2) \; ,

            where :math:`Q(\chi^2, \nu)` is the probability that a
            :math:`\chi^2` distribution with :math:`\nu` d.o.f. exceeds
            :math:`\chi^2`.

        Parameters
        ----------
        Pn : float or array
            Power threshold.

        Returns
        -------
        Probability : float or array
            The probability to obtain a power larger than the specified
            threshold from noise.
        """
        # Integral from Pn to infinity of the chi-square pdf with 2 d.o.f.
        return scipy.special.chdtrc(2, Pn)

    def Prob(self, Pn):
        """
        Outdated -- use "prob" instead.
        """
        return self.prob(Pn)

    def probInv(self, Prob):
        """
        Returns the minimum power for a given probability level `Prob`.

        This is the inverse of `prob(Pn)`.

        Parameters
        ----------
        Prob : float or array
            Probability.

        Returns
        -------
        Power threshold : float or array
            The power exceeded by noise with probability `Prob`.
        """
        return scipy.special.chdtri(2, Prob)
