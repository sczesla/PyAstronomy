# -*- coding: utf-8 -*-

# Copyright (c) 2012 Sebastian Schröter, Stefan Czesla, and Mathias Zechmeister

# The MIT License (MIT)

# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from __future__ import print_function, division
from . import periodBase
from numpy import sum, cos, sin, arange, min, max, pi, \
                  exp, log, zeros, argmax, sqrt, arctan2, dot
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE 



class Gls:
    """
    Compute the Generalized Lomb-Scargle (GLS) periodogram.

    The *Gls* class computes the error-weighted Lomb-Scargle periodogram as
    developed by [ZK09]_ using various possible normalizations.

    The constructor of *Gls* takes a *TimeSeries* instance (i.e., a light curve)
    as first argument. The constructor allows to pass keywords to adjust the
    `freq` array, which will be used to calculate the periodogram.

    The main result of the calculation, i.e., the power, are stored in the
    class property `power`.

    Parameters
    ----------
    lc : tuple or list or TimeSeries object
        The light curve data either in the form of a TimeSeries object (or any
        object providing the attributes time, flux, and error) or a tuple or list
        providing time as first element, flux as second element, and optionally,
        the error as third element.
    fbeg, fend : float, optional
        The beginning and end frequencies for the periodogram
        (inverse units of time axis).
    Pbeg, Pend : float, optional
        The beginning and end periods for the periodogram
        (same units as for time axis).
    ofac : int
        Oversampling factor of frequency grid (default=10).
    hifac : float
        Maximum frequency `freq` = `hifac` * (average Nyquist frequency)
        (default=1).
    freq : array, optional
        Contains the frequencies at which to calculate the periodogram.
        If given, fast and verbose option are not available.
        If not given, a frequency array will be automatically generated.
    norm : string, optional
        The normalization; either of "ZK", "Scargle", "HorneBaliunas", "Cumming", "wrms", "chisq".
        The default is unity ("ZK").
    ls : boolean, optional
        If True, the conventional Lomb-Scargle periodogram will be computed
        (default is False).
    fast : boolean, optional
        If True, recursive relations for trigonometric functions will be used
        leading to faster evaluation (default is False).
    verbose : boolean, optional
        Set True to obtain some statistical output (default is False).

    Attributes
    ----------
    power : array
        The normalized power of the GLS.
    freq : array
        The frequency array.
    ofac : int
        The oversampling factor of frequency grid.
    hifac : float
        The maximum frequency.
    t : array
        The abscissa data values.
    y : array
        The ordinate data values.
    e_y : array
        The errors of the data values.
    norm : string, {'ZK', 'Scargle', 'HorneBaliunas', 'Cumming', 'wrms', 'chisq'}
        The used normalization.

    Examples
    --------
    Create 1000 unevenly sampled data points with frequency=0.1,
    measurement error and Gaussian noise
    
    >>> time = np.random.uniform(54000., 56000., 1000)
    >>> flux = 0.15 * np.sin(2. * np.pi * time / 10.)

    Add some noise
    
    >>> error = 0.5 * np.ones(time.size)
    >>> flux += np.random.normal(0, error)

    Compute the full error-weighted Lomb-Periodogram
    in 'ZK' normalization and calculate the significance
    of the maximum peak.
    
    >>> gls = Gls((time, flux, error), verbose=True)

    >>> maxPower = gls.pmax
    >>> print("GLS maximum power: ", maxPower)
    >>> print("GLS statistics of maximum power peak: ", gls.stats(maxPower))
    >>> gls.plot(block=True)

    """
    # Available normalizations
    norms = ['ZK', 'Scargle', 'HorneBaliunas', 'Cumming', 'wrms', 'chisq', 'lnL', 'dlnL']

    def __init__(self, lc, fbeg=None, fend=None, Pbeg=None, Pend=None, ofac=10, hifac=1, freq=None, norm="ZK", ls=False, fast=False, verbose=False, **kwargs):

        self.freq = freq
        self.fbeg = fbeg
        self.fend = fend
        self.Pbeg = Pbeg
        self.Pend = Pend
        self.ofac = ofac
        self.hifac = hifac
        self.ls = ls
        self.norm = norm
        self.fast = fast
        self.label = {'title': 'Generalized Lomb Periodogram',
                      'xlabel': 'Frequency'}
        if "stats" in kwargs:
            print("Warning: 'stats' option is outdated. Please use 'verbose' instead.")
            verbose = kwargs["stats"]

        self._normcheck(norm)

        self._assignTimeSeries(lc)
        self._buildFreq()
        self._calcPeriodogram()
        self.pnorm(norm)
        self._peakPeriodogram()

        # Output statistics
        if verbose:
            self.info()

    def _sete_y(self, ey):
        self.yerr = ey
    
    def _gete_y(self):
        return self.yerr

    """ Treat e_y and yerr synonymously for backward compatibility """
    e_y = property(_gete_y, _sete_y)
    

    def _assignTimeSeries(self, lc):
        """
        A container class that holds the observed light curve.

        Parameters
        ----------
        time : array
            The time array.
        flux : array
            The observed flux/data.
        error : array, optional
            The error of the data values.

        """
        if isinstance(lc, (tuple, list)):
            # t, y[, e_y] were given as list or tuple.
            if len(lc) in (2, 3):
                self.t = np.ravel(lc[0])
                self.y = np.ravel(lc[1])
                self.e_y = None
                if len(lc) == 3 and lc[2] is not None:
                    # Error has been specified.
                    self.e_y = np.ravel(lc[2])
            else:
                raise(ValueError("lc is a list or tuple with " + str(len(lc)) + " elements. Needs to have 2 or 3 elements." + \
                                   " solution=Use 2 or 3 elements (t, y[, e_y]) or an instance of TimeSeries"))
        else:
            # Assume lc is an instance of TimeSeries.
            self.t, self.y, self.e_y = lc.time, lc.flux, lc.error

        self.th = self.t - self.t.min()
        self.tbase = self.th.max()
        self.N = len(self.y)

        # Re-check array length compatibility
        if (len(self.th) != self.N) or ((self.e_y is not None) and (len(self.e_y) != self.N)):
            raise(ValueError("Incompatible dimensions of input data arrays (time and flux [and error]). Current shapes are: " + \
                             ', '.join(str(np.shape(x)) for x in (self.t, self.y, self.e_y))))

    def _buildFreq(self):
        """
        Build frequency array (`freq` attribute).

        Attributes
        ----------
        fnyq : float
            Half of the average sampling frequency of the time series.

        """
        self.fstep = 1 / self.tbase / self.ofac   # frequency sampling depends on the time span, default for start frequency
        self.fnyq = 0.5 / self.tbase * self.N     # Nyquist frequency
        self.f = self.freq

        if self.freq is None:
            if self.fast:
                raise(ValueError("freq and fast cannot be used together."))
            # Build frequency array if not present.
            if self.fbeg is None:
                self.fbeg = self.fstep if self.Pend is None else 1 / self.Pend
            if self.fend is None:
                self.fend = self.fnyq * self.hifac if self.Pbeg is None else 1 / self.Pbeg

            if self.fend <= self.fbeg:
                raise(ValueError("fend is smaller than (or equal to) fbeg but it must be larger." + \
                               "Choose fbeg and fend so that fend > fbeg."))

            self.freq = arange(self.fbeg, self.fend, self.fstep)
        else:
            self.fbeg, self.fend = np.min(self.freq), np.max(self.freq)

        self.nf = len(self.freq)

        # An ad-hoc estimate of the number of independent frequencies (Eq. (24) in ZK_09).
        self.M = (self.fend-self.fbeg) * self.tbase

    def _calcPeriodogram(self):

        if self.e_y is None:
            w = np.ones(self.N)
        else:
            w = 1 / (self.e_y * self.e_y)
        self.wsum = w.sum()
        w /= self.wsum

        self._Y = dot(w, self.y)       # Eq. (7)
        wy = self.y - self._Y          # Subtract weighted mean
        self._YY = dot(w, wy**2)       # Eq. (10)
        wy *= w                        # attach errors

        C, S, YC, YS, CC, CS = np.zeros((6, self.nf))

        if self.fast:
            # Prepare trigonometric recurrences.
            eid = exp(2j * pi * self.fstep * self.th)  # cos(dx)+i sin(dx)

        for k, omega in enumerate(2.*pi*self.freq):
            # Circular frequencies.
            if self.fast:
                if k % 1000 == 0:
                    # init/refresh recurrences to stop error propagation
                    eix = exp(1j * omega * self.th)  # exp(ix) = cos(x) + i*sin(x)
                cosx = eix.real
                sinx = eix.imag
                eix *= eid              # increase freq for next loop
            else:
                x = omega * self.th
                cosx = cos(x)
                sinx = sin(x)

            C[k] = dot(w, cosx)         # Eq. (8)
            S[k] = dot(w, sinx)         # Eq. (9)

            YC[k] = dot(wy, cosx)       # Eq. (11)
            YS[k] = dot(wy, sinx)       # Eq. (12)
            wcosx = w * cosx
            CC[k] = dot(wcosx, cosx)    # Eq. (13)
            CS[k] = dot(wcosx, sinx)    # Eq. (15)

        SS = 1. - CC
        if not self.ls:
            CC -= C * C            # Eq. (13)
            SS -= S * S            # Eq. (14)
            CS -= C * S            # Eq. (15)
        D = CC*SS - CS*CS          # Eq. (6)

        self._a = (YC*SS-YS*CS) / D
        self._b = (YS*CC-YC*CS) / D
        self._off = -self._a*C - self._b*S

        # power
        self.p = (SS*YC*YC + CC*YS*YS - 2.*CS*YC*YS) / (self._YY*D)   # Eq. (5) in ZK09

    def _normcheck(self, norm):
        """
        Check normalization

        Parameters
        ----------
        norm : string
            Normalization string

        """
        if norm not in self.norms:
            raise(ValueError("Unknown norm: " + str(norm) + ". " + \
                "Use either of " + ', '.join(self.norms)))

    def pnorm(self, norm="ZK"):
        """
        Assign or modify normalization (can be done afterwards).

        Parameters
        ----------
        norm : string, optional
            The normalization to be used (default is 'ZK').

        Examples
        --------
        >>> gls.pnorm('wrms')

        """
        self._normcheck(norm)
        self.norm = norm
        p = self.p
        power = p   # default ZK
        self.label["ylabel"] = "Power ("+norm+")"

        if norm == "Scargle":
            popvar = input('pyTiming::gls - Input a priori known population variance:')
            power = p / float(popvar)
        elif norm == "HorneBaliunas":
            power = (self.N-1)/2. * p
        elif norm == "Cumming":
            power = (self.N-3)/2. * p / (1.-self.p.max())
        elif norm == "chisq":
            power = self._YY *self.wsum * (1.-p)
            self.label["ylabel"] = "chisq"
        elif norm == "wrms":
            power = sqrt(self._YY*(1.-p))
            self.label["ylabel"] = "wrms"
        elif norm == "lnL":
            chi2 = self._YY *self.wsum * (1.-p)
            power = -0.5*chi2 - 0.5*np.sum(np.log(2*np.pi * self.e_y))
            self.label["ylabel"] = "lnL"
        elif norm == "dlnL":
            # dlnL = lnL - lnL0 = -0.5 chi^2 + 0.5 chi0^2 = 0.5 (chi0^2 - chi^2) = 0.5 chi0^2 p
            power = 0.5 * self._YY * self.wsum * p
            self.label["ylabel"] = "dlnL"

        self.power = power

    def _peakPeriodogram(self):
        """
        Analyze the highest periodogram peak.
        """
        # Index with maximum power
        k = argmax(self.p)
        # Maximum power
        self.pmax = pmax = self.p[k]
        self.rms = rms = sqrt(self._YY*(1.-pmax))

        # Statistics of highest peak
        self.hpstat = p = {}

        # Best parameters
        p["fbest"] = fbest = self.freq[k]
        p["amp"] = amp = sqrt(self._a[k]**2 + self._b[k]**2)
        p["ph"] = ph = arctan2(self._a[k], self._b[k]) / (2.*pi)
        p["T0"]  = self.t.min() - ph/fbest
        p["offset"] = self._off[k] + self._Y            # Re-add the mean.

        # Error estimates
        p["amp_err"] = sqrt(2./self.N) * rms
        p["ph_err"] = ph_err = sqrt(2./self.N) * rms/amp/(2.*pi)
        p["T0_err"] = ph_err / fbest
        p["offset_err"] = sqrt(1./self.N) * rms

        # Get the curvature in the power peak by fitting a parabola y=aa*x^2
        if 1 < k < self.nf-2:
            # Shift the parabola origin to power peak
            xh = (self.freq[k-1:k+2] - self.freq[k])**2
            yh = self.p[k-1:k+2] - pmax
            # Calculate the curvature (final equation from least square)
            aa = dot(yh, xh) / dot(xh, xh)
            p["f_err"] = e_f = sqrt(-2./self.N / aa * (1.-self.pmax))
            p["Psin_err"] = e_f / fbest**2
        else:
            self.hpstat["f_err"] = np.nan
            self.hpstat["Psin_err"] = np.nan
            print("WARNING: Highest peak is at the edge of the frequency range.\nNo output of frequency error.\nIncrease frequency range to sample the peak maximum.")

    def sinmod(self, t):
        """
        Calculate best-fit sine curve.

        The parameters of the best-fit sine curve can be accessed via
        the dictionary attribute `hpstat`. Specifically, "amp" holds the
        amplitude, "fbest" the best-fit frequency, "T0" the reference time
        (i.e., time of zero phase), and "offset" holds the additive offset
        of the sine wave. 

        Parameters
        ----------
        t : array
            Time array at which to calculate the sine.

        Returns
        -------
        Sine curve : array
            The best-fit sine curve (i.e., that for which the
            power is maximal).
        """
        try:
            p = self.hpstat
            return p["amp"] * sin(2*np.pi*p["fbest"]*(t-p["T0"])) + p["offset"]
        except Exception as e:
            print("Failed to calculate best-fit sine curve.")
            raise(e)

    def info(self):
        """
        Prints some basic statistical output screen.
        """
        print("Generalized LS - statistical output")
        print("-----------------------------------")
        print("Number of input points:      %-6d" % self.N)
        print("Weighted mean of dataset:   % f"  % self._Y)
        print("Weighted rms of dataset:    % f"  % sqrt(self._YY))
        print("Time base:                  % f"  % self.tbase)
        print("Number of frequency points:  %-6d" % self.nf)
        print()
        print("Maximum power p [%s]: % f" % (self.norm, self.power.max()))
        print("RMS of residuals:     % f" % self.rms)
        if self.e_y is not None:
            print("  Mean weighted internal error:  %f" % (sqrt(self.N/sum(1./self.e_y**2))))
        print("Best sine frequency:  % f +/- % f" % (self.hpstat["fbest"], self.hpstat["f_err"]))
        print("Best sine period:     % f +/- % f" % (1./self.hpstat["fbest"], self.hpstat["Psin_err"]))
        print("Amplitude:            % f +/- % f" % (self.hpstat["amp"], self.hpstat["amp_err"]))
        print("Phase (ph):           % f +/- % f" % (self.hpstat["ph"], self.hpstat["ph_err"]))
        print("Phase (T0):           % f +/- % f" % (self.hpstat["T0"], self.hpstat["T0_err"]))
        print("Offset:               % f +/- % f" % (self.hpstat["offset"], self.hpstat["offset_err"]))
        print("-----------------------------------")

    def plot(self, block=False, period=False):
        """
        Create a plot.
        """
        try:
            import matplotlib
            if (matplotlib.get_backend() != "TkAgg"):
                matplotlib.use("TkAgg")
            import matplotlib.pylab as plt
            from matplotlib.ticker import FormatStrFormatter
        except ImportError:
            raise(ImportError("Could not import matplotlib.pylab."))

        fig = plt.figure()
        fig.subplots_adjust(hspace=0.15, wspace=0.08, right=0.97, top=0.95)
        ax = fig.add_subplot(3, 1, 1)
        ax.set_title("Normalized periodogram")
        if period:
            ax.set_xscale("log")
            ax.set_xlabel("Period")
        else:
            ax.set_xlabel("Frequency")
        ax.set_ylabel(self.label["ylabel"])
        ax.plot(1/self.freq if period else self.freq, self.power, 'b-')

        fbest, T0 = self.hpstat["fbest"], self.hpstat["T0"]
        # Data and model
        datstyle = {'yerr':self.e_y, 'fmt':'r.', 'capsize':0}
        tt = arange(self.t.min(), self.t.max(), 0.01/fbest)
        ymod = self.sinmod(tt)
        yfit = self.sinmod(self.t)
        ax1 = fig.add_subplot(3, 2, 3)
        # ax1.set_xlabel("Time")
        ax1.set_ylabel("Data")
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.errorbar(self.t, self.y, **datstyle)
        ax1.plot(tt, ymod, 'b-')

        tt = arange(T0, T0+1/fbest, 0.01/fbest)
        yy = self.sinmod(tt)
        ax2 = fig.add_subplot(3, 2, 4, sharey=ax1)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax2.get_yticklabels(), visible=False)
        # ax2.set_xlabel("Time")
        # ax2.set_ylabel("Data")
        ax2.errorbar(self.t*fbest % 1, self.y, **datstyle)
        xx = tt*fbest % 1
        ii = np.argsort(xx)
        ax2.plot(xx[ii], yy[ii], 'b-')

        # Residuals
        yres = self.y - yfit
        ax3 = fig.add_subplot(3, 2, 5, sharex=ax1)
        ax3.set_xlabel("Time")
        ax3.set_ylabel("Residuals")
        ax3.errorbar(self.t, yres, **datstyle)
        ax3.plot([self.t.min(), self.t.max()], [0,0], 'b-')

        ax4 = fig.add_subplot(3, 2, 6, sharex=ax2, sharey=ax3)
        # ax4.set_title("Data")
        ax4.set_xlabel("Phase")
        # ax4.set_ylabel("Data")
        plt.setp(ax4.get_yticklabels(), visible=False)
        ax4.errorbar(self.t*fbest % 1, yres, **datstyle)
        ax4.plot([0,1], [0,0], 'b-')

        if hasattr(plt.get_current_fig_manager(), 'toolbar'):
            # check seems not needed when "TkAgg" is set
            plt.get_current_fig_manager().toolbar.pan()
        #t = fig.canvas.toolbar
        #plt.ToggleTool(plt.wx_ids['Pan'], False)
        if block:
            print("Close the plot to continue.")
        else:
            plt.ion()
        plt.show()
        # plt.show(block=block) # unexpected keyword argument 'block' in older matplotlib
        return plt

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
        self._normcheck(self.norm)
        if self.norm == "ZK": return (1.-Pn)**((self.N-3.)/2.)
        if self.norm == "Scargle": return exp(-Pn)
        if self.norm == "HorneBaliunas": return (1-2*Pn/(self.N-1)) ** ((self.N-3)/2)
        if self.norm == "Cumming": return (1+2*Pn/(self.N-3)) ** (-(self.N-3)/2)
        if self.norm == "wrms": return (Pn**2/self._YY) ** ((self.N-3)/2)
        if self.norm == "chisq": return (Pn/self._YY/self.wsum) ** ((self.N-3)/2)
        if self.norm == "ZK":
            p = Pn
        if self.norm == "dlnL":
            p = 2 * Pn / self._YY / self.wsum
        if self.norm == "lnL":
            chi2 = -2*Pn - np.sum(np.log(2*np.pi * self.e_y**2))
            p = 1 - chi2/self._YY/self.wsum
        return (1-p) ** ((self.N-3)/2)


    def probInv(self, Prob):
        """
        Calculate minimum power for given probability.

        This function is the inverse of `Prob(Pn)`.
        Returns the minimum power for a given probability threshold `Prob`.

        Parameters
        ----------
        Prob : float
            Probability threshold.

        Returns
        -------
        Power threshold : float
            The minimum power for the given false-alarm probability threshold.

        """
        self._normcheck(self.norm)
        if self.norm == "ZK": return 1.-Prob**(2./(self.N-3.))
        if self.norm == "Scargle": return -log(Prob)
        if self.norm == "HorneBaliunas": return (self.N-1) / 2 * (1-Prob**(2/(self.N-3)))
        if self.norm == "Cumming": return (self.N-3) / 2 * (Prob**(-2./(self.N-3))-1)
        if self.norm == "wrms": return sqrt(self._YY * Prob**(2/(self.N-3)))
        if self.norm == "chisq": return self._YY * self.wsum * Prob**(2/(self.N-3))
        p = 1 - Prob**(2/(self.N-3))
        if self.norm == "ZK": return p
        if self.norm == "lnL": return -0.5*self._YY*self.wsum*(1.-p) - 0.5*np.sum(np.log(2*np.pi * self.e_y**2))
        if self.norm == "dlnL": return 0.5 * self._YY * self.wsum * p


    def FAP(self, Pn):
        """
        Obtain the false-alarm probability (FAP).

        The FAP denotes the probability that at least one out of M independent
        power values in a prescribed search band of a power spectrum computed
        from a white-noise time series is as large as or larger than the
        threshold, `Pn`. It is assessed through

        .. math:: FAP(Pn) = 1 - (1-Prob(P>Pn))^M \\; ,

        where "Prob(P>Pn)" depends on the type of periodogram and normalization
        and is calculated by using the *prob* method; *M* is the number of
        independent power values and is computed internally.

        Parameters
        ----------
        Pn : float
            Power threshold.

        Returns
        -------
        FAP : float
            False alarm probability.

        """
        prob = self.M * self.prob(Pn)
        if prob > 0.01:
            return 1. - (1.-self.prob(Pn))**self.M
        return prob

    def powerLevel(self, FAPlevel):
        """
        Power threshold for FAP level.

        Parameters
        ----------
        FAPlevel : float or array_like
              "False Alarm Probability" threshold

        Returns
        -------
        Threshold : float or array
            The power threshold pertaining to a specified false-alarm
            probability (FAP). Powers exceeding this threshold have FAPs
            smaller than FAPlevel.

        """
        Prob = 1. - (1.-FAPlevel)**(1./self.M)
        return self.probInv(Prob)

    def stats(self, Pn):
        """
        Obtain basic statistics for power threshold.

        Parameters
        ----------
        Pn : float
            Power threshold.

        Returns
        -------
        Statistics : dictionary
            A dictionary containing {'Pn': *Pn*, 'Prob': *Prob(Pn)* ,
            'FAP': *FAP(Pn)*} for the specified power threshold, *Pn*.

        """
        return {'Pn': Pn, 'Prob': self.prob(Pn), 'FAP': self.FAP(Pn)}

    def toFile(self, ofile, header=True):
        """
        Write periodogram to file.

        The output file is a standard text file with two columns,
        viz., frequency and power. 

        Parameters
        ----------
        ofile : string
            Name of the output file.

        """
        with open(ofile, 'w') as f:
            if header:
                f.write("# Generalized Lomb-Scargle periodogram\n")
                f.write("# Parameters:\n")
                f.write("#    Data file: %s\n" % self.df)
                f.write("#    ofac     : %s\n" % self.ofac)
                f.write("#    norm     : %s\n" % self.norm)
                f.write("# 1) Frequency, 2) Normalized power\n")
            for line in zip(self.freq, self.power):
                f.write("%f  %f\n" % line)

        print("Results have been written to file: ", ofile)
