# Copyright 2016 Stefan Czesla, Tanja Molle

# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in the
# Software without restriction, including without limitation the rights to use, copy,
# modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so, subject to the
# following conditions:

# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE
# AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from __future__ import print_function, division
import numpy as np
import six.moves as smo
import scipy.special as ss


class SampCorr(object):

    def __init__(self):
        pass

    def getak(self, N):
        # Calculate the required coefficients (a_k)
        ak = np.zeros(N + 2)
        for k in smo.range(N + 2):
            ak[k] = (-1.0)**k * ss.binom(N + 1, k)
        return ak

    def get_rhok(self, N):
        ak = self.getak(N)
        rho = np.zeros(len(ak))
        a = np.concatenate((ak, np.zeros(len(ak))))
        b = np.concatenate((ak, np.zeros(len(ak))))
        r0 = np.sum(a * b)
        for k in smo.range(len(ak)):
            rho[k] = np.sum(a * b) / r0
            b = np.roll(b, 1)
        return rho

    def get_r(self, N, ss=None):
        rk = self.get_rhok(N)
        if ss is None:
            # Sample size not given
            r = 1.0 + 2.0 * np.sum(rk[1:])
        else:
            # Sample size given
            k = np.arange(len(rk))
            r = 1.0 + 2.0 * np.sum((ss - k[1:]) / ss * rk[1:])
        return r

    def get_nb(self, N, ss):
        """
        BH46 Eq. 14
        """
        rk = self.get_rhok(N)
        # Sample size given

        onb = 1. / ss

        k = np.arange(len(rk), dtype=float)
        #print("k: ", k[1:], (1.-k[1:]/ss), rk )
        onb += (2.0 / ss) * np.sum((1. - k[1:] / ss) * rk[1:])
        nb = 1. / onb
        return nb, onb

    def get_nvr(self, N):
        rho = self.get_rhok(N)
        rf = np.sum(rho**2)
        return rf


class BSBase(object):
    """
    Functionality to estimate noise from equidistantly and arbitrarily sampled data.
    """

    def __init__(self):
        # Factor to convert MAD into std estimate (0.67448975019608171)
        self._madfac = np.sqrt(2.) * ss.erfinv(2. * 0.75 - 1.)

    def _checkJP(self, j):
        """
        Check validity of specified jump parameter.
        """
        if j < 1:
            raise(ValueError(
                "Jump parameter must be an integer equal or larger than 1. Currently, a value of " + str(j) + " is specified."))

    def _checkN(self, N):
        """
        Check validity of specified order of approximation.
        """
        if N < 0:
            raise(ValueError(
                "Order of approximation (N) must be an integer larger than zero. Currently, a value of " + str(N) + " is specified."))

    def subsetIndexDissection(self, ndp, N, j):
        """
        Find array indices of (N+2)-length subsets.

        Here, (N+2)-sized sets of indices are constructed, which
        are required in the construction of the beta sample.

        Parameters
        ----------
        ndp : int
            Number of available data points
        N : int
            Last order of the Taylor series taken into account. Chunk
            length will be N+2.
        j : int
            Jump parameter (>0)

        Returns
        -------
        k-indices : list
            A list holding N+2 elements. Each element of the list is an
            array holding the indices of the data point no. one, two, ..., N+2
            of the subsets required to construct the beta sample.
        """
        self._checkJP(j)

        # Range of available indices
        x = np.arange(ndp, dtype=np.int)

        # Chunk length from which to extract subsamples
        step = (2 + N) * j

        # Temporarily store indices of subsets
        tmp = [list() for _ in smo.range(N + 2)]

        for k in smo.range(N + 2):
            for cj in smo.range(j):
                # Loop over jump parameter (cj = current jump)

                # Collect all data points being first, second, etc. in
                # individual chunks

                # The middle expression achieves that only complete subsets are considered.
                # The "or None" statement effectively makes 0 indicate the entire index range.
                # From point with no. k, the distance to the last point is l(k) = (N+1)*j+1 - k*j.
                s = x[k * j + cj:-((N + 1) * j + 1 - k * j) + 1 or None:step]
                tmp[k].extend(s.copy())

        # Convert into arrays
        result = [np.array(x) for x in tmp]

        return result

    def subsetIndices(self, cd):
        """
        Compose indices sets defining the subsets to construct the beta sample.

        Parameters
        ----------
        cd : list
            The output of subsetIndexDissection

        Returns
        -------
        Subset indices : 2d array
            A 2d-array holding the indices of all subsets, arranged
            so that result[i,::] holds the N+2 indices pertaining
            to the i-th subset.
        """
        # Length of individual subsets (N+2)
        cl = len(cd)
        # Number of chunks
        nc = len(cd[0])

        result = np.zeros((nc, cl), dtype=np.int)
        for i in smo.range(cl):
            result[::, i] = cd[i]
        return result

    def stdc4(self, n):
        """
        Calculate c4 factor.

        The c4 factor is required to obtain an unbiased estimator
        for the standard deviation.

        It is proportional to the factor B used in Kenney 1940, who
        started from a slightly different definition of the sample
        variance.

        Parameters
        ----------
        n : int
            Number of points in the sample

        Returns
        -------
        c4, ln(c4) : float
            The c4 factor and its natural logarithm.
        """
        lnc4 = 0.5 * (np.log(2.0) - np.log(n - 1.0)) + \
            ss.gammaln(n / 2.0) - ss.gammaln((n - 1.) / 2.)
        c4 = np.exp(lnc4)
        return c4, lnc4

    def stdUnbiased(self, y):
        """
        Get unbiased estimate of the standard deviation and its standard deviation.

        Parameters
        ----------
        y : array
            Sample from which to estimate standard deviation.

        Returns
        -------
        Std : float
            Unbiased estimate of the standard deviation.
        Std of std : float
            Standard deviation of the unbiased standard deviation estimator.
        """
        s = np.std(y, ddof=1)
        c4 = self.stdc4(len(y))[0]
        return s / c4, s / c4 * np.sqrt(c4**(-2) - 1.0)

    def meanUnbStd(self, x):
        """
        Mean and (unbiased) standard deviation of the mean.

        Parameters
        ----------
        x : array
            Sampe from which to calculate mean and std

        Returns
        -------
        mean, std : floats
            Mean and unbiased standard deviation of the mean
        """
        m = np.mean(x)
        s = self.stdUnbiased(x)[0] / np.sqrt(len(x))
        return m, s

    def estimateStdMAD(self, x, mode):
        """
        Estimate standard deviation based on median absolute deviation (MAD)

        Parameters
        ----------
        x : array
            Sample from which to determine the standard deviation.
        mode : string, {zm, em}
            If 'zm', the population median is assumed to be zero.
            If 'em', the population median is estimated as the
            sample median.

        Returns
        -------
        std : float
            Estimate of the standard deviation
        """
        if mode == "zm":
            mad = np.median(np.abs(x))
        elif mode == "em":
            mad = np.median(np.abs(x - np.median(x)))
        sigest = mad / self._madfac
        return sigest

    def variance(self, x, mode):
        """
        Estimate variance from sample

        Parameters
        ----------
        x : array
            sample
        mode : string, {n, nmo}
            Estimator to use: 'n' for 1/n version with
            zero mean (not estimated) and 'nmo' for 1/(n-1)
        """
        if mode == "n":
            # 1/n version
            return np.sum(x**2) / len(x)
        elif mode == "nmo":
            # 1/(n-1) version
            m = np.mean(x)
            return np.sum((x - m)**2) / (len(x) - 1.)
        else:
            raise(ValueError("Unknown mode for variance: " + str(mode)))


class BSArbSamp(BSBase):
    """
    Estimate noise in equidistantly sampled data.
    """

    def __init__(self):
        super(BSArbSamp, self).__init__()
        self.betaSample = None
        self.betaSig, self.stdBetaSig = None, None
        self.mbeta, self.stdmbeta = None, None

    def getCoeffsArbSamp(self, t, gamma=1.):
        """
        Calculate coefficients (ak) for arbitrary sampling.

        Parameters
        ----------
        t : array
            Sampling instants of the subset.
        gamma : float, optional
            Scaling of the coefficients (default is one).

        Returns
        -------
        ak : array
            Set of coefficients.
        """
        np2 = len(t)
        Tp = np.matrix(np.zeros((np2 - 1, np2 - 1)))
        ti = np.zeros(np2 - 1)
        for i in smo.range(np2 - 1):
            Tp[i, ::] = t[1:]**i
            ti[i] = t[0]**i
        ap = -np.array(np.dot(Tp.I, np.matrix(ti).T)).ravel()
        a = np.ones(np2)
        a[1:] = ap
        return gamma * a

    def getBetaSampleShift(self, x, y, N, j):
        """
        Construct beta sample using shifting procedure.

        Parameters
        ----------  
        x, y : array
            Values from which to estimate noise
        N : int
            Last order of the Taylor series taken into account
        j : int
            Jump parameter (>0)

        Returns
        -------
        Beta sample : array
            An array holding all available beta values.
        """

        # Required chunk length to obtain realization
        # of beta = (N+2) + (N-1)*(j-1) = (N+1)*j + 1
        cl = (N + 1) * j + 1

        b = np.zeros(len(y) - cl + 1)

        for i in smo.range(len(y) - cl + 1):
            # Realizations of the beta sample
            ak = self.getCoeffsArbSamp(x[i:i + cl:j])
            b[i] = np.sum(y[i:i + cl:j] * ak)
            # Normalize to obtain width of sigma_0
            b[i] /= np.sqrt(np.sum(ak**2))

        return b

    def getBetaSample(self, x, y, N, j):
        """
        Combine data points to calculate beta values.

        Parameters
        ----------
        x : array
            Sampling of data
        y : array
            Values from which to estimate noise
        N : int
            Last order of the Taylor series taken into account
        j : int
            Jump parameter (>0)

        Returns
        -------
        Betas : array
            An array holding all available beta values.
        """
        self._checkJP(j)
        self._checkN(N)

        self.karr = self.subsetIndexDissection(len(y), N, j)
        # Indices of subsets required to construct beta sample
        ssindices = self.subsetIndices(self.karr)

        result = np.zeros(len(ssindices))

        for c in smo.range(len(ssindices)):
            chunkSamp = x[ssindices[c, ::]].copy()
            # Subtract first value for numerical reasons
            chunkSamp -= chunkSamp[0]
            ak = self.getCoeffsArbSamp(chunkSamp)

            for k in smo.range(N + 2):
                result[c] += ak[k] * y[ssindices[c, k]]

            # Normalize to obtain width of sigma_0
            result[c] /= np.sqrt(np.sum(ak**2))
        return result

    def betaSigma(self, x, y, N, j, ignoreNaN=True, returnMAD=False, ibs=False):
        """
        Estimate standard deviation of noise term in data set y.

        In this implementation, arbitrary sampling is taken into account.

        The method assigns the following attributes, which may be accessed after execution
        to work with the result:

        ============  ========  ===========================================
        Attribute     Type      Meaning
        ------------  --------  -------------------------------------------
        betaSample    array     The beta sample

        estimates     dict      Summary of the estimates obtained
                                from the beta sample (bs).
                                
                                "s2E": variance estimate of the bs
                                (expectation value of zero),
                                "s2Evar": variance of s2E
                                
                                "sE": Estimate of std of bs,
                                "sEstd": Std of sE
                                
                                "s2": variance estimate of the bs
                                (mean estimated from bs),
                                "s2var": variance of s2,
                                "s": sqrt(s2)
                                
                                "sME": Standard deviation based on MAD
                                with zero expectation,
                                "sMEstd": Estimation of std of sME
                                
                                "sM": Standard deviation based on MAD
                                with median estimated from bs,
                                "sMEstd": Estimation of std of sM
        ============  ========  ===========================================

        Parameters
        ----------   
        x : array
            Sampling of the data.
        y : array
            Data values from which to estimate standard deviation of noise.
        N : int
            Last order of the Taylor series to be taken into account.
        j : int, optional
            Jump parameter (default is one, i.e., consecutive data points are combined
            to estimate the noise).
        ignoreNaN : boolean, optional
            If True (default), NaN values in the beta sample are ignored in the
            calculation.
        returnMAD : boolean, optional
            If True, the estimate obtained using the MAD is returned instead of that
            of the MVU estimator (default is False). The standard error is estimated by
            scaling the standard deviation of the MVU estimator by a factor of
            1.64.
        ibs : boolean, optional
            If True, an independent beta sample is constructed. Default is False.

        Returns
        -------
        Estimate of STD in beta sample and the STD of the estimate: float, float
            The standard deviation determined in the beta sample. Non-robust
            estimates sE and sEstd if `returnMAD` is False (default) or robust
            estimates sME and sMEstd if `returnMAD` is True. The `estimates`
            attribute holds a more comprehensive summary of the estimates.
        """
        self._checkJP(j)
        self._checkN(N)

        if ibs:
            # Independent beta sample
            self.betaSample = self.getBetaSample(x, y, N, j)
        else:
            # Correlated beta sample
            self.betaSample = self.getBetaSampleShift(x, y, N, j)

        if len(self.betaSample) < 2:
            raise(ValueError(
                "Size of beta sample is too small (" + str(len(self.betaSample)) + ")"))

        if ignoreNaN:
            # Consider only non-NaN values
            bs = self.betaSample[np.isnan(self.betaSample) == False]
        else:
            # Use entire beta sample
            bs = self.betaSample

        # Get variance of beta sample
        self.bvar = self.variance(bs, mode="n")
        if ibs:
            # Independent beta sample
            self.bvarvar = 2. * self.bvar**2 / len(bs)
        else:
            # Dependent beta sample with correlation
            # Use estimate from independent sample to obtain estimate
            # because correlation structure is not exactly known
            self.bvarvar = 2. * self.bvar**2 / len(x) * (N + 2)

        # Estimate STD and its STD by propagation
        self.bstd = np.sqrt(self.bvar)
        self.bstdstd = np.sqrt(self.bvarvar / (4. * self.bvar))

        # Get standard deviation of beta sample and its standard error
        #self.betaSig, self.stdBetaSig = self.stdUnbiased(bs)
        # Get mean value of beta sample and its standard deviation
        self.mbeta, self.stdmbeta = self.meanUnbStd(bs)
        # Estimate the standard deviation based on MAD
        self.madstd = self.estimateStdMAD(bs, mode='zm')

        if ibs:
            self.stdMadstd = np.sqrt(
                1. / (2. * len(bs))) * self.estimateStdMAD(bs, mode='zm') * 1.64
        else:
            # Dependent beta sample with correlation
            # Use estimate assuming independent sample (for size of beta sample)
            # to obtain estimate of correlation
            lbs = len(x) / (N + 2)
            self.stdMadstd = np.sqrt(1. / (2. * lbs)) * \
                self.estimateStdMAD(bs, mode='zm') * 1.64

        # Estimates and variances
        # svar : Assumed to be equal to sEvar (more accurate formula given
        #        by Bayley and Hammersley 1946)
        self.estimates = {"s2E": self.bvar, "s2Evar": self.bvarvar,
                          "sE": self.bstd, "sEstd": self.bstdstd,
                          "s2": self.variance(bs, mode="nmo"),
                          "s2var": self.bstd,
                          "s": np.sqrt(self.variance(bs, mode="nmo")),
                          "sME": self.madstd, "sMEstd": self.stdMadstd,
                          "sM": self.estimateStdMAD(bs, mode='em'),
                          "sMEstd": self.stdMadstd}

        if not returnMAD:
            return self.bstd, self.bstdstd
        else:
            return self.madstd, self.stdMadstd


class BSEqSamp(BSBase):

    def __init__(self):
        super(BSEqSamp, self).__init__()
        self.betaSample = None
        self.betaSig, self.stdBetaSig = None, None
        self.mbeta, self.stdmbeta = None, None

    def get_ak(self, N):
        """
        Calculate the required coefficients (a_k)

        Parameters
        ----------
        N : int
            Order of approximation

        Returns
        -------
        ak : array
            The coefficients
        """
        ak = np.zeros(N + 2)
        for k in smo.range(N + 2):
            ak[k] = (-1.0)**k * ss.binom(N + 1, k)
        return ak

    def get_rhok(self, N):
        """
        Calculate (auto)correlation function for beta sample

        Parameters
        ----------
        N : int
            Order of approximation

        Returns
        -------
        correlation function : array
            Correlation function
        """
        ak = self.get_ak(N)
        rho = np.zeros(len(ak))
        a = np.concatenate((ak, np.zeros(len(ak))))
        b = np.concatenate((ak, np.zeros(len(ak))))
        r0 = np.sum(a * b)
        for k in smo.range(len(ak)):
            rho[k] = np.sum(a * b) / r0
            b = np.roll(b, 1)
        return rho

    def getBetaSample(self, y, N, j):
        """
        Construct the beta sample.

        Parameters
        ----------  
        y : array
            Values from which to estimate noise
        N : int
            Last order of the Taylor series taken into account
        j : int
            Jump parameter (>0)

        Returns
        -------
        Beta sample : array
            An array holding all available beta values.
        """
        self.karr = self.subsetIndexDissection(len(y), N, j)

        # Calculate the required coefficients (a_k)
        self.ak = self.get_ak(N)

        result = np.zeros(self.karr[0].size)

        for k in smo.range(N + 2):
            result += self.ak[k] * y[self.karr[k]]

        # Normalize to obtain width of sigma_0
        result /= np.sqrt(ss.binom(2 * N + 2, N + 1))
        return result

    def getBetaSampleShift(self, y, N, j):
        """
        Construct beta sample using shifting procedure.

        Parameters
        ----------  
        y : array
            Values from which to estimate noise
        N : int
            Last order of the Taylor series taken into account
        j : int
            Jump parameter (>0)

        Returns
        -------
        Beta sample : array
            An array holding all available beta values.
        """
        # Calculate the required coefficients (a_k)
        self.ak = self.get_ak(N)

        # Required chunk length to obtain realization
        # of beta = (N+2) + (N-1)*(j-1) = (N+1)*j + 1
        cl = (N + 1) * j + 1

        b = np.zeros(len(y) - cl + 1)

        for i in smo.range(len(y) - cl + 1):
            # Realizations of the beta sample
            b[i] = np.sum(y[i:i + cl:j] * self.ak)

        # Normalize to obtain width of sigma_0
        b /= np.sqrt(ss.binom(2 * N + 2, N + 1))
        return b

    def betaSigma(self, y, N, j, ignoreNaN=True, returnMAD=False, ibs=False):
        """
        Estimate standard deviation of noise term in data set y.

        It is explicitly assumed that the data are equidistantly sampled.

        ============  ========  ===========================================
        Attribute     Type      Meaning
        ------------  --------  -------------------------------------------
        betaSample    array     The beta sample

        estimates     dict      Summary of the estimates obtained
                                from the beta sample (bs).
                                
                                "s2E": variance estimate of the bs
                                (expectation value of zero),
                                "s2Evar": variance of s2E
                                
                                "sE": Estimate of std of bs,
                                "sEstd": Std of sE
                                
                                "s2": variance estimate of the bs
                                (mean estimated from bs),
                                "s2var": variance of s2,
                                "s": sqrt(s2)
                                
                                "sME": Standard deviation based on MAD
                                with zero expectation,
                                "sMEstd": Estimation of std of sME
                                
                                "sM": Standard deviation based on MAD
                                with median estimated from bs,
                                "sMEstd": Estimation of std of sM
        ============  ========  ===========================================

        Parameters
        ----------   
        y : array
            Data values from which to estimate standard deviation of noise.
        N : int
            Last order of the Taylor series to be taken into account.
        j : int
            Jump parameter
        ignoreNaN : boolean, optional
            If True (default), NaN values in the beta sample are ignored in the
            calculation.
        returnMAD : boolean, optional
            If True, the estimate obtained using the MAD is returned instead of that
            of the MVU estimator (default is False). The standard error is estimated by
            scaling the standard deviation of the MVU estimator by a factor of
            1.64.
        ibs : boolean, optional
            If True, an independent beta sample is constructed. Default is False.

        Returns
        -------
        Estimate of STD in beta sample and the STD of the estimate: float, float
            The standard deviation determined in the beta sample. Non-robust
            estimates sE and sEstd if `returnMAD` is False (default) or robust
            estimates sME and sMEstd if `returnMAD` is True. The `estimates`
            attribute holds a more comprehensive summary of the estimates.
        """
        self._checkJP(j)
        self._checkN(N)

        if ibs:
            # Get independent beta sample
            self.betaSample = self.getBetaSample(y, N, j)
        else:
            # Get correlated sample from shifting procedure
            self.betaSample = self.getBetaSampleShift(y, N, j)

        if len(self.betaSample) < 2:
            raise(ValueError(
                "Size of beta sample is too small (" + str(len(self.betaSample)) + ")"))

        if ignoreNaN:
            # Consider only non-NaN values
            bs = self.betaSample[np.isnan(self.betaSample) == False]
        else:
            # Use entire beta sample
            bs = self.betaSample

        # Get variance of beta sample
        self.bvar = self.variance(bs, mode="n")
        if ibs:
            # Independent beta sample
            self.bvarvar = 2. * self.bvar**2 / len(bs)
        else:
            # Dependent beta sample with correlation
            rhok = self.get_rhok(N)
            self.bvarvar = 2. * self.bvar**2 / \
                len(bs) * (1. + 2. * np.sum(rhok[1:]**2))

        # Estimate STD and its STD by propagation
        self.bstd = np.sqrt(self.bvar)
        self.bstdstd = np.sqrt(self.bvarvar / (4. * self.bvar))

        # Get standard deviation of beta sample and its standard error
        #self.betaSig, self.stdBetaSig = self.stdUnbiased(bs)
        # Get mean value of beta sample and its standard deviation
        self.mbeta, self.stdmbeta = self.meanUnbStd(bs)
        # Estimate the standard deviation based on MAD
        self.madstd = self.estimateStdMAD(bs, mode='zm')

        if ibs:
            rfac = 1.0
        else:
            # Dependent beta sample with correlation
            rhok = self.get_rhok(N)
            rfac = (1. + 2. * np.sum(rhok[1:]**2))
            # Estimate the STD of the estimate
        self.stdMadstd = np.sqrt(rfac / (2. * len(bs))) * \
            self.estimateStdMAD(bs, mode='zm') * 1.64

        # Estimates and variances
        # svar : Assumed to be equal to sEvar (more accurate formula given
        #        by Bayley and Hammersley 1946)
        self.estimates = {"s2E": self.bvar, "s2Evar": self.bvarvar,
                          "sE": self.bstd, "sEstd": self.bstdstd,
                          "s2": self.variance(bs, mode="nmo"),
                          "s2var": self.bstd,
                          "s": np.sqrt(self.variance(bs, mode="nmo")),
                          "sME": self.madstd, "sMEstd": self.stdMadstd,
                          "sM": self.estimateStdMAD(bs, mode='em'),
                          "sMEstd": self.stdMadstd}

        if not returnMAD:
            return self.bstd, self.bstdstd
        else:
            return self.madstd, self.stdMadstd

