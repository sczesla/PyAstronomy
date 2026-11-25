from __future__ import print_function, division
import unittest


class TestStringlengthSanity(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testsanity_slExample(self):
        """Checking stringlength example"""
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy import pyTiming as pyt

        period = 1.75

        x = np.linspace(0, 10, 100)
        y = 20 * np.sin(2 * np.pi * x / period)
        y += np.random.normal(0, 1, len(x))

        # Trial periods to be tested (200 trial periods between 0.5 and 4.5;
        # same units as x-axis)
        tps = (0.5, 4.5, 200)

        # Calculate string length
        p, sl = pyt.stringlength_dat(x, y, tps)

        # Show the string length. An alias at the double period
        # (half frequency) is obvious.

    #         plt.plot(p, sl, 'b.-')
    #         plt.ylabel("String length")
    #         plt.xlabel("Trial period")
    #         plt.show()

    def testsanity_stringlength(self):
        """Checking sanity of stringlength"""

        import numpy as np
        from PyAstronomy import pyTiming as pyt

        x = np.array([0.0, 0.5, 0.999999])
        y = np.array([0.0, 0.5, 1.0])

        # Calculate string length
        p, sl = pyt.stringlength_dat(x, y, tps=np.array([1]), norm="no")
        self.assertAlmostEqual(
            sl[0],
            np.sqrt(2.0) + 1.0,
            delta=1e-6,
            msg="closed stringlength does not match",
        )

        p, sl = pyt.stringlength_dat(x, y, tps=np.array([1]), norm="no", closed=False)
        self.assertAlmostEqual(
            sl[0],
            np.sqrt(2.0),
            delta=1e-6,
            msg="stringlength (not closed) does not match",
        )


class TestGlsSanity(unittest.TestCase):
    def testsanity_fastMode(self):
        """Checking that GLS fast mode reproduces the direct computation"""
        import numpy as np
        from PyAstronomy.pyTiming import pyPeriod

        rng = np.random.default_rng(1)
        t = np.sort(100 * rng.random(80))
        e = 0.2 + 0.3 * rng.random(80)
        y = 1.5 * np.sin(2 * np.pi * t / 7.3 + 0.4) + e * rng.normal(size=80)

        for ls in (False, True):
            slow = pyPeriod.Gls((t, y, e), ls=ls)
            fast = pyPeriod.Gls((t, y, e), ls=ls, fast=True)
            # Fast mode used to shift the power by one frequency step
            np.testing.assert_allclose(fast.freq, slow.freq)
            np.testing.assert_allclose(fast.power, slow.power, atol=1e-9)
            self.assertEqual(fast.hpstat["fbest"], slow.hpstat["fbest"])

    def testsanity_lnL(self):
        """Checking the lnL normalization against the direct log-likelihood"""
        import numpy as np
        from PyAstronomy.pyTiming import pyPeriod

        rng = np.random.default_rng(2)
        t = np.sort(100 * rng.random(60))
        e = 0.2 + 0.3 * rng.random(60)
        y = 3 + 1.5 * np.sin(2 * np.pi * t / 7.3 + 0.4) + e * rng.normal(size=60)

        for err in (e, None):
            lc = (t, y) if err is None else (t, y, err)
            sig = np.ones_like(t) if err is None else err
            gls = pyPeriod.Gls(lc, norm="lnL")
            # Gaussian log-likelihood of the best-fit sine
            res = y - gls.sinmod(t)
            lnL = -0.5 * np.sum(res**2 / sig**2) - 0.5 * np.sum(np.log(2 * np.pi * sig**2))
            self.assertAlmostEqual(gls.power.max(), lnL, delta=1e-8 * abs(lnL))
            # prob and probInv are inverse to each other
            for prob in (1e-3, 0.1):
                self.assertAlmostEqual(gls.prob(gls.probInv(prob)), prob, delta=1e-10)

    def testsanity_offset(self):
        """Checking that the best-fit sine reproduces the chi^2 of the power"""
        import numpy as np
        from PyAstronomy.pyTiming import pyPeriod

        rng = np.random.default_rng(3)
        t = np.sort(100 * rng.random(60))
        e = 0.2 + 0.3 * rng.random(60)
        y = 3 + 1.5 * np.sin(2 * np.pi * t / 7.3 + 0.4) + e * rng.normal(size=60)

        for ls in (False, True):
            gls = pyPeriod.Gls((t, y, e), ls=ls)
            # chi^2 of the best-fit sine from the power ...
            chi2p = gls._YY * gls.wsum * (1 - gls.p.max())
            # ... and directly from the model, including its offset
            chi2 = np.sum((y - gls.sinmod(t))**2 / e**2)
            self.assertAlmostEqual(chi2, chi2p, delta=1e-9 * chi2p)
            if ls:
                # Without fitted offset, the offset is the weighted mean
                w = 1 / e**2
                self.assertAlmostEqual(gls.hpstat["offset"], np.sum(w * y) / np.sum(w), delta=1e-12)

    def testsanity_info(self):
        """Checking that info returns numbers"""
        import numpy as np
        from PyAstronomy.pyTiming import pyPeriod

        rng = np.random.default_rng(4)
        t = np.sort(100 * rng.random(40))
        e = 0.2 + 0.3 * rng.random(40)
        y = np.sin(2 * np.pi * t / 7.3) + e * rng.normal(size=40)

        info = pyPeriod.Gls((t, y, e)).info(noprint=True)
        mwie = info["mean_weighted_internal_error"]
        self.assertIsInstance(mwie, float)
        self.assertAlmostEqual(mwie, np.sqrt(len(e) / np.sum(1 / e**2)), delta=1e-12)

    def testsanity_infoMaximumPower(self):
        """Checking that info reports the power of the highest peak"""
        import numpy as np
        from PyAstronomy.pyTiming import pyPeriod

        rng = np.random.default_rng(11)
        t = np.sort(100 * rng.random(40))
        e = 0.2 + 0.3 * rng.random(40)
        y = np.sin(2 * np.pi * t / 7.3) + e * rng.normal(size=40)

        for norm in pyPeriod.Gls.norms:
            if norm == "Scargle":
                continue
            gls = pyPeriod.Gls((t, y, e), norm=norm)
            k = np.argmax(gls.p)
            self.assertEqual(gls.info(noprint=True)["maximum_power"], gls.power[k])
        # In chisq (and wrms), the highest peak is the smallest value
        gls.pnorm("chisq")
        self.assertEqual(gls.info(noprint=True)["maximum_power"], gls.power.min())

    def testsanity_FAP(self):
        """Checking FAP and powerLevel, also for small probabilities"""
        import numpy as np
        from PyAstronomy.pyTiming import pyPeriod

        rng = np.random.default_rng(5)
        t = np.sort(100 * rng.random(50))
        e = 0.2 + 0.3 * rng.random(50)
        y = np.sin(2 * np.pi * t / 7.3) + e * rng.normal(size=50)
        gls = pyPeriod.Gls((t, y, e))

        # powerLevel is the inverse of FAP, also far below 1e-16
        for fap in (1e-20, 1e-8, 0.009, 0.011, 0.5):
            self.assertAlmostEqual(gls.FAP(gls.powerLevel(fap)) / fap, 1, delta=1e-9)

        # Exact value 1 - (1-q)^M = Mq - M(M-1)/2 q^2 + ..., not only Mq
        q = 0.005 / gls.M
        Pn = gls.probInv(q)
        expected = gls.M * q - gls.M * (gls.M - 1) / 2 * q**2 \
            + gls.M * (gls.M - 1) * (gls.M - 2) / 6 * q**3
        self.assertAlmostEqual(gls.FAP(Pn) / expected, 1, delta=1e-6)

    def testsanity_ScargleDeprecated(self):
        """Checking that the 'Scargle' normalization warns"""
        import numpy as np
        from unittest import mock
        from PyAstronomy.pyTiming import pyPeriod

        rng = np.random.default_rng(6)
        t = np.sort(100 * rng.random(30))
        y = rng.normal(size=30)

        # The normalization asks for the population variance via input()
        with mock.patch("builtins.input", return_value="1.0"):
            with self.assertWarns(FutureWarning):
                pyPeriod.Gls((t, y), norm="Scargle")
            gls = pyPeriod.Gls((t, y))
            with self.assertWarns(FutureWarning):
                gls.pnorm("Scargle")

    def testsanity_LombScargleVariance(self):
        """Checking that LombScargle is normalized by the sample variance"""
        import numpy as np
        from PyAstronomy.pyTiming import pyPeriod

        rng = np.random.default_rng(7)
        t = np.sort(100 * rng.random(50))
        y = np.sin(2 * np.pi * t / 7.3) + rng.normal(size=50)

        ls = pyPeriod.LombScargle(pyPeriod.TimeSeries(t, y), ofac=10, hifac=1)
        # Same periodogram in the normalization of Horne & Baliunas
        gls = pyPeriod.Gls((t, y), freq=ls.freq, ls=True, norm="HorneBaliunas")
        np.testing.assert_allclose(ls.power, gls.power, rtol=1e-8)
        # Thus, also the same probabilities
        Pn = ls.power.max()
        self.assertAlmostEqual(ls.prob(Pn) / gls.prob(Pn), 1, delta=1e-8)

    def testsanity_toFileDeprecated(self):
        """Checking that toFile warns and writes the periodogram"""
        import os
        import tempfile
        import numpy as np
        from PyAstronomy.pyTiming import pyPeriod

        rng = np.random.default_rng(10)
        t = np.sort(100 * rng.random(30))
        y = np.sin(2 * np.pi * t / 7.3) + rng.normal(size=30)
        gls = pyPeriod.Gls((t, y))

        with tempfile.TemporaryDirectory() as d:
            fn = os.path.join(d, "gls.txt")
            # Default call (with header) used to fail
            with self.assertWarns(FutureWarning):
                gls.toFile(fn)
            dat = np.loadtxt(fn)
        self.assertEqual(dat.shape, (gls.nf, 2))
        np.testing.assert_allclose(dat[:, 0], gls.freq, atol=1e-6)
        np.testing.assert_allclose(dat[:, 1], gls.power, atol=1e-6)


class TestFourierSanity(unittest.TestCase):
    def testsanity_probabilities(self):
        """Checking prob, probInv, FAP, and powerLevel of Fourier"""
        import numpy as np
        from PyAstronomy.pyTiming import pyPeriod

        rng = np.random.default_rng(8)
        lc = pyPeriod.TimeSeries(np.arange(1000.), rng.poisson(100, 1000))
        fft = pyPeriod.Fourier(lc)

        # Chi^2 distribution with 2 d.o.f.
        for Pn in (0.5, 2., 10.):
            self.assertAlmostEqual(fft.prob(Pn), np.exp(-Pn / 2), delta=1e-14)
            self.assertAlmostEqual(fft.Prob(Pn), fft.prob(Pn), delta=1e-14)
        for prob in (1e-6, 0.1, 0.5):
            self.assertAlmostEqual(fft.prob(fft.probInv(prob)) / prob, 1, delta=1e-10)
        for fap in (0.001, 0.05):
            self.assertAlmostEqual(fft.FAP(fft.powerLevel(fap)) / fap, 1, delta=1e-10)
        # Pure Poisson noise: fraction of powers above Pn is about exp(-Pn/2)
        self.assertAlmostEqual(np.mean(fft.power > 4), np.exp(-2), delta=0.04)

    def testsanity_input(self):
        """Checking that Fourier rejects uneven sampling and unknown norms"""
        import numpy as np
        from PyAstronomy.pyTiming import pyPeriod
        from PyAstronomy.pyaC import pyaErrors as PE

        y = np.ones(10)
        t = np.arange(10.)
        t[5] += 0.5
        with self.assertRaises(PE.PyAValError):
            pyPeriod.Fourier(pyPeriod.TimeSeries(t, y))
        with self.assertRaises(PE.PyAValError):
            pyPeriod.Fourier(pyPeriod.TimeSeries(np.arange(10.), y), norm="other")


class TestPeriodBaseSanity(unittest.TestCase):
    def testsanity_FAP(self):
        """Checking FAP and powerLevel of PeriodBase (via LombScargle)"""
        import numpy as np
        from PyAstronomy.pyTiming import pyPeriod

        rng = np.random.default_rng(9)
        t = np.sort(100 * rng.random(50))
        y = np.sin(2 * np.pi * t / 7.3) + rng.normal(size=50)
        ls = pyPeriod.LombScargle(pyPeriod.TimeSeries(t, y), ofac=10, hifac=1)

        # powerLevel is the inverse of FAP, also for arrays and far below 1e-16
        faps = np.array([1e-20, 1e-8, 0.009, 0.011, 0.5])
        np.testing.assert_allclose(ls.FAP(ls.powerLevel(faps)), faps, rtol=1e-9)

        # Exact value 1 - (1-q)^M = Mq - M(M-1)/2 q^2 + ..., not only Mq
        q = 0.005 / ls.M
        Pn = ls.probInv(q)
        expected = ls.M * q - ls.M * (ls.M - 1) / 2 * q**2 \
            + ls.M * (ls.M - 1) * (ls.M - 2) / 6 * q**3
        self.assertAlmostEqual(ls.FAP(Pn) / expected, 1, delta=1e-6)
