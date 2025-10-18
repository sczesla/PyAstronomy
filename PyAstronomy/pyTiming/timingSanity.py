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
