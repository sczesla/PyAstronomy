from __future__ import print_function, division
import unittest
from PyAstronomy import funcFit as fuf
import os
import numpy as np


class TestFuncFitSanity(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        if os.path.isfile("saveState.tmp"):
            os.remove("saveState.tmp")

    def testsanity_saveState(self):
        gf = fuf.GaussFit1d()
        gf.assignValue({"A": 1, "mu": 2, "off": 3, "sig": 4, "lin": 5.0})
        gf.setRestriction({"A": [0, 10], "off": [None, 7.8]})
        gf.thaw(["A", "mu", "lin"])
        dat = gf.saveState("saveState.tmp", clobber=True)

        gf2 = fuf.GaussFit1d()
        gf2.restoreState("saveState.tmp")
        self.assertEqual(gf2.parameters(), gf.parameters())
        self.assertEqual(gf2.getRestrictions(), gf.getRestrictions())
        self.assertEqual(gf2.frozenParameters(), gf.frozenParameters())

        gf2.restoreState(dat)
        self.assertEqual(gf2.parameters(), gf.parameters())
        self.assertEqual(gf2.getRestrictions(), gf.getRestrictions())
        self.assertEqual(gf2.frozenParameters(), gf.frozenParameters())

    def testsanity_coordinateGrid(self):
        x = np.arange(5)
        y = np.arange(10) + 10.0
        z = np.arange(2) + 100.0
        g = fuf.coordinateGrid(x, y, z)
        self.assertEqual(g[0, 0, 0, 0], 0.0)
        self.assertEqual(g[0, 0, 0, 1], 10.0)
        self.assertEqual(g[0, 0, 0, 2], 100.0)
        self.assertEqual(g[1, 0, 1, 0], 1.0)
        self.assertEqual(g[1, 2, 1, 1], 12.0)
        self.assertEqual(g[1, 2, 1, 2], 101.0)

    def testsanity_fminArgs(self):
        """
        Test handling of fminArgs
        """
        np.random.seed(101)

        fwhm = 10833/100_000
        snr = 100
        
        # Wvl and 'data' (independent, without signal)
        w = np.arange(10830, 10836, fwhm)
        d = np.random.normal(1, 1/snr, len(w))
        e = np.ones_like(w)/snr
        
        gf = fuf.GaussFit1d()
        gf["mu"] = 10833
        gf["off"] = 1.0
        gf["A"] = -0.1
        # Assumption
        gf["sig"] = fwhm*3
        
        # 'Real' estimate from the data
        gf.thaw(["A"])
        gf.fit(w, d, yerr=e, minAlgo='spfmp', fminArgs={'disp':0})
