from __future__ import print_function, division
import unittest
import numpy
from PyAstronomy import funcFit as fuf
import six


class TestFuncFitSanity(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testsanity_parameterAssignment1(self):
        gf = fuf.GaussFit1d()
        origVars = ["A", "mu", "sig", "off", "lin"]
        vals = [1.0, 2.0, 3.0, 4.0, 5.0]
        for k, v in zip(origVars, vals):
            gf[k] = v
        for k, v in zip(origVars, vals):
            self.assertEqual(v, gf[k])
        gf.assignValue(dict(zip(origVars, numpy.zeros(len(origVars)))))
        self.assertEqual(numpy.sum(list(gf.parameters().values())), 0.0)

    def testsanity_combine1(self):
        gf = fuf.GaussFit1d()
        gff = gf + gf + gf
        for p in six.iterkeys(gff.parameters()):
            gff[p] = numpy.random.random()
            gff.thaw(p)
            gff.thaw([p, p])
            gff.freeze(p)
            gff.freeze([p, p])
            gff.setRestriction({p: [None, None]})
        for prop in ["A", "mu", "sig", "off", "lin"]:
            for c in [1, 2, 3]:
                gff[prop, "Gaussian", c] = numpy.random.random()
                s = (prop, "Gaussian", c)
                gff.thaw(s)
                gff.thaw([s, s])
                gff.freeze(s)
                gff.freeze([s, s])
                gff.setRestriction({s: [None, 10.0]})

    def testsanity_description(self):
        """
        Check sanity of 'description'
        """
        gf = fuf.GaussFit1d()
        gff = gf + gf
        self.assertEqual(
            gff.description(),
            "Gaussian(No. 1) + Gaussian(No. 2)",
            "Wrong description: " + gff.description(),
        )

    def testsanity_gaussfit1d_parameterization(self):
        """ Test parameterization of GaussFit1d """
        import numpy as np
        
        gf_hs = fuf.GaussFit1d(prm=("h", "sig"))
        gf_as = fuf.GaussFit1d(prm=("A", "sig"))
        gf_af = fuf.GaussFit1d(prm=("A", "FWHM"))
        gf_hf = fuf.GaussFit1d(prm=("h", "FWHM"))
        
        fsig = 2*np.sqrt(2*np.log(2))
        print(f"fsig = {fsig}")
        
        gfs = [gf_hs, gf_as, gf_af, gf_hf]
        
        gf_as["A"] = 1
        gf_as["sig"] = 1
        
        gf_af["A"] = 1
        gf_af["FWHM"] = fsig
        
        gf_hs["h"] = 1/np.sqrt(2*np.pi)
        gf_hs["sig"] = 1
        
        gf_hf["h"] = 1/np.sqrt(2*np.pi)
        gf_hf["FWHM"] = fsig
        
        x = np.linspace(-3,3,100)
        ys = []
        for g in gfs:
            ys.append( g.evaluate(x) )
            
        for i in range(1, len(ys)):
            print(f"Checking y({i}) against 0")
            np.testing.assert_allclose(ys[i], ys[0], rtol=0, atol=1e-8)

    def testsanity_gaussfit1d_parameterization_fit(self):
        """ Test parameterization of GaussFit1d (fitting) """
        import numpy as np
        
        np.random.seed(1919)
        
        gf_hs = fuf.GaussFit1d(prm=("h", "sig"))
        gf_as = fuf.GaussFit1d(prm=("A", "sig"))
        gf_af = fuf.GaussFit1d(prm=("A", "FWHM"))
        gf_hf = fuf.GaussFit1d(prm=("h", "FWHM"))
        
        fsig = 2*np.sqrt(2*np.log(2))
        print(f"fsig = {fsig}")
        
        gfs = [gf_hs, gf_as, gf_af, gf_hf]
        
        gf_as["A"] = 1
        gf_as["sig"] = 1
        
        x = np.linspace(-3,3,100)
        y = gf_as.evaluate(x)
        y += np.random.normal(0, 0.001, len(x))
        
        gf_hf.thaw("h")
        gf_hf.thaw("FWHM")
        gf_hf.thaw("mu")
        
        gf_hf["h"] = 0.1
        gf_hf["FWHM"] = 2
        
        gf_hf.fit(x, y, minAlgo="spfmp")
        gf_hf.parameterSummary()
        
        np.testing.assert_almost_equal(gf_hf["h"], 1/np.sqrt(2*np.pi), decimal=3)
        np.testing.assert_almost_equal(gf_hf["FWHM"], fsig, decimal=3)
        np.testing.assert_almost_equal(gf_hf["mu"], 0, decimal=3)


    def testsanity_modeGrenander(self):
        """ Test sanity of mode Grenander """
        import numpy as np
        from PyAstronomy import funcFit as fuf
        
        np.random.seed(1919)
        x = np.random.normal(0.9, 2, 10000)
        
        k = 20
        for p in [1,2,3,4,5,6,7]:
            m = fuf.modeGrenander(x, k, p)
            print(p, m)
            self.assertAlmostEqual(m, 0.9, delta=0.1, msg="Problem with modeGrenander")



class MultiVoigtSanity(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testsanity_oneVsMulti(self):
        """
        Checking MultiVoigt1d vs. Voigt1d
        """
        from PyAstronomy import funcFit as fuf
        import numpy as np

        v1 = fuf.Voigt1d()
        vm = fuf.MultiVoigt1d(3)

        x = np.linspace(-10.0, 10.0, 200)

        ps = {
            "lin": 0.0067,
            "off": -12.753,
            "A": 9.0,
            "ad": 2.17,
            "al": 1.78,
            "mu": -0.771,
        }

        v1.assignValues(ps)

        vm["lin"] = v1["lin"]
        vm["off"] = v1["off"]

        for i in range(1, 4):
            num = str(i)
            for p in ["ad", "al", "mu"]:
                vm[p + num] = v1[p]
            vm["A" + num] = v1["A"] / 3.0

        self.assertAlmostEqual(
            np.max(np.abs(vm.evaluate(x) - v1.evaluate(x))),
            0.0,
            delta=1e-12,
            msg="MultiVoigt and Voigt1d deviate with 'amplitude separation'.",
        )

        v1 = fuf.Voigt1d()
        vm = fuf.MultiVoigt1d(3)

        x = np.linspace(-10.0, 10.0, 200)

        ps = {
            "lin": 0.0067,
            "off": -12.753,
            "A": 9.0,
            "ad": 2.17,
            "al": 1.78,
            "mu": -0.771,
        }

        v1.assignValues(ps)

        vm["lin"] = v1["lin"]
        vm["off"] = v1["off"]

        for i in range(1, 4):
            num = str(i)
            for p in ["ad", "al", "mu"]:
                vm[p + num] = v1[p]
            if i == 1:
                vm["A" + num] = v1["A"]

        self.assertAlmostEqual(
            np.max(np.abs(vm.evaluate(x) - v1.evaluate(x))),
            0.0,
            delta=1e-12,
            msg="MultiVoigt and Voigt1d deviate with one nonvanishing profile in multiVoigt.",
        )
