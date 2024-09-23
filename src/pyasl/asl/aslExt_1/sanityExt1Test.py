from __future__ import print_function, division
import unittest
from PyAstronomy import pyasl


class TestSanityOfPyaslExt1(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testsanity_pizzolatoExample(self):
        """
        Example of Pizzolato 2003 relations
        """
        from PyAstronomy import pyasl
        import matplotlib.pylab as plt
        import numpy as np

        p = pyasl.Pizzolato2003()

        # Define array of rotation periods [days]
        prot = np.arange(0.2, 30, 0.1)

        lx = np.zeros(prot.size)
        lxlbol = np.zeros(prot.size)

        # B-V color of star
        bv = 0.7

        # Obtain ...
        for i in range(prot.size):
            # ... log10 of X-ray luminosity
            lx[i] = p.log10lxbv(bv, prot[i])[0]
            # ... and log10(Lx/Lbol)
            lxlbol[i] = p.log10lxlbolbv(bv, prot[i])[0]

        # Plot result
        plt.subplot(2, 1, 1)
        plt.plot(prot, lx, 'bp-')
        plt.subplot(2, 1, 2)
        plt.plot(prot, lxlbol, 'bp-')
#     plt.show()

    def testsanity_pizzolato(self):
        """
        Pizzolato 2003 sanity
        """
        from PyAstronomy import pyasl
        import matplotlib.pylab as plt
        import numpy as np

        p = pyasl.Pizzolato2003()

        x = p.log10lxbv(0.87, 1.0)
        self.assertAlmostEqual(
            x[0], 29.9, msg="Lx saturation level for b-v=0.87 does not match.", delta=1e-7)
        self.assertAlmostEqual(
            x[1], 0.3, msg="Lx error for saturation level for b-v=0.87 does not match.", delta=1e-7)

        x = p.log10lxlbolmass(1.15, 3.0)
        self.assertAlmostEqual(
            x[0], -4.1694, msg="Pizzolato relation (m=1.15, pr=3.0) failed (value=" + str(x[0]) + ").", delta=1e-4)

    def testsanity_expCorrRN_Example1(self):
        """
        Sanity of example 1 for expCorrRN
        """
        import matplotlib.pylab as plt
        from PyAstronomy import pyasl

        # Generate 200 exponentially correlated Gaussian
        # random numbers with a decay time of 5
        c1 = pyasl.expCorrRN(200, 5)

        # Generate 200 exponentially correlated Gaussian
        # random numbers with decay time 10, mean 4, and
        # standard deviation of 2.3.
        #
        # The results are: The correlated random numbers,
        # the uncorrelated numbers used as input, and the
        # correlated coefficient (exp(-1/tau)).
        c2, g, f = pyasl.expCorrRN(200, 10, mean=4.0, std=2.3, fullOut=True)

        plt.subplot(2, 1, 1)
        plt.plot(range(200), c1, 'bp-')
        plt.subplot(2, 1, 2)
        plt.plot(range(200), c2, 'bp-')
        plt.plot(range(200), g, 'g.')
#     plt.show()

    def testsanity_expCorrRN_Example2(self):
        """
        Sanity of example 2 for expCorrRN
        """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy import pyasl

        # Generate n exponentially correlated Gaussian
        # random numbers with a decay time, tau
        n = 500
        tau = 5.
        c1 = pyasl.expCorrRN(n, tau)

        # Obtain autocorrelation function
        ac = np.correlate(c1, c1, mode="full")[n - 1:]

        # Plot correlated random numbers and autocorrelation
        # function along with exponential model.
        x = np.arange(float(n))
        plt.subplot(2, 1, 1)
        plt.plot(x, c1, 'bp-')
        plt.subplot(2, 1, 2)
        plt.plot(x, ac, 'b.')
        plt.plot(x, np.exp(-x / tau) * ac.max(), 'r--')
    #   plt.show()

    def testsanity_ramirez2005Example(self):
        """
        Check sanity of Ramirez/Melendez 2005 example
        """
        from PyAstronomy import pyasl

        # Create class instance
        r = pyasl.Ramirez2005()

        # Which color bands are available
        print("Available color bands: ", r.availableBands())

        # Convert B-V to effective temperature and back
        bv = 0.75
        feh = 0.0
        teff = r.colorToTeff("B-V", bv, feh)
        bv1 = r.teffToColor("B-V", teff, feh)
        # Watch out for differences between input bv and the output bv1
        print("B-V = ", bv, ", Teff = ", teff,
              ", bv1 = ", bv1, ", bv-bv1 = ", bv - bv1)

    def testsanity_ramirez2005Example2(self):
        """
        Checking Ramirez, Ballesteros example
        """
        from PyAstronomy import pyasl

        b = pyasl.BallesterosBV_T()
        r = pyasl.Ramirez2005()

        # Convert B-V to effective temperature and back
        for bv in [0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45]:
            tr = r.colorToTeff("B-V", bv, 0.0)
            tb = b.bv2T(bv)

            print(("B-V [mag] = {3:4.2f} : Teff (R05) = {0:4.0f} K, " +
                   "Teff (B12) = {1:4.0f} K, dTeff = {2: 4.0f} K").format(tr, tb, tr - tb, bv))

    def testsanity_ramirez2005(self):
        """
        Check sanity Ramirez/Melendez 2005 conversions
        """
        from PyAstronomy import pyasl
        import numpy as np
        from itertools import product

        feh = 0.0
        r = pyasl.Ramirez2005()
        # Compare with numbers given in Allen (astrophysical quantities, 4th edition)
        for teff, bv in zip([7000., 5950., 5310., 4410.], [0.35, 0.58, 0.74, 1.15]):
            self.assertAlmostEqual(bv, r.teffToColor("B-V", teff, feh), delta=0.03,
                                   msg="B-V conversion (Teff " + str(teff) + ") out of sync with Allen")

        for teff, bv in zip([6500., 6000., 5000.], [0.286, 0.36, 0.535]):
            self.assertAlmostEqual(bv, r.teffToColor("b-y", teff, feh), delta=0.01,
                                   msg="B-V conversion (Teff " + str(teff) + ") out of sync with Allen")

        for teff, bv in zip([5000., 4500., 4000.], [0.433, 0.51, 0.735]):
            self.assertAlmostEqual(bv, r.teffToColor("R_C-I_C", teff, feh, 'g'), delta=0.04,
                                   msg="R_C-I_C conversion (Teff " + str(teff) + ") out of sync with Allen (giants)")

        for teff, bv in zip([5050., 4660., 4050.], [0.86, 1., 1.5]):
            self.assertAlmostEqual(bv, r.teffToColor("B-V", teff, feh, 'g'), delta=0.1,
                                   msg="B-V conversion (Teff " + str(teff) + ") out of sync with Allen (giants)")

        # Check consistency
        bv = np.arange(0.4, 1.01, 0.05)
        feh = np.arange(-2., 0.3, 0.2)
        for x in product(bv, feh):
            teff = r.colorToTeff("B-V", x[0], x[1], 'ms')
            col = r.teffToColor("B-V", teff, x[1], 'ms')
            self.assertAlmostEqual(x[0], col, delta=1e-2,
                                   msg="Two-way color-temperature conversion is inconsistent. bv = %.4e, feh = %.4e" % x)

    def testsanity_BallesterosExample(self):
        """
        Checking example for Ballesteros
        """
        from PyAstronomy import pyasl

        b = pyasl.BallesterosBV_T()

        bv = 0.65

        # Convert B-V into effective temperature
        teff = b.bv2T(0.65)
        print("B-V = {0:4.2f} mag -> Teff = {1:4.0f} K".format(bv, teff))

        # Convert effective temperature into B-V color
        teff = 4568.0
        bv = b.t2bv(teff)
        print("Teff = {0:4.0f} K -> B-V = {1:4.2f} mag".format(teff, bv))

    def testsanity_gyroAgeBarnes(self):
        """
        Check Gyro age of Barnes 2007
        """
        # Eq 17, 18 (B-V, period, relative error)
        dat = [[0.5, 7, 20], [0.65, 12, 15], [1.0, 20, 13], [1.5, 30, 13]]
        for d in dat:
            a, e = pyasl.gyroAgeBarnes(d[1], d[0])
            self.assertAlmostEqual(e / a * 100., d[2], delta=1.0,
                                   msg="Relative errors of Barnes 2007 cannot be reproduced.")

    def testsanity_sanity_gyroAgeBarnesExample(self):
        """
        Example gyro age (Barnes 2007)
        """
        from PyAstronomy import pyasl

        # Parameters of the Sun (Barnes 2007, p 1174)
        bv = 0.642
        p = 26.09

        # Obtain solar age ...
        age = pyasl.gyroAgeBarnes(p, bv)
        # ... and print it
        print("Solar age: {0:4.2f} +/- {1:4.2f} Ga".format(*age))

    def testsanity_chromoAgeRHK(self):
        """
        Checking sanity of chromospheric age (Donahue)
        """
        import numpy as np
        # R'HK vs log10(age) (estimated from Fig. 1 in Donahue 1998)
        dat = [[-4.25, 7.], [-4.6, 9.1], [-5.0, 9.7]]
        for d in dat:
            a = np.log10(pyasl.chromoAgeRHK(d[0]) * 1e9)
            print(a, d[1])
            self.assertAlmostEqual(a, d[1], delta=0.3,
                                   msg="Cannot reproduce age estimate from Donahue fig 1.")

    def testsanity_chromoAgeRHKExample(self):
        """
        Example chromospheric age (Donahue)
        """
        from PyAstronomy import pyasl

        # Approximate chromospheric age of the Sun
        print("Solar age: {0:4.2f} Ga".format(pyasl.chromoAgeRHK(-4.95)))

    def testsanity_abundancepatternsExample(self):
        """
        Checking example of AbundancePatterns
        """
        from PyAstronomy import pyasl

        ap = pyasl.AbundancePatterns()

        print("Names of the available abundance patterns:")
        print(ap.availablePatterns())

        print()
        print("Get the Asplund et al. pattern (aspl) as a dictionary using")
        print("atomic number as a key:")
        print(ap.pattern("aspl", form="dict", key="number"))

        print()
        print("Get (relative) number abundance of oxygen using elemental symbol:")
        print(ap.abundance("O", pat="wilm"))
        print("or atomic number")
        print(ap.abundance(8, pat="wilm"))
        
        print()
        for el, ab in ap.patternByMass("aspl", key="symbol").items():
            print(f"Element: {el} with mass fraction of {ab}")

    def testsanity_abundancePatterns(self):
        """
        Checking abundance pattern
        """
        from PyAstronomy import pyasl

        ap = pyasl.AbundancePatterns()

        self.assertAlmostEqual(2.82e-07, ap.abundance("P", pat="feld"), delta=1e-12,
                               msg="Abundance of P (feld) incorrect.")
        p = ap.pattern("wilm", form="dict", key="symbol")
        self.assertAlmostEqual(p["Ti"], 6.46e-08, delta=1e-12,
                               msg="Ti abundance broken (wilm)")

    def testsanity_SpecTypeDeJagerExample1(self):
        """
        Sanity of example 1 (SpecTypeDeJager)
        """
        from PyAstronomy import pyasl

        # Instantiate class object
        sdj = pyasl.SpecTypeDeJager()

        llum, lteff = sdj.lumAndTeff("K0", "V")

        print("Luminosity = {0:4.2f} Lsun".format(10.0**llum))
        print("Effective temperature = {0:6.1f} K".format(10.0**lteff))

    def testsanity_SpecTypeDeJagerExample2(self):
        """
        Sanity of example 2 (SpecTypeDeJager)
        """
#     from __future__ import print_function
        from PyAstronomy import pyasl
        import matplotlib.pylab as plt

        # Instantiate class object
        sdj = pyasl.SpecTypeDeJager()

        # Set luminosity class
        lk = "V"

        # Save spectral types, log(teff), and log(luminosity)
        spts = []
        lteffs = []
        llums = []

        # Save information to annotate abscissa
        xt = []
        xtl = []

        for t in "OBAFGKM":
            for n in range(10):
                if (t == "O") and (n == 0):
                    # Skip the invalid "O0" type
                    continue

                # Save the spectral type
                spts.append(t + str(n))

                # Get log10 of luminosity and effective temperature
                ll, lt = sdj.lumAndTeff(spts[-1], lk)
                # and save to lists
                llums.append(ll)
                lteffs.append(lt)

                # Save location (i.e., number in the list) and
                # spectral for annotating the abscissa
                if (n == 0) or (n == 5):
                    xt.append(len(spts) - 1)
                    xtl.append(spts[-1])


#     ax1 = plt.subplot(2,1,1)
#     # Plot log10(effective temperature)
#     plt.plot(lteffs)
#     plt.ylabel("$\log_{10}$(T$_{eff}$)")
#     plt.setp(ax1, xticks=xt,xticklabels=xtl)
#     ax2 = plt.subplot(2,1,2)
#     # Plot log10(luminosity)
#     plt.plot(llums)
#     plt.ylabel("$\log_{10}$($L/L_{\odot}$)")
#     plt.setp(ax2, xticks=xt,xticklabels=xtl)
#     plt.xlabel("Spectral type")
#     plt.show()

    def testsanity_SpecTypeDeJagerSanity(self):
        """
        Sanity of SpecTypeDeJager
        """
        from PyAstronomy import pyasl

        # Instantiate class object
        sdj = pyasl.SpecTypeDeJager()

        llum, lteff = sdj.lumAndTeff("K0", "V")
        self.assertAlmostEqual(10.0**(lteff - 3.712), 1.,
                               delta=0.05, msg="Mismatch K0 V (teff)")
        self.assertAlmostEqual(10.0**(llum + 0.258), 1.,
                               delta=0.2, msg="Mismatch K0 V (lum)")

        llum, lteff = sdj.lumAndTeff("B7", "Ia")
        self.assertAlmostEqual(10.0**(lteff - 4.072), 1.,
                               delta=0.05, msg="Mismatch B7 Ia (teff)")
        self.assertAlmostEqual(10.0**(llum - 5.123), 1.,
                               delta=0.2, msg="Mismatch B7 Ia (lum)")

        llum, lteff = sdj.lumAndTeff("G4", "IV")
        self.assertAlmostEqual(10.0**(lteff - 3.723), 1.,
                               delta=0.05, msg="Mismatch G4 IV (teff)")
        self.assertAlmostEqual(10.0**(llum - 0.802), 1.,
                               delta=0.2, msg="Mismatch G4 IV (lum)")

        llum1, lteff1 = sdj.lumAndTeff("G4", "V")
        llum2, lteff2 = sdj.lumAndTeff("G4.5", "V")
        llum3, lteff3 = sdj.lumAndTeff("G5", "V")
        self.assertTrue(llum1 > llum2 > llum3,
                        msg="Float subtype, wrong order of luminosities")
        self.assertTrue(lteff1 > lteff2 > lteff3,
                        msg="Float subtype, wrong order of luminosities")

    def testsanity_MMSCETSTableExample(self):
        """
        Sanity of example for MMSCETS table
        """
        from PyAstronomy import pyasl
        import matplotlib.pylab as plt

        # Instantiate class
        m = pyasl.MMSCETSTable()

        # Print the entire data file
        for l in m.getContent():
            print(l, end='')

        print()
        print("Available columns: ", ", ".join(m.availableColumns()))

        # Get the tabulated data as an ASCII table
        td = m.getTable()

        # Plot absolute visual brightness vs. effective temperature
        plt.plot(td["Teff"], td["Mv"], 'b.-')
        plt.xlabel("Teff [K]")
        plt.ylabel("Mv [mag]")
#     plt.show()

    def testsanity_Roche(self):
        """
        Check Roche
        """
        from PyAstronomy import pyasl
        import numpy as np
        
        # From Monacki 1984 (Table 1)
        qdat = ((0.01, 3.16664, 3.15345, 0.85853, 1.15622),
                (0.23, 3.78419, 3.54726, 0.64581, 1.45853),
                (0.84, 3.99643, 3.48615, 0.51794, 1.66962),
                (0.54, 3.95665, 3.54129, 0.56301, 1.59529)
            )
        
        for qd in qdat:
            q = qd[0]
            x1, p1 = pyasl.get_lagrange_1(q)
            x2, p2 = pyasl.get_lagrange_2(q)
            self.assertAlmostEqual(qd[1], p1, delta=1e-4, msg="Roche: problem with Monacki " + str(qd))
            self.assertAlmostEqual(qd[2], p2, delta=1e-4, msg="Roche: problem with Monacki " + str(qd))
            self.assertAlmostEqual(qd[3], x1, delta=1e-4, msg="Roche: problem with Monacki " + str(qd))
            self.assertAlmostEqual(qd[4], x2, delta=1e-4, msg="Roche: problem with Monacki " + str(qd))

    def testsanity_roche_example(self):
        """ Cehcking example of Roche """
        from PyAstronomy import pyasl
        import numpy as np
        import matplotlib.pylab as plt
        
        x, y = np.linspace(-1.5,2,300), np.linspace(-1.6,1.6,300)
        xx, yy = np.meshgrid(x, y)
        # Coordinates in orbital plain
        z = 0
        
        # Mass ratio
        q = 0.2
        
        # Get dimensional values of Roche potential
        p = pyasl.rochepot_dl(xx, yy, z, q)
        
        # Positions (and potentials) of Lagrange points
        l1, l1pot = pyasl.get_lagrange_1(q)
        l2, l2pot = pyasl.get_lagrange_2(q)
        l3, l3pot = pyasl.get_lagrange_3(q)
        l4, l5 = pyasl.get_lagrange_4(), pyasl.get_lagrange_5()
        l4pot = pyasl.rochepot_dl(l4[0], l4[1], l4[2], q)
        l5pot = pyasl.rochepot_dl(l5[0], l5[1], l5[2], q)
        
        print("Effective (dimensionless) radii of first and second mass")
        print("According to the approximation of Eggleton 1983:")
        r1eff = pyasl.roche_lobe_radius_eggleton(q, 1)
        r2eff = pyasl.roche_lobe_radius_eggleton(q, 2)
        print("    Reff1: %5.3f" % r1eff)
        print("    Reff2: %5.3f" % r2eff)
        print()
        print("Roche volume and effective radius from Monte Carlo integration:")
        mcvol1 = pyasl.roche_vol_MC(q,1)
        mcvol2 = pyasl.roche_vol_MC(q,2)
        print("    MC Roche lobe volume 1: %6.4f +/- %6.4f" % (mcvol1[0:2]))
        print("    MC Roche lobe volume 2: %6.4f +/- %6.4f" % (mcvol2[0:2]))      
        print("    MC effective radius 1: %6.4f +/- %6.4f" % (mcvol1[2:]))
        print("    MC effective radius 2: %6.4f +/- %6.4f" % (mcvol2[2:]))   
        
#         plt.contour(p, [l5pot*1.02, l3pot, l2pot, l1pot], colors=['g', 'c', 'b', 'r'], extent=[-1.5,2,-1.6,1.6])

    def testsanity_roche(self):
        """ Sanity checks of Roche and root finding """
        from PyAstronomy import pyasl
        from PyAstronomy import constants as PC
        
        def yp(x):
            return (x-0.1)*(x-0.5)*(x-0.7)*(x-0.5000002)
         
        rro = [0.1, 0.5, 0.7]
        roots = sorted(pyasl.bisect_root_find(yp, 0, 1, 4, 1e-4, bqargs={}))
        self.assertEqual(len(rro), len(roots), msg=f"Problem with the number of roots in bisect root finding. Roots are {roots}")
        for i, r in enumerate(roots):
            self.assertAlmostEqual(r, rro[i], 5, msg=f"Problem with bisect root finding. Roots are {roots}")
            
        pc = PC.PyAConstants()
        pc.setSystem("SI")
        q = pc.MEarth/pc.MSun
        rads = pyasl.get_epradius_ss_polar_side(q, pot=None)
        rads = [r*pc.AU/1e3 for r in rads]
        print(rads)
        self.assertAlmostEqual(rads[0], 1.5e6, delta=1e4, msg=f"Problem with L1 height of Earth. Radii are {rads}")

    def testsanity_roche_example2(self):
        """
        Roche lobe example 2 (radii of Earth and hot Jupiter)
        """
        from PyAstronomy import pyasl
        from PyAstronomy import constants as PC
        
        # Shape of the Earth and its Roche lobe (no rotation of Earth)
        pc = PC.PyAConstants()
        pc.setSystem("SI")
        
        # Mass ratio (m2/m1)
        q = pc.MEarth/pc.MSun
        print(f"Earth/Sun mass ratio = {q}")
        
        # Radii of Roche lobe along Earth--Sun connecting line, polar (out-of-plane),
        # and side (in plane)
        rads = pyasl.get_epradius_ss_polar_side(q, pot=None)
        # Convert into km
        rads = [r*pc.AU/1e3 for r in rads]
        print(f"Roche lobe radii (substellar, polar, side) [1e6 km] = " + \
            ", ".join(["%g"%(r/1e6) for r in rads]))
        
        reff_earth = pc.REarth/pc.RJ
        # Small epsilon because rocky Earth is a small body on a large
        # orbit (by Roche standards)
        erads = pyasl.get_radius_ss_polar_side(q, 1, reff_earth, eps=1e-10)
        
        # Convert into km
        erads = [r*pc.RJ/1e3 for r in erads]
        print(f"Earth radii (substellar, polar, side) [km] = " + \
            ", ".join(["%g"%(r) for r in erads[0:3]]))
        
        print()
        print("Roche lobe and planetary shape of typical hot Jupiter")
        # Typical hot Jupiter
        q = 1e-3
        # sma in au
        sma = 0.02
        # Effective (circular disk) transit radius [RJ]
        reff = 1.4
        
        rads = pyasl.get_epradius_ss_polar_side(q, pot=None)
        # Convert into km
        rads = [r*pc.AU/1e3 for r in rads]
        print(f"Roche lobe radii (substellar, polar, side) [1e6 km] = " + \
            ", ".join(["%g"%(r/1e6) for r in rads]))
        
        jrads = pyasl.get_radius_ss_polar_side(q, sma, reff, eps=1e-10)
        print(f"Hot Jupiter radii (substellar, polar, side) [RJ] = " + \
            ", ".join(["%g"%(r) for r in jrads[0:3]]))

    def testsanity_roche_limit_fluid(self):
        """
        Test sanity of roche_limit_fluid
        """
        from PyAstronomy import constants as PC
        from PyAstronomy import pyasl
        pc = PC.PyAConstants()
        pc.setSystem("SI")
        a_r_jup = pyasl.roche_limit_fluid(pc.MJ/pc.MSun, pc.RJ)
        self.assertAlmostEqual(a_r_jup/1e3, 1_713_025, delta=60_000, msg=f"Roche lobe limit for Jupiter does not fit ({a_r_jup} km).")

    def testsanity_sysrem_example(self):
        """
        Checking SysRem example
        """
        from PyAstronomy import pyasl
        import numpy as np
        import matplotlib.pylab as plt
        
        # No. of data sets (e.g., light curves)
        nds = 50
        # No. of data points per data set
        nobs = 200
        
        obs, sigs = [], []
        
        # Generate mock data
        x = np.linspace(0,1,nobs)-0.5
        for i in range(nds):
            # Add noise
            y = np.random.normal(0,0.01+0.0001*i,nobs)
            # Add 3rd degree polynomial
            p = (i, 0.2*i, 0.3*i)
            y += np.polyval(p,x)
            # Add moving Gaussian signal
            y -= 0.02 * np.exp(-(x-(i/nds-0.5))**2/(2*0.02**2))
            # Add growing but stationary Gaussian signal
            y -= (i**2*0.05) * np.exp(-x**2/(2*0.06**2))
            obs.append(y)
            sigs.append(np.ones_like(y)*0.01+0.0001*i)
            
            #plt.plot(x, y, '.-')
        #plt.title("Mock data")
        #plt.show()
        
        sr = pyasl.SysRem(obs, sigs)
        # First iteration
        r, a, c = sr.iterate()
        #plt.subplot(2,1,1)
        #plt.title("Residuals after first (top) and second iteration (bottom)")
        #plt.imshow(r, origin='lower', aspect="auto")
        # Second iteration
        r, a, c = sr.iterate()
        #plt.subplot(2,1,2)
        #plt.imshow(r, origin='lower', aspect="auto")
        #plt.show()

    def testsanity_sysrem_lindata(self):
        """
        Checking SysRem sanity for linear data
        """
        from PyAstronomy import pyasl
        import numpy as np
        import matplotlib.pylab as plt
        
        np.random.seed(123)
        
        # No. of data sets (e.g., light curves)
        nds = 50
        # No. of data points per data set
        nobs = 200
        
        obs, sigs = [], []
        
        a = np.linspace(0,1,200)
        c = np.random.random(50)
        
        m = np.outer(a,c)
        
        # Generate mock data
        for i in range(nds):
            obs.append(m[::,i])
            sigs.append(np.ones_like(a)*0.01+0.0001*i)
            
        
        sr = pyasl.SysRem(obs, sigs)
        # First iteration
        print(np.std(m))
        r, a, c = sr.iterate()
        print(np.std(r))
        self.assertAlmostEqual(np.std(r), 0, msg="SysRem problem with linear data", delta=1e-10)

    def testsanity_parallacticAngle_example(self):
        """ Checking example of parallactic angle """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy import pyasl
        
        # Geo-latitude of Keck [deg]
        gl = 19.83
        
        # Declinations to be considered
        decs = np.arange(-80, 81, 20)
        # Hour angles to be considered
        t = np.linspace(0.0, 8*15, 1000)
        
        for dec in decs:
            q = pyasl.parallacticAngle(t, gl, dec)
            plt.plot(t/15, q, '-', label="$\delta = $" + f"{dec} deg")
        
        plt.ylabel("Parallactic angle [deg]")
        plt.xlabel("Hour angle [h]")
        plt.legend()
        #plt.show()
        
    def testsanity_sysrem_pca(self):
        """ Check sanity of SYSREM with PCA """
        from PyAstronomy import pyasl
        import numpy as np
        from sklearn.decomposition import PCA
        # n observations (e.g., light curves) with m data points each
        n = 4
        m = 7
        
        # Some arbitrary observations (with observations is COLUMNS)
        obs = np.zeros( (m,n) )
        for i in range(0, n):
            for j in range(m):
                obs[j,i] = j*i**3+j*i**2+(j+i+1)
        # Equal error for all data points
        sigs = np.ones_like(obs)

        pca = PCA()
        res = pca.fit(obs.T)
        
        p1 = pca.components_[0,::]
        p2 = pca.components_[1,::]
        
        sr = pyasl.SysRem(obs, sigs, ms_obs=False, ms_feat=True)
        r1, a1, c1 = sr.iterate()
        r2, a2, c2 = sr.iterate()
        a1 = a1/np.linalg.norm(a1)
        a2 = a2/np.linalg.norm(a2)
        
        sign = lambda x: -1 if x < 0 else 1

        # Get sign right
        if sign(a1[0]) != sign(p1[0]):
            a1 *= -1
        if sign(a2[0]) != sign(p2[0]):
            a2 *= -1    
        
        np.testing.assert_almost_equal(a1, p1, decimal=5, err_msg=f'SYSREM: Problem with PCA component 1. Compared {p1} and {a1} (p1 and a1)', verbose=True)
        np.testing.assert_almost_equal(a2, p2, decimal=5, err_msg=f'SYSREM: Problem with PCA component 2. Compared {p2} and {a2} (p2 and a2)', verbose=True)
                
    def testsanity_sysrem_example(self):
        """
        Checking sanity of SYSREM example
        """
        from PyAstronomy import pyasl
        import numpy as np
        import matplotlib.pylab as plt
        
        # No. of data sets (e.g., light curves)
        nds = 50
        # No. of data points per data set
        nobs = 200
        
        obs, sigs = [], []
        
        # Generate mock data
        x = np.linspace(0,1,nobs)-0.5
        for i in range(nds):
            # Add noise
            y = np.random.normal(0,0.01+0.0001*i,nobs)
            # Add 3rd degree polynomial
            p = (i, 0.2*i, 0.3*i)
            y += np.polyval(p,x)
            # Add moving Gaussian signal
            y -= 0.02 * np.exp(-(x-(i/nds-0.5))**2/(2*0.02**2))
            # Add growing but stationary Gaussian signal
            y -= (i**2*0.05) * np.exp(-x**2/(2*0.06**2))
            obs.append(y)
            sigs.append(np.ones_like(y)*0.01+0.0001*i)
        
            plt.plot(x, y, '.-')
        plt.title("Mock data")
        # plt.show()
        
        sr = pyasl.SysRem(obs, sigs)
        # First iteration
        r, a, c = sr.iterate()
        # plt.subplot(2,1,1)
        # plt.title("Residuals after first (top) and second iteration (bottom)")
        # plt.imshow(r, origin='lower', aspect="auto")
        # Second iteration
        r, a, c = sr.iterate()
        # plt.subplot(2,1,2)
        # plt.imshow(r, origin='lower', aspect="auto")
        # plt.show()

    def testsanity_scaleHeight(self):
        """
        Checking sanity of scale height
        """
        from PyAstronomy import pyasl
        T, mu, g = 290, 28.97, 9.8
        she = pyasl.atmosphericScaleHeight(T, mu, g)
        she2 = pyasl.atmosphericScaleHeight_MR(290, 28.97, 1, 1, "E")
        self.assertAlmostEqual(she/she2, 1, delta=1e-3, msg="Scale height of Earth cross check failed")
        
        T, mu, mp, rp = 165, 2.2, 1, 1
        shj = pyasl.atmosphericScaleHeight_MR(T, mu, mp, rp, "J")
        self.assertAlmostEqual(shj, 25.152, delta=1e-3, msg=f"Scale height of Jupiter does not match. result obtained is {shj}")
    
    def testsanity_scaleHeight_example(self):
        """
        Checking example of scale height
        """
        from PyAstronomy import pyasl
    
        T, mu, g = 290, 28.97, 9.8
        she = pyasl.atmosphericScaleHeight(T, mu, g) 
        
        print("Earth")
        print(f"T, mu, g = {T} K, {mu}, {g} m/s**2")
        print(f"Scale height = {she:4.1f} [km]")
        
        
        T, mu, mp, rp = 165, 2.2, 1, 1
        shj = pyasl.atmosphericScaleHeight_MR(T, mu, mp, rp, "J")
        
        print("Jupiter")
        print(f"T, mu, mp, rp = {T} K, {mu}, {mp} [MJ], {rp} [RJ]")
        print(f"Scale height = {shj:4.1f} [km]")

    def testsanity_solarDirectFluxMeinelExample(self):
        """
        Checking example of solarDirectFluxMeinel
        """
        from PyAstronomy import pyasl
        import numpy as np
        # import matplotlib.pylab as plt
        
        z = np.linspace(0,89.9,200)
        
        for h in [0,2,4]:
        
            f = pyasl.solarDirectFluxMeinel(z, height=h)
        
            # plt.plot(z, f, '-', label=f"Height = {h} km")
            
        # plt.legend()
        # plt.xlabel("Zenith angle [deg]")
        # plt.ylabel("Direct solar flux [W/m**2]")
        # plt.show()