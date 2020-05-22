from __future__ import print_function, division
from .airtovac import airtovac, vactoair, vactoair2, airtovac2
import numpy
import unittest
from PyAstronomy.pyaC import SaniBase
from .aitoffLegacy import inverseAitoff, aitoff
from .astroTimeLegacy import daycnv, juldate
from .keplerOrbit import KeplerEllipse, MarkleyKESolver
import datetime
from .dopplerShift import dopplerShift
from .fluxConversion import flux2photons, photons2flux
from .rotBroad import rotBroad, fastRotBroad
from .cardinalPoint import getCardinalPoint
from .posAngle import positionAngle
from .coordinates import hmsToDeg, degToHMS, degToDMS, dmsToDeg, coordsSexaToDeg, coordsDegToSexa
import six
import six.moves as smo
from PyAstronomy.pyaC import pyaErrors as PE


class SanityOfPyasl(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_airtovacExample1(self):
        """
        Checking conversion example for `airtovac2`
        """
        from PyAstronomy import pyasl
        import numpy as np

        # Define wavelength array
        wvl = np.arange(10) + 5000.0
        print("Input wavelengths: ", wvl)

        # Convert wavelength in air to wavelength
        # in vacuum. By default, the conversion
        # specified by Ciddor 1996 are used.
        wvlVac = pyasl.airtovac2(wvl)
        print("Wavelength in vacuum: ", wvlVac)

        # Convert wavelength from vacuum to air
        wvlAir = pyasl.vactoair2(wvlVac)
        print("Wavelength in air: ", wvlAir)

    def sanity_airtovacExample2(self):
        """
        Checking difference example for `airtovac2`
        """
        from PyAstronomy import pyasl
        import numpy
        import matplotlib.pylab as plt

        # Compare the Edlen and Ciddor conversions
        wvl = numpy.arange(3000, 10000, 1.0)
        wvlVacCiddor = pyasl.airtovac2(wvl, mode="ciddor")
        wvlVacEdlen = pyasl.airtovac2(wvl, mode="edlen53")

        plt.subplot(2, 1, 1)
        plt.title("Difference in air wavelength (Ciddor-Edlen)")
        plt.ylabel("dWvl [A]")
        plt.xlabel("Vacuum wvl [A]")
        plt.plot(wvl, wvlVacCiddor - wvlVacEdlen, 'b.-')
        plt.subplot(2, 1, 2)
        plt.title("Difference in air wavelength (Ciddor-Edlen, in RV)")
        plt.ylabel("dRV [m/s]")
        plt.xlabel("Vacuum wvl [A]")
        plt.plot(wvl, (wvlVacCiddor - wvlVacEdlen) /
                 wvlVacCiddor * 299792458., 'b.-')
#    plt.show()

    def sanity_airtovac(self):
        """
        Checking `airtovac` and back...
        """
        wvl_vac = numpy.arange(700.0) * 10.0 + 3001.0
        for mode in ["ciddor", "edlen53", "peckReeder"]:
            wvl_air = vactoair2(wvl_vac)
            wvl_vac2 = airtovac2(wvl_air)
            self.assertTrue(self.mrd(wvl_vac, wvl_vac2) < 1e-10,
                            msg="Problem with vac<->air conversion (mode = " + str(mode) + ")")

    def sanity_specAirVacConvertExample(self):
        """
        Checking example of "specAirVacConvert".
        """
        from PyAstronomy import pyasl
        from PyAstronomy import funcFit as fuf
        import numpy as np
        import matplotlib.pylab as plt

        g = fuf.GaussFit1d()
        g["A"] = -0.1
        g["sig"] = 0.03
        g["mu"] = 5004.4752
        g["off"] = 1.0

        # Create a "spectrum" ...
        wvl = np.linspace(5000., 5010., 1000)
        flux = g.evaluate(wvl)
        # ... and add some noise
        flux += np.random.normal(0., 0.02, len(wvl))

        # Assume that this spectrum is one observed
        # or modeled in vacuum conditions. We want
        # to convert it into a spectrum in air, but
        # retain the old wavelength axis.
        airflux, vind = pyasl.specAirVacConvert(wvl, flux,
                                                direction="vactoair")

        # Plot the result
        plt.plot(wvl, flux, "b.-")
        plt.plot(wvl, airflux, "r.-")
#    plt.show()

    def sanity_refractiveIndexExample(self):
        """
        Checking example of refractive index.
        """
        from PyAstronomy import pyasl
        import numpy as np

        # Define wavelength array
        wvl = np.arange(10) + 5000.0

        # Obtain refractive index according to Ciddor 1996
        ri = pyasl.RefractiveIndexAV(mode="ciddor")
        n = ri.refractiveIndex(wvl)

        print("Wavelength and 1.0 - Refractive index of 'standard air':")
        for w, nc in zip(wvl, n):
            print("{0:5.1f}  {1:10.8e}".format(w, nc - 1.0))

    def sanity_aitoff(self):
        """
        Checking Aitoff projection.
        """
        for i in smo.range(100):
            l = numpy.random.random() * 360.0
            b = numpy.random.random() * 180.0 - 90.0
            x, y = aitoff(l, b)
            l2, b2 = inverseAitoff(x, y)
            self.assertTrue(self.mrd(l2, l) < 1e-6)
            self.assertTrue(self.mrd(b2, b) < 1e-6)

    def sanity_aitoffExample(self):
        """
        Checking example of aitoff projection
        """
        from PyAstronomy import pyasl

        # Define longitude and latitude in degrees
        l = 130.
        b = -35.
        print("Input - Longitude: %4d deg, Latitude: %4d deg" % (l, b))

        print("Aitoff project them and ...")
        x, y = pyasl.aitoff(l, b)
        print(x, y)

        print("... get them back.")
        l2, b2 = pyasl.inverseAitoff(x, y)
        print(l2, b2)

    def sanity_daycnv(self):
        """
        Checking `daycnv`.
        """
        y, m, d, h = tuple(daycnv(2451596.812500))
        self.assertEqual(y, 2000)
        self.assertEqual(m, 2)
        self.assertEqual(d, 22)
        self.assertEqual(h, 7.5)

        jd = 2451596.812500
        # Subtract 30 min
        jd -= 30. / (24. * 60.)
        for i in range(1, 59):
            jd2 = jd + float(i) / (24. * 60.)
            y, m, d, h, minu, sec, ms = tuple(daycnv(jd2, mode='dtlist'))
            minu += round(sec / 60.)
            self.assertAlmostEqual(minu, i, delta=1e-7,
                                   msg="Problem with minutes in daycnv")

        import numpy as np
        jds = [2345678., 2345679.]
        dts = daycnv(jds, mode="dt")
        self.assertEqual((dts[0] - datetime.datetime(1710, 2, 23, 12)).total_seconds(), 0.0,
                         msg="Wrong date returned for JD 2345678 (mode dt)")
        self.assertEqual((dts[1] - dts[0]).total_seconds(), 86400.,
                         msg="Wrong distance between JDs 2345678 and 2345679 (mode dt)")

        ll = daycnv(jds, mode="dtlist")
        self.assertEqual(np.sum(np.abs(np.array(ll[1]) - np.array(ll[0]))), 1,
                         msg="Wrong distance between JDs 2345678 and 2345679 (mode dtlist)")

        fa = daycnv(jds, mode="idl")
        self.assertEqual(np.sum(np.abs(np.array(fa[1]) - np.array(fa[0]))), 1,
                         msg="Wrong distance between JDs 2345678 and 2345679 (mode idl)")

        jds = np.array([2345678., 2345678. + 19876.0 / 86400.])
        dts = daycnv(jds, mode="dt")
        self.assertEqual(np.round((dts[1] - dts[0]).total_seconds()), 19876.,
                         msg="Distance of 19876 s not reproduced by daycvn.")

    def sanity_daycnvExample(self):
        """
        Example of daycnv
        """
        from PyAstronomy import pyasl

        # Convert JD to calendar date
        jd = 2440000.0 + 18614. / (24. * 3600.)
        print("year = %4d, month = %2d, day = %2d, hour = %5.3f"
              % tuple(pyasl.daycnv(jd)))
        print()
        print("year = %4d, month = %2d, day = %2d, hour = %2d, minute = %2d, seconds = %2d, microseconds = %6d"
              % tuple(pyasl.daycnv(jd, mode='dtlist')))
        print()
        dt = pyasl.daycnv(jd, mode='dt')
        print("Datetime object: ", dt)

    def sanity_daycnv_etc(self):
        """
        Checking example of daycnv etc...
        """
        from PyAstronomy import pyasl
        import datetime

        # Convert JD to calendar date
        jd = 2440000.0
        print("year = %4d, month = %2d, day = %2d, hour = %5.3f"
              % tuple(pyasl.daycnv(jd)))
        print()

        # Convert calendar date to JD
        dt = datetime.datetime(1968, 5, 23, 12)
        print("Input date: ", dt)
        print("Corresponding Julian date: ", pyasl.jdcnv(dt))
        print("Corresponding reduced Julian date: ", pyasl.juldate(dt))
        print("Corresponds to weekday: ", pyasl.weekday(dt))
        print()

        # Get current Julian date and convert to heliocentric date
        cjd = pyasl.get_juldate()
        print("Current Julian date: ", cjd)
        print("Current (reduced) heliocentric Julian date (ra=100 deg, dec=37 deg): ",
              pyasl.helio_jd(cjd - 2.4e6, 100.0, 37.0))

    def sanity_juldate(self):
        """
        Checking `juldate`.
        """
        dt = datetime.datetime(2000, 2, 22, 7, 30)
        jd = juldate(dt)
        self.assertEqual(jd, 2451596.812500 - 2400000.0)

    def sanity_keplerOrbitPos(self):
        """
        Checking Kepler orbit.
        """
        time = numpy.linspace(0, 6.78, 200)
        ke = KeplerEllipse(4.5, 3.14159, 0.314)
        pos = ke.xyzPos(time)
        for i in smo.range(len(time)):
            pos1 = ke.xyzPos(time[i])
            self.assertEqual(numpy.sum(numpy.abs((pos[i, ::] - pos1))), 0.0)

    def sanity_keplerAngMom(self):
        """
        Checking Kepler orbit: Angular momentum
        """
        import numpy as np
        ke = KeplerEllipse(1, 2.0 * np.pi, e=0.0, i=90.0,
                           Omega=0., w=0.0, tau=0.0)
        am1 = ke.orbAngMomentum()
        am2 = ke.orbAngMomentum(t=101.7)
        d = np.max(np.abs(am1 - am2))
        self.assertAlmostEqual(d, 0.0, delta=1e-12,
                               msg="Angular momentum is not constant in time")

        d = np.max(np.abs(am1 - np.array([0, -1, 0])))
        self.assertAlmostEqual(d, 0.0, delta=1e-12,
                               msg="Angular momentum wrong for i=90 deg")

        ke = KeplerEllipse(1, 2.0 * np.pi, e=0.0, i=0.0,
                           Omega=0., w=0.0, tau=0.0)
        am1 = ke.orbAngMomentum()
        d = np.max(np.abs(am1 - np.array([0, 0, 1])))
        self.assertAlmostEqual(d, 0.0, delta=1e-12,
                               msg="Angular momentum wrong for i=0 deg")

        for w in [10.0, 87.0, 277.8]:
            ke = KeplerEllipse(1, 2.0 * np.pi, e=0.0, i=0.0,
                               Omega=0., w=w, tau=0.0)
            ame = ke.orbAngMomentum()
            d = np.max(np.abs(ame - am1))
            self.assertAlmostEqual(
                d, 0.0, delta=1e-12, msg="Angular momentum changes with argument of periapsis")

    def sanity_keplerOrbitNodes(self):
        """
        Checking node position of Kepler ellipse.
        """
        import numpy as np
        ke = KeplerEllipse(1, 1, e=0.0, i=90.0, Omega=0., w=0.0, tau=0.0)

        n1 = ke.xyzNodes_LOSZ(los='+z')
        self.assertAlmostEqual(max(np.abs(n1[0] - np.array([1.0, 0.0, 0.0]))), 0.0, delta=1e-15,
                               msg="Location of ascending node does not match.")
        self.assertAlmostEqual(max(np.abs(n1[1] - np.array([-1.0, 0.0, 0.0]))), 0.0, delta=1e-15,
                               msg="Location of descending node does not match.")

        n2 = ke.xyzNodes_LOSZ(los='-z')
        for i in range(2):
            md = np.max(np.abs(n1[i] - n2[(i + 1) % 2]))
            self.assertAlmostEqual(
                md, 0.0, delta=1e-15, msg="Nodes are not exchanged when LOS is reversed.")

        # None of this should have an impact on the node
        for w in smo.range(0, 360, 80):
            for tau in smo.range(0, 100, 25):
                for i in smo.range(5, 85, 20):
                    ke2 = KeplerEllipse(
                        1, 1, e=0.0, i=i, Omega=0., w=w, tau=tau)
                    n3 = ke2.xyzNodes_LOSZ(los='+z')
                    for j in smo.range(2):
                        md = np.max(np.abs(n1[j] - n3[j]))
                        self.assertAlmostEqual(md, 0.0, delta=1e-12,
                                               msg="Unexpected impact on position of node.")

        # Check node rotation
        ke3 = KeplerEllipse(1, 1, e=0.0, i=90., Omega=90., w=0.0, tau=0.0)
        n4 = ke3.xyzNodes_LOSZ(los='+z')
        self.assertAlmostEqual(np.max(np.abs(n4[0] - np.array([0., 1., 0.]))), 0.0, delta=1e-12,
                               msg="Wrong ascending node for Omega = 90")
        ke3 = KeplerEllipse(1, 1, e=0.0, i=90., Omega=180., w=0.0, tau=0.0)
        n4 = ke3.xyzNodes_LOSZ(los='+z')
        self.assertAlmostEqual(np.max(np.abs(n4[0] - np.array([-1., 0., 0.]))), 0.0, delta=1e-12,
                               msg="Wrong ascending node for Omega = 180")
        ke3 = KeplerEllipse(1, 1, e=0.0, i=90., Omega=270., w=0.0, tau=0.0)
        n4 = ke3.xyzNodes_LOSZ(los='+z')
        self.assertAlmostEqual(np.max(np.abs(n4[0] - np.array([0., -1., 0.]))), 0.0, delta=1e-12,
                               msg="Wrong ascending node for Omega = 270")

    def sanity_MarkleyKESolver_precision(self):
        """
        Checking precision of Markley solver for Kepler's equation.
        """
        ks = MarkleyKESolver()
        p = ks.precisionTest()
        self.assertLess(p, 1e-14)

    def sanity_dmstodeg(self):
        """
        Check sanity of dmstodeg
        """
        try:
            dmsToDeg(0, 12, 34)
        except PE.PyAValError:
            # This should occur
            pass
        except Exception as e:
            # Another exception -- not good
            raise(e)
        else:
            # This is a mistake: no exception was raised...
            raise(Exception(
                "dmstodeg should have raised an exception but did not: dmsToDeg(0, 12, 34)"))

        try:
            dmsToDeg(-10, 12, 34, esign=+1)
        except PE.PyAValError:
            # This should occur
            pass
        except Exception as e:
            # Another exception -- not good
            raise(e)
        else:
            # This is a mistake: no exception was raised...
            raise(Exception(
                "dmstodeg should have raised an exception but did not: dmsToDeg(-10, 12, 34, esign=+1)"))

        self.assertAlmostEqual(dmsToDeg(0, 1, 0, -1), -1. / 60., delta=1e-12,
                               msg="dmsToDeg: (0,1,0,-1) not correctly converted")

    def sanity_degToSexaCoordConversion(self):
        """
        Checking degrees-sexagesimal coordinate conversion.
        """
        # From SIMBAD
        hd1s = "00 05 08.83239 +67 50 24.0135"
        hd1d = (1.28680161, +67.84000375)
        hd123456s = "14 08 24.78626 -26 20 56.6438"
        hd123456d = (212.10327608, -26.34906773)

        sexa = coordsDegToSexa(92.5, -21.8)
        self.assertEqual(sexa, "06 10 00.000  -21 48 00.000", msg="Error in conversion of '92.5, -21.8' into sexagesimal coordinates. Obtained: " + str(sexa))

        hdhmss = [(0.0, 5., 08.83239), (14., 8., 24.78626)]
        hddmss = [(67., 50., 24.0135), (-26., 20., 56.6438)]

        hdss = [hd1s, hd123456s]
        hdds = [hd1d, hd123456d]

        for i, hds in enumerate(hdss):
            d = hmsToDeg(*hdhmss[i])
            self.assertAlmostEqual(
                d, hdds[i][0], delta=1e-7, msg="hmsToDeg: Incorrect conversion")
            hms = degToHMS(d)
            for j in range(3):
                self.assertAlmostEqual(
                    hms[j], hdhmss[i][j], delta=1e-5, msg="degToHMS: Incorrect conversion")

            d = dmsToDeg(*hddmss[i])
            self.assertAlmostEqual(
                d, hdds[i][1], delta=1e-7, msg="dmsToDeg: Incorrect conversion")
            dms = degToDMS(d)
            for j in range(3):
                self.assertAlmostEqual(
                    dms[j], hddmss[i][j], delta=1e-5, msg="degToDMS: Incorrect conversion")

            r, d, hms, dms = coordsSexaToDeg(hdss[i], fullOut=True)
            self.assertAlmostEqual(
                r, hdds[i][0], delta=1e-7, msg="coordsSexaToDeg: ra incorrect")
            self.assertAlmostEqual(
                d, hdds[i][1], delta=1e-7, msg="coordsSexaToDeg: dec incorrect")
            for j in range(3):
                self.assertAlmostEqual(
                    hms[j], hdhmss[i][j], delta=1e-5, msg="degToHMS: Incorrect conversion")
            for j in range(3):
                self.assertAlmostEqual(
                    dms[j], hddmss[i][j], delta=1e-5, msg="degToDMS: Incorrect conversion")

            rat, dect = coordsDegToSexa(r, d, asString=False)
            print(rat, dect, r, d)
            for j in range(3):
                self.assertAlmostEqual(
                    rat[j], hdhmss[i][j], delta=1e-5, msg="degToHMS: Incorrect conversion")
            for j in range(3):
                self.assertAlmostEqual(
                    dect[j], hddmss[i][j], delta=1e-5, msg="degToDMS: Incorrect conversion")

    def sanity_degToSexaCoordConversionExample(self):
        """
        Checking example of sexagesimal-decimal conversion
        """
        from PyAstronomy import pyasl

        # Coordinates of HD 1 from SIMBAD
        hd1 = "00 05 08.83239 +67 50 24.0135"

        print("Coordinates of HD 1 (SIMBAD): ", hd1)

        # Obtain decimal representation
        ra, dec = pyasl.coordsSexaToDeg(hd1)
        print("Coordinates of HD 1 [deg]: %010.6f  %+09.6f" % (ra, dec))

        # Convert back into sexagesimal representation
        sexa = pyasl.coordsDegToSexa(ra, dec)
        print("Coordinates of HD 1 [sexa]: ", sexa)

    def sanity_PlanckExample(self):
        """
        Checking sanity of Planck example
        """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy.pyasl import planck

        # Define wavelength in meters
        lam = np.arange(1000.0 * 1e-10, 20000. * 1e-10, 20e-10)

        # Get the Planck spectrum in [W/(m**2 m)] for a temperature of 7000 K
        s7 = planck(7000., lam=lam)
        # Get the Planck spectrum in [W/(m**2 m)] for a temperature of 5000 K
        s5 = planck(5000., lam=lam)

        # Convert into erg/(cm**2 * A * s)
        s5erg = s5 * 1e-7
        s7erg = s7 * 1e-7

        # Integrate the spectrum and compare with Stefan-Boltzmann law
        i5 = np.sum(s5) * (lam[1] - lam[0])
        i7 = np.sum(s7) * (lam[1] - lam[0])

        print("5000 K integral: %.3e W/m**2 (Stefan-Boltzmann predicts %.3e W/m**2)" %
              (i5, (5.67e-8 * 5000.**4)))
        print("7000 K integral: %.3e W/m**2 (Stefan-Boltzmann predicts %.3e W/m**2)" %
              (i7, (5.67e-8 * 7000.**4)))

#     plt.xlabel("Wavelength [$\AA$]")
#     plt.ylabel("Flux [erg/cm$^2$/A/s]")
#     plt.plot(lam*1e10, s5erg, 'r-')
#     plt.plot(lam*1e10, s7erg, 'b-')
#     plt.show()

    def sanity_magToFluxDensity_bessel98(self):
        """
        Checking consistency of magToFluxDensity_bessel98
        """
        import numpy as np
        from PyAstronomy import pyasl
        
        fnu = np.array("1.790 4.063 3.636 3.064 2.416 1.589 1.021 0.640 0.676 0.285 0.238".split(), dtype=np.float)
        flam = np.array("417.5 632 363.1 217.7 112.6 31.47 11.38 3.961 4.479 0.708 0.489".split(), dtype=np.float)
        
        for i, b in enumerate("U, B, V, R, I, J, H, K, Kp, L, L*".split(",")):
            a, lame = pyasl.magToFluxDensity_bessel98(b.strip(), 0.0, "nu")
            self.assertAlmostEqual((a*1e20 - fnu[i])/fnu[i], 0.0, delta=8e-3, msg="fnu does not match for zero magnitude")
            b, lame = pyasl.magToFluxDensity_bessel98(b.strip(), 0.0, "lam")
            print()       
            self.assertAlmostEqual((b*1e11 - flam[i])/flam[i], 0.0, delta=8e-3, msg="flam does not match for zero magnitude")

    def sanity_magToFluxDensity_bessel98_example(self):
        """
        Checking example of magToFluxDensity_bessel98
        """
        from PyAstronomy import pyasl
        import numpy as np
        
        mag_R = 15.5
        
        fd_nu, le = pyasl.magToFluxDensity_bessel98("R", mag_R, "nu")
        fd_lam, _ = pyasl.magToFluxDensity_bessel98("R", mag_R, "lam")
        
        print("R-band magnitude: ", mag_R)
        print("R-band flux density [erg/cm**2/s/Hz]: ", fd_nu)
        print("R-band flux density [erg/cm**2/s/A]: ", fd_lam)
        
        print("Effective wavelength of filter [A]: ", le)
        print("Convert f_nu into f_lam [erg/cm**2/s/A] by multiplication with (c/lam**2): ", \
            fd_nu * (299792458e2/(le/1e8)**2) / 1e8 )


class SanityOfMoonpos(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_exampleMoonposition(self):
        """
        Checking example of moonpos (position)
        """
        import datetime
        from PyAstronomy import pyasl
        import numpy as np

        # Convert calendar date to JD
        # using the datetime package
        jd = datetime.datetime(2013, 4, 16)
        jd = pyasl.jdcnv(jd)
        jd = np.arange(jd, jd + 20, 1)
        # Calculate Moon positions
        res = pyasl.moonpos(jd)

        print("%15s  %8s  %8s  %11s  %8s  %8s" %
              ("JD", "RA", "DEC", "DIST", "GEOLON", "GEOLAT"))
        print("%15s  %8s  %8s  %11s  %8s  %8s" %
              ("[d]", "[deg]", "[deg]", "[km]", "[deg]", "[deg]"))
        for i in range(jd.size):
            print("%15.4f  %8.4f  %8.4f  %11.4f  %8.4f  %8.4f" %
                  (jd[i], res[0][i], res[1][i], res[2][i], res[3][i], res[4][i]))

    def sanity_moonphaseExample(self):
        """
        Checking exmaple of moonphase
        """
        import datetime
        from PyAstronomy import pyasl
        import numpy as np

        # Convert calendar date to JD
        # using the datetime package
        jd = datetime.datetime(2013, 4, 16)
        jd = pyasl.jdcnv(jd)
        jd = np.arange(jd, jd + 20, 1)
        mp = pyasl.moonphase(jd)

        print("%15s  %3s" % ("JD", "Phase"))
        for i in range(jd.size):
            print("%15.4f  %3d%%" % (jd[i], mp[i] * 100.))


class SanityOfKeplerOrbitExamples(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example1(self):
        """
        Checking Markley solver example.
        """
        from PyAstronomy import pyasl

        # Instantiate the solver
        ks = pyasl.MarkleyKESolver()

        # Solves Kepler's Equation for a set
        # of mean anomaly and eccentricity.
        # Uses the algorithm presented by
        # Markley 1995.
        M = 0.75
        e = 0.3
        print("Eccentric anomaly: ", ks.getE(M, e))

    def sanity_example2(self):
        """
        Checking Kepler orbit calculation example.
        """
        import numpy as np
        from PyAstronomy import pyasl
        import matplotlib.pylab as plt

        # Instantiate a Keplerian elliptical orbit with
        # semi-major axis of 1.3 length units,
        # a period of 2 time units, eccentricity of 0.5,
        # longitude of ascending node of 70 degrees, an inclination
        # of 10 deg, and a periapsis argument of 110 deg.
        ke = pyasl.KeplerEllipse(1.3, 2., e=0.5, Omega=70., i=10.0, w=110.0)

        # Get a time axis
        t = np.linspace(0, 1.9, 200)

        # Calculate the orbit position at the given points
        # in a Cartesian coordinate system.
        pos = ke.xyzPos(t)
        print("Shape of output array: ", pos.shape)

        # x, y, and z coordinates for 50th time point
        print("x, y, z for 50th point: ", pos[50, ::])

        # Calculate orbit radius as a function of the
        radius = ke.radius(t)

        # Calculate velocity on orbit
        vel = ke.xyzVel(t)

        # Find the nodes of the orbit (Observer at -z)
        ascn, descn = ke.xyzNodes_LOSZ()

        # Plot x and y coordinates of the orbit
        plt.subplot(2, 1, 1)
        plt.title(
            "Periapsis (red diamond), Asc. node (green circle), desc. node (red circle)")
        plt.xlabel("East ->")
        plt.ylabel("North ->")
        plt.plot([0], [0], 'k+', markersize=9)
        plt.plot(pos[::, 1], pos[::, 0], 'bp')
        # Point of periapsis
        plt.plot([pos[0, 1]], [pos[0, 0]], 'rd')
        # Nodes of the orbit
        plt.plot([ascn[1]], [ascn[0]], 'go', markersize=10)
        plt.plot([descn[1]], [descn[0]], 'ro', markersize=10)
        # Plot RV
        plt.subplot(2, 1, 2)
        plt.xlabel("Time")
        plt.ylabel("Radial velocity [length/time]")
        plt.plot(t, vel[::, 2], 'r.-')
#     plt.show()


class SanityOfBinnin(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example_binningx0dt(self):
        """
        Checking `binningx0dt` example 1.
        """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy.pyasl import binningx0dt

        # Generate some data
        x = np.arange(999)
        y = np.sin(x / 100.)
        y += np.random.normal(0, 0.1, len(x))

        # Bin using fixed number of bins and start at x0 = -10.
        # Use beginning of bin as starting value.
        r1, dt1 = binningx0dt(x, y, nbins=50, x0=-10, useBinCenter=False)
        # Use fixed bin width. Specify another (wrong) error estimate and
        # use bin center.
        r2, dt2 = binningx0dt(x, y, yerr=np.ones(len(x)) * 0.2, dt=dt1,
                              x0=-10, useBinCenter=True, removeNoError=True)

        print("dt1, dt2: ", dt1, dt2)
        print("Input data points in last bin: ", r2[-1, 3])

        # Use the reducedBy flag to indicate the binning. In this case, x0
        # will be set to the lowest x value in the data, and the number of
        # bins will be calculated as: int(round(len(x)/float(reduceBy))).
        # Here, we will, thus, obtain 100 bins.
        r3, dt3 = binningx0dt(x, y,
                              useBinCenter=True, removeNoError=True, reduceBy=10)

        print("dt3: ", dt3)
        print("Number of bins in third version: ", len(r3[::, 0]))

        # Plot the output
        plt.plot(x, y)
        plt.errorbar(r1[::, 0], r1[::, 1], yerr=r1[::, 2], fmt='kp--')
        plt.errorbar(r2[::, 0], r2[::, 1], yerr=r2[::, 2], fmt='rp--')
        plt.errorbar(r3[::, 0], r3[::, 1], yerr=r3[::, 2], fmt='gp--')
#     plt.show()

    def sanity_example_binningx0dt_example2(self):
        """
        Checking `binningx0dt` example 2.
        """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy.pyasl import binningx0dt

        # Generate some data
        x = np.arange(-100, 999)
        # Create some holes in the data
        x = np.delete(x, list(range(340, 490)))
        x = np.delete(x, list(range(670, 685)))
        x = np.delete(x, list(range(771, 779)))
        y = np.sin(x / 100.)
        y += np.random.normal(0, 0.1, len(x))

        # Bin using bin width of 27 and starting at minimum x-value.
        # Use beginning of bin as starting value.
        r1, dt1 = binningx0dt(x, y, dt=27, x0=min(x), useBinCenter=True)

        # As previously, but use the mean x-value in the bins to produce the
        # rebinned time axis.
        r2, dt2 = binningx0dt(x, y, dt=27, x0=min(x), useMeanX=True)

        print("Median shift between the time axes: ",
              np.median(r1[::, 0] - r2[::, 0]))
        print(" -> Time bins are not aligned due to 'forced' positioning of")
        print("    the first axis.")

        # Plot the output
        plt.plot(x, y, 'b.-')
        plt.errorbar(r1[::, 0], r1[::, 1], yerr=r1[::, 2], fmt='kp--')
        plt.errorbar(r2[::, 0], r2[::, 1], yerr=r2[::, 2], fmt='rp--')
#     plt.show()

    def sanity_example_binningx0dt_example3(self):
        """
        Chacking example 3 for binningx0dt (nanHandling)
        """
        from PyAstronomy.pyasl import binningx0dt
        import matplotlib.pyplot as plt
        import numpy as np

        # Set up figures
        fig = plt.figure()
        ax0 = fig.add_subplot(411)
        ax1 = fig.add_subplot(412, sharex=ax0, sharey=ax0)
        ax2 = fig.add_subplot(413, sharex=ax0, sharey=ax0)
        ax3 = fig.add_subplot(414, sharey=ax0)

        # Set up data
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        y = [0.3, 0.5, 0.7, 0.2, 0.5, 0.9, 0.2, 0.7, 0.8, 0.6]
        yerr = [0.1] * len(x)

        r, dt = binningx0dt(x, y, yerr=yerr, x0=0.5, dt=2)
        ax0.plot(x, y, marker='o')
        ax0.plot(r[::, 0], r[::, 1], linestyle='--',
                 drawstyle='steps-mid', marker='s')
        ax0.set_title("Normal (w/o NaNs)")
        ax0.set_xticklabels([])

        y = [0.3, 0.5, np.nan, 0.2, 0.5, 0.9, np.nan, np.nan, 0.8, 0.6]
        x, y = np.array(x), np.array(y)
        r, dt = binningx0dt(x, y, yerr=yerr, x0=0.5, dt=2)
        ax1.plot(x, y, marker='o')
        ax1.plot(r[::, 0], r[::, 1], linestyle='--',
                 drawstyle='steps-mid', marker='s')
        ax1.set_title("With NaNs and nanHandling='None' (default)")
        # ax1.set_xticklabels([])

        r, dt = binningx0dt(x, y, yerr=yerr, x0=0.5,
                            dt=2, nanHandling="ignore")
        ax2.plot(x, y, marker='o')
        ax2.plot(r[::, 0], r[::, 1], linestyle='--',
                 drawstyle='steps-mid', marker='s')
        ax2.set_title("With NaNs and nanHandling='ignore'")

        r, dt = binningx0dt(x, y, x0=0.5, dt=2, nanHandling=0.5)
        ax3.plot(x, y, marker='o')
        ax3.plot(r[::, 0], r[::, 1], linestyle='--',
                 drawstyle='steps-mid', marker='s')
        ax3.set_title("With NaNs and nanHandling=0.5")

        ax0.set_xlim(0, 11.5)
        ax3.set_xlim(0, 11.5)
        ax0.set_ylim(0, 1.1)

#     plt.show()


class SanityOfPhotonConversion(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_fluxConversion(self):
        """
        Checking `flux2photons` and `photons2flux`.
        """
        wvl = numpy.linspace(300., 20000., 100)
        flux = 1. / wvl

        photons = flux2photons(wvl, flux)
        flux2 = photons2flux(wvl, photons)
        self.assertLess(max(numpy.abs(flux - flux2)), 1e-12)

    def sanity_example(self):
        """
        Checking sanity of flux conversion example.
        """
        from PyAstronomy import pyasl

        # Wavelength in Angstrom
        wvl = 4000.
        # Flux in erg/s
        flux = 1.5e-14

        # Convert into photons
        photons = pyasl.flux2photons(wvl, flux)

        # How many photons is this?
        print("%g erg/s at %g A correspond to %g photons/s"
              % (flux, wvl, photons))

        # Converting back
        flux2 = pyasl.photons2flux(wvl, photons)

        print("%g photons/s at %g A correspond to %g erg/s"
              % (photons, wvl, flux2))


class SanityOfFolding(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_foldAtExample(self):
        """
        Checking `foldAt` example.
        """
        from PyAstronomy.pyasl import foldAt
        import matplotlib.pylab as plt
        import numpy as np

        # Generate some data ...
        time = np.random.random(1000) * 100.
        flux = 0.05 * np.sin(time * (2. * np.pi / 21.5) + 15)
        # ... and add some noise
        flux += np.random.normal(0, 0.02, len(flux))

        # Obtain the phases with respect to some
        # reference point (in this case T0=217.4)
        phases = foldAt(time, 21.5, T0=217.4)

        # Sort with respect to phase
        # First, get the order of indices ...
        sortIndi = np.argsort(phases)
        # ... and, second, rearrange the arrays.
        phases = phases[sortIndi]
        flux = flux[sortIndi]

        # Plot the result
        plt.plot(phases, flux, 'bp')
        # plt.show()

    def sanity_foldAtEpoch(self):
        """
        Checking epoch calculation of `foldAt`
        """
        from PyAstronomy.pyasl import foldAt
        import numpy as np
        time = np.array([-1.5, -0.5, 0.5, 1.5, 9.5])
        p, e = foldAt(time, 1.0, T0=0.0, getEpoch=True)
        self.assertFalse(np.any(np.abs(p - 0.5) > 1e-5), msg="Bad phase")
        pe = np.array([-2., -1., 0.0, 1.0, 9.0])
        self.assertFalse(np.any(np.abs(pe - e) > 1e-5), msg="Bad epoch")

    def sanity_foldAtBackAndForth(self):
        """
        Checking relations from foldAt documentation
        """
        from PyAstronomy.pyasl import foldAt
        import numpy as np
        t = np.linspace(-10, 10, 137)
        per = 2.35
        T0 = -1.55
        p, e = foldAt(t, per, T0=T0, getEpoch=True)
        t2 = T0 + (e+p)*per
        self.assertFalse(np.any(np.abs(t2 - t) > 1e-5), msg="foldAt back and forth: Bad conversion!")


class SanityOfDopplerShift(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_backAndForth(self):
        """
        Doppler shift -- shifting back and forth
        """
        wvl = numpy.linspace(6000., 6100., 10000)
        flux = numpy.ones(len(wvl)) + (wvl / wvl.min()) * 0.05
        flux -= numpy.exp(-(wvl - 6050.)**2 / (2. * 0.5**2)) * 0.05

        nflux, wlprime = dopplerShift(wvl, flux, 20., edgeHandling="firstlast")
        nflux2, wlprime = dopplerShift(
            wvl, nflux, -20., edgeHandling="firstlast")

        self.assertLess(max(numpy.abs(flux - nflux2)), 1e-5)

    def sanity_checkExample(self):
        """
        Doppler shift -- Checking whether example works
        """
        from PyAstronomy import pyasl
        import matplotlib.pylab as plt
        import numpy as np

        # Create a "spectrum" with 0.01 A binning ...
        wvl = np.linspace(6000., 6100., 10000)
        # ... a gradient in the continuum ...
        flux = np.ones(len(wvl)) + (wvl / wvl.min()) * 0.05
        # ... and a Gaussian absoption line
        flux -= np.exp(-(wvl - 6050.)**2 / (2. * 0.5**2)) * 0.05

        # Shift that spectrum redward by 20 km/s using
        # "firstlast" as edge handling method.
        nflux1, wlprime1 = pyasl.dopplerShift(
            wvl, flux, 20., edgeHandling="firstlast")

        # Shift the red-shifted spectrum blueward by 20 km/s, i.e.,
        # back on the initial spectrum.
        nflux2, wlprime = pyasl.dopplerShift(wvl, nflux1, -20.,
                                             edgeHandling="fillValue", fillValue=1.0)

        # Check the maximum difference in the central part
        indi = np.arange(len(flux) - 200) + 100
        print("Maximal difference (without outer 100 bins): ",
              max(np.abs(flux[indi] - nflux2[indi])))

        # Plot the outcome
        plt.title(
            "Initial (blue), shifted (red), and back-shifted (green) spectrum")
        plt.plot(wvl, flux, 'b.-')
        plt.plot(wvl, nflux1, 'r.-')
        plt.plot(wvl, nflux2, 'g.-')
#     plt.show()


class SanityOfrotBroad(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example(self):
        """
        Sanity of "rotBroad" example.
        """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy import funcFit as fuf
        from PyAstronomy import pyasl
        import scipy.integrate as sci

        # Create a spectrum with a single Gaussian
        # line using funcFit's GaussFit1d object.
        # Note that this object is not used for
        # fitting here, but only a calculate a
        # Gaussian.
        g = fuf.GaussFit1d()
        g["mu"] = 5005.
        g["A"] = -0.1
        g["sig"] = 0.1
        g["off"] = 1.0

        # Evaluate the spectrum with 0.01 A bin size
        wvl = np.linspace(5003., 5007., 400)
        flux = g.evaluate(wvl)

        # Obtain the broadened spectrum using
        # vsini = 13.3 km/s and no limb-darkening
        rflux = pyasl.rotBroad(wvl, flux, 0.0, 13.3)

        # Obtain the broadened spectrum using
        # vsini = 13.3 km/s and strong limb-darkening
        lflux = pyasl.rotBroad(wvl, flux, 0.9, 13.3)

        # Check that the area of the line did not change
        # in response to the broadening
        print("Initial EW [A]: ", 4. - sci.trapz(flux, wvl))
        print("After broadening without LD: ", 4. - sci.trapz(rflux, wvl))
        print("After broadening with LD: ", 4. - sci.trapz(lflux, wvl))

        # Plot the results
        plt.title("Rotational broadening")
        plt.xlabel("Wavelength [A]")
        plt.ylabel("Normalized flux")
        plt.plot(wvl, flux, 'b-')
        plt.plot(wvl, rflux, 'r-')
        plt.plot(wvl, lflux, 'g-')
#     plt.show()

    def sanity_rotBroadEW(self):
        """
        Checking whether "rotBroad" changes the EW.
        """
        import numpy as np
        from PyAstronomy import funcFit as fuf
        import scipy.integrate as sci

        # Create a spectrum with a single Gaussian
        # line using funcFit's GaussFit1d object.
        # Note that this object is not used for
        # fitting here, but only a calculate a
        # Gaussian.
        g = fuf.GaussFit1d()
        g["mu"] = 5005.
        g["A"] = -0.1
        g["sig"] = 0.1
        g["off"] = 1.0

        # Evaluate the spectrum with 0.01 A bin size
        wvl = np.linspace(5003., 5007., 400)
        flux = g.evaluate(wvl)

        refEW = abs(g["A"])

        for broad in np.linspace(0.1, 100., 10):
            for eps in np.arange(0., 1., 0.2):
                nflux = rotBroad(wvl, flux, epsilon=eps, vsini=broad)
                ew = 4. - sci.trapz(nflux, wvl)
                self.assertAlmostEqual(ew, refEW, 5, "The EW in rotbroad changes for: vsini = " +
                                       str(broad) + ", and eps = " + str(eps))

    def sanity_slowVsFast(self):
        """
        Comparing 'slow' and 'fast' algorithm.
        """
        wvl = numpy.arange(4999., 5011., 0.04)
        flux = numpy.zeros(len(wvl))

        flux[150] = 0.5

        bfast = fastRotBroad(wvl, flux, 0.37, 37.)
        bslow = rotBroad(wvl, flux, 0.37, 37.)

        indi = numpy.where(numpy.logical_and(wvl > 5000., wvl < 5010.))[0]
        self.assertAlmostEqual(
            numpy.max(numpy.abs(bfast[indi] - bslow[indi])), 0.0, delta=1e-5)

        flux[150] = 1.0
        A = -0.05
        s = 0.1
        mu = 5004.1635788
        flux += A / numpy.sqrt(2. * numpy.pi * s**2) * \
            numpy.exp(-(wvl - mu)**2 / (2. * s**2))

        bfast = fastRotBroad(wvl, flux, 0.81, 11.37)
        bslow = rotBroad(wvl, flux, 0.81, 11.37)

        indi = numpy.where(numpy.logical_and(wvl > 5000., wvl < 5010.))[0]
        self.assertAlmostEqual(
            numpy.max(numpy.abs(bfast[indi] - bslow[indi])), 0.0, delta=1e-5)


class SanityOfBaryvel(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example(self):
        """
        Sanity of "baryvel" example
        """
        from PyAstronomy import pyasl

        jd = 2.476468576e6

        heli, bary = pyasl.baryvel(jd, deq=2000.0)

        print("Earth's velocity at JD: ", jd)
        print("Heliocentric velocity [km/s]: ", heli)
        print("Barycentric velocity [km/s] : ", bary)

        # Coordinates of Sirius
        ra = 101.28715535
        dec = -16.71611587

        vh, vb = pyasl.baryCorr(jd, ra, dec, deq=2000.0)
        print("Barycentric velocity of Earth toward Sirius: ", vb)

    def sanity_baryCorrAltair(self):
        """
        Checking barycentric correction for Altair (IDL example)
        """
        from PyAstronomy import pyasl
        vh, vb = pyasl.baryCorr(
            2449398.5000000000, 297.69488, 8.8676389, deq=2000)
        self.assertAlmostEqual(vb, 10.583091, delta=1e-5)

    def sanity_helcorrExample(self):
        """
        Checking sanity of helcorr example.
        """
        from PyAstronomy import pyasl

        # Coordinates of European Southern Observatory
        # (Coordinates of UT1)
        longitude = 289.5967661
        latitude = -24.62586583
        altitude = 2635.43

        # Coordinates of HD 12345 (J2000)
        ra2000 = 030.20313477
        dec2000 = -12.87498346

        # (Mid-)Time of observation
        jd = 2450528.2335

        # Calculate barycentric correction (debug=True show
        # various intermediate results)
        corr, hjd = pyasl.helcorr(longitude, latitude, altitude,
                                  ra2000, dec2000, jd, debug=True)

        print("Barycentric correction [km/s]: ", corr)
        print("Heliocentric Julian day: ", hjd)


class SanityOfSVD(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_deltaExample(self):
        """
        The delta function example from SVD.
        """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy import pyasl

        # Get some "data"
        wvl = np.arange(5000., 5010., 0.02)
        template = np.ones(len(wvl))

        # There are two sharp lines in the template
        template[100] = 0.5
        template[115] = 0.3
        # Apply rotational broadening to the delta spectrum
        nflux = pyasl.rotBroad(wvl, template, 0.5, 23.45,
                               edgeHandling="firstlast")

        # Carry out decomposition
        svd = pyasl.SVD()
        svd.decompose(template, 51)
        # Obtain the broadening function needed to
        # recover "observed" spectrum. Note that the
        # edges (51/2 bins) will actually be neglected.
        b = svd.getBroadeningFunction(nflux)
        # Get the model, which results from the broadening
        # function and the template; obtain the indices
        # where it applies, too.
        m, mind = svd.getModel(b, modelIndices=True)

        # Plot the outcome
        plt.plot(b, 'bp-')
        plt.plot(mind, m, 'r.')
        plt.plot(nflux, 'g--')
#    plt.show()

    def sanity_NoiseExample(self):
        """
        The noise example from SVD
        """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy import pyasl

        # Get some "data"
        wvl = np.arange(5000., 5010., 0.02)
        template = np.ones(len(wvl))

        # There are two sharp lines in the template
        template[100] = 0.5
        template[115] = 0.3
        # Apply rotational broadening to the delta spectrum
        nflux = pyasl.rotBroad(wvl, template, 0.5, 23.45,
                               edgeHandling="firstlast")
        nflux += np.random.normal(0., 0.005, len(nflux))

        # Carry out decomposition
        svd = pyasl.SVD()
        svd.decompose(template, 51)

        # Access the singular values
        sv = svd.getSingularValues()

        # Calculate the reduced chi square as a function of the number
        # of singular values neglected in the calculation of the
        # model.
        chi = []
        for i in range(1, len(sv), 5):
            b = svd.getBroadeningFunction(nflux, wlimit=sorted(sv)[i])
            m, mind = svd.getModel(b, modelIndices=True, asarray=True)
            chi.append(((nflux[mind] - m)**2 / 0.005**2).sum() / len(mind))

        plt.title("Reduced $\chi^2$ vs. number of neglected singular values")
        plt.plot(range(1, len(sv), 5), chi, 'bp-')
#    plt.show()


class SanityOfUnred(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example(self):
        """
        Check the example
        """
        from PyAstronomy import pyasl
        import numpy as np
        import matplotlib.pylab as plt

        # Approximate a solar spectrum using a Planck
        # function with a temperature of 5778 K between
        # 3000 A and 8000 A.
        wvl = np.arange(3000., 8000., 1.0)
        flux = pyasl.planck(T=5778., lam=wvl * 1e-10)

        # Deredden the spectrum assuming ebv=0.1
        fluxUnred = pyasl.unred(wvl, flux, ebv=0.1, R_V=3.1)

        # Plot the result
        plt.title("Reddened flux (red) and dereddened flux (blue)")
        plt.plot(wvl, flux, 'r--')
        plt.plot(wvl, fluxUnred, 'b--')
#    plt.show()


class SanityOfEstimateSNR(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example(self):
        """
        Check the example for estimateSNR
        """
        from PyAstronomy import pyasl
        import numpy as np

        # Number of data points
        N = 10000
        # Signal to noise ratio
        SNR = 50.0

        # Create some data with noise and a sinusoidal
        # variation.
        x = np.arange(N)
        y = np.random.normal(0.0, 1.0 / SNR, N) + 1.0
        y += np.sin(x / 500.0 * 2. * np.pi) * 0.1

        # Estimate the signal to noise ratio. Check whether the
        # estimate fits the input...
        # Use a chunk length of 20 data points, a polynomial of degree
        # one, and produce a "control plot".
        snrEsti = pyasl.estimateSNR(x, y, 20, deg=1, controlPlot=False)
        print("Estimate of the SNR: ", snrEsti["SNR-Estimate"])

        # Use a chunks with a length of 27, a polynomial of degree
        # two, and produce a "control plot".
        snrEsti = pyasl.estimateSNR(
            x, y, 27, deg=2, controlPlot=False, xlenMode="excerpt")
        print("Estimate of the SNR: ", snrEsti["SNR-Estimate"])


class SanityOfOutlier(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_ESDexample(self):
        """
        Check the example for the generalized ESD
        """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy import pyasl

        # Convert data given at:
        # http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm
        # to array.
        x = np.array([float(x) for x in "-0.25 0.68 0.94 1.15 1.20 1.26 1.26 1.34 1.38 1.43 1.49 1.49 \
              1.55 1.56 1.58 1.65 1.69 1.70 1.76 1.77 1.81 1.91 1.94 1.96 \
              1.99 2.06 2.09 2.10 2.14 2.15 2.23 2.24 2.26 2.35 2.37 2.40 \
              2.47 2.54 2.62 2.64 2.90 2.92 2.92 2.93 3.21 3.26 3.30 3.59 \
              3.68 4.30 4.64 5.34 5.42 6.01".split()])

        # Apply the generalized ESD
        r = pyasl.generalizedESD(x, 10, 0.05, fullOutput=True)

        print("Number of outliers: ", r[0])
        print("Indices of outliers: ", r[1])
        print("        R      Lambda")
        for i in range(len(r[2])):
            print("%2d  %8.5f  %8.5f" % ((i + 1), r[2][i], r[3][i]))

        # Plot the "data"
        plt.plot(x, 'b.')
        # and mark the outliers.
        for i in range(r[0]):
            plt.plot(r[1][i], x[r[1][i]], 'rp')
#     plt.show()

    def sanity_pointDistGESDExample(self):
        """
        Sanity of pointDistGESD example
        """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy import pyasl

        # Get some data
        x = np.random.normal(0., 0.1, 50)

        # Introduce outliers
        x[27] = 1.0
        x[43] = -0.66

        # Run distance based outlier detection
        r = pyasl.pointDistGESD(x, 5)

        print("Number of outliers detected: ", r[0])
        print("Indices of these outliers: ", r[1])

        plt.plot(x, 'b.')
        for i in range(len(r[1])):
            plt.plot(r[1][i], x[r[1][i]], 'rp')
#     plt.show()

    def sanity_polyResOutlierExample(self):
        """
        Sanity of polyResOutlier example
        """
        from PyAstronomy import pyasl
        import numpy as np
        import matplotlib.pylab as plt

        # Generate some "data"
        x = np.arange(100)
        y = np.random.normal(x * 0.067, 1.0, len(x))

        # Introduce an outliers
        y[14] = -5.0
        y[67] = +9.8

        # Find outliers based on a linear (deg = 1) fit.
        # Assign outlier status to all points deviating by
        # more than 3.0 standard deviations from the fit,
        # and show a control plot.
        iin, iout = pyasl.polyResOutlier(
            x, y, deg=1, stdlim=3.0, controlPlot=False)

        # What about the outliers
        print("Number of outliers: ", len(iout))
        print("Indices of outliers: ", iout)

        # Remove outliers
        xnew, ynew = x[iin], y[iin]

        # Plot result (outlier in red)
        plt.plot(x, y, 'r.')
        plt.plot(xnew, ynew, 'bp')
#     plt.show()

    def sanity_slidingPolyResOutlier(self):
        """
        Checking slidingPolyResOutlier
        """
        import numpy as np
        from PyAstronomy import pyasl

        N = 50
        x = np.arange(N)
        y = np.random.normal(0., 1., N)

        y[5] = 30
        y[6] = 10

        iin, iout = pyasl.slidingPolyResOutlier(x, y, 15, controlPlot=False)

        self.assertEqual(len(
            iout), 2, msg="slidingPolyResOutlier: Wrong number of outliers detected (" + str(len(iout)) + ").")
        self.assertTrue(np.all(iout == np.array(
            [5, 6])), msg="slidingPolyResOutlier: Wrong outliers detected (indices " + str(iout) + ").")

    def sanity_slidingPolyResOutlierExample(self):
        """
        Checking sliding PolyResOutlierExample
        """
        from PyAstronomy import pyasl
        import numpy as np
        import matplotlib.pylab as plt

        # Generate some "data"
        x = np.arange(100)
        y = np.random.normal(x * 0.067 + 0.01 * x**2, 1.0, len(x))

        # Introduce an outliers
        y[14] = -5.0
        y[67] = +9.8

        # Find outliers based on a linear (deg = 1) fit.
        # Assign outlier status to all points deviating by
        # more than 3.0 standard deviations from the fit,
        # and show a control plot.
        iin, iout = pyasl.slidingPolyResOutlier(
            x, y, 20, deg=1, stdlim=3.0, controlPlot=False)

        # What about the outliers
        print("Number of outliers: ", len(iout))
        print("Indices of outliers: ", iout)

        # Remove outliers
        xnew, ynew = x[iin], y[iin]

        # Plot result (outlier in red)
        plt.plot(x, y, 'r.')
        plt.plot(xnew, ynew, 'bp')
#     plt.show()

    def sanity_generalizedESD_NISTR_compatibility(self):
        """ Checking generalizedESD vs R implementation of NIST """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy import pyasl
        
        # Convert data given at:
        # http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm
        # to array.
        x = np.array([float(x) for x in "-0.25 0.68 0.94 1.15 1.20 1.26 1.26 1.34 1.38 1.43 1.49 1.49 \
                  1.55 1.56 1.58 1.65 1.69 1.70 1.76 1.77 1.81 1.91 1.94 1.96 \
                  1.99 2.06 2.09 2.10 2.14 2.15 2.23 2.24 2.26 2.35 2.37 2.40 \
                  2.47 2.54 2.62 2.64 2.90 2.92 2.92 2.93 3.21 3.26 3.30 3.59 \
                  3.68 4.30 4.64 5.34 5.42 6.01".split()])
        
        # Apply the generalized ESD
        r = pyasl.generalizedESD(x, 10, 0.05, fullOutput=True, ubvar=True)
        
        # NIST example
        # https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.r
        r1 = r[2][0] - 3.118906
        r2 = r[2][6] - 2.279327
        self.assertAlmostEqual(r1, 0.0, delta=1e-7, msg="generalizedESD incompatible with NIST's R implementation (r1 = " + str(r1) + ")")
        self.assertAlmostEqual(r2, 0.0, delta=1e-7, msg="generalizedESD incompatible with NIST's R implementation (r2 = " + str(r2) + ")")


class SanityOfMagnitudes(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_absMagToPowerExample(self):
        """
        Sanity of example for `absMagToPower`.
        """
        from PyAstronomy import pyasl

        absMagSun = 4.75
        print("Absolute bolometric magnitude of the Sun: ", absMagSun)
        print("  Absolute luminosity [erg/s]: ",
              pyasl.absMagToPower(absMagSun))

    def sanity_absModuleToDistExample(self):
        """
        Sanity of distance module example.
        """
        from PyAstronomy import pyasl

        # Apparent magnitude
        appMag = 11.37
        # Absolute (bolometric) magnitude of Sun
        absMagSun = 4.75

        print("Distance of a sun-like star with apparent bolometric ", end=' ')
        print("brightness of 11.37 mag: %5.2f pc" %
              (pyasl.absModuleToDist(appMag, absMagSun)))

    def sanity_absModuleToDist(self):
        """
        Sanity checks for absModuleToDist.
        """
        from PyAstronomy import pyasl
        x = pyasl.absModuleToDist(5.37, 5.37)
        self.assertAlmostEqual(x, 10., delta=1e-7)
        x = pyasl.absModuleToDist(5.37 + 5., 5.37)
        self.assertAlmostEqual(x, 100., delta=1e-7)


class SanityOfSunpos(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_sunposExample(self):
        """
        Sanity of example for `sunpos`.
        """
        import numpy as np
        from PyAstronomy import pyasl
        import datetime

        # Convert calendar date into JD
        # use the datetime package
        jd = datetime.datetime(2013, 4, 16)
        jd = pyasl.jdcnv(jd)
        print("JD = " + str(jd))
        pos = pyasl.sunpos(jd, full_output=True)
        print("Coordinates of the Sun (ra, dec): %g, %g" % (pos[1], pos[2]))
        print("Solar elongation = %g and obliquity = %g" % (pos[3], pos[4]))

        # Get the Sun's RA and DEC values for a period of time.
        startjd = datetime.datetime(2013, 4, 16)
        endjd = datetime.datetime(2013, 6, 16)
        # Convert into Julian dates
        startjd = pyasl.jdcnv(startjd)
        endjd = pyasl.jdcnv(endjd)
        print()
        pos = pyasl.sunpos(startjd, end_jd=endjd, jd_steps=10,
                           plot=False, full_output=True)

        for i in range(len(pos[0])):
            print("At JD = %g: ra = %g, dec = %g" %
                  (pos[0][i], pos[1][i], pos[2][i]))


class SanityOfNutation(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_nutationExample(self):
        """
        Sanity of example for `nutate`
        """
        from PyAstronomy import pyasl
        import datetime
        import numpy as np

        # Convert calendar date into JD
        # use the datetime package
        jd = datetime.datetime(2013, 4, 16)
        jd = pyasl.jdcnv(jd)
        print("Nutation for the date.")
        res = pyasl.nutate(jd)
        print("JD = " + str(jd) + ", Longitude = " + str(res[0]) +
              ", Obliquity = " + str(res[1]))

        # Get nutation for an array of JDs.
        startjd = datetime.datetime(2013, 4, 16)
        endjd = datetime.datetime(2013, 6, 16)
        startjd = pyasl.jdcnv(startjd)
        endjd = pyasl.jdcnv(endjd)
        jds = np.arange(startjd, endjd, .5)
        print()
        print("Plot the results")
        res = pyasl.nutate(jds, plot=False)

        print("Longitude: ", res[0])
        print("Obliquity: ", res[1])

    def sanity_conutateExample(self):
        """
        Sanity of example for `co_nutate`.
        """
        from PyAstronomy import pyasl
        import datetime
        import numpy as np

        # Convert calendar date into JD
        # use the datetime package
        jd = datetime.datetime(2013, 4, 16)
        jd = pyasl.jdcnv(jd)
        # Specify RA and DEC (degrees)
        ra = 10.
        dec = 30.
        print("Get change in RA and DEC due to Earth's nutation for JD = "
              + str(jd))
        print(pyasl.co_nutate(jd, ra, dec))

        print()
        print("Get change for several RAs and DECs for the same JD")
        ra = np.arange(0., 160., 20.)
        dec = np.arange(-80., 80., 20.)
        res = pyasl.co_nutate(np.repeat(jd, ra.size), ra, dec)
        print(res[0], res[1])

        print()
        print("Get change for several RAs and DECs for different JDs")
        jds = np.arange(jd, jd + ra.size, 1)
        res = pyasl.co_nutate(jds, ra, dec)
        print("JD             delta(RA)   delta(DEC)")
        for i in range(ra.size):
            print("%12.5f   %8.5f   %8.5f" % (jds[i], res[0][i], res[1][i]))


class SanityOfAberration(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_aberrationExample(self):
        """
        Sanity of example for `aberration`.
        """
        from PyAstronomy import pyasl
        import datetime
        import numpy as np

        # Convert calendar date to JD
        # use the datetime package
        jd = datetime.datetime(2013, 4, 16)
        jd = pyasl.jdcnv(jd)
        # Specify RA and DEC
        ra = 10.
        dec = 30.
        print("Get change in RA and DEC due to annual aberration")
        print("  for JD = " + str(jd) + ":",
              np.ravel(pyasl.co_aberration(jd, ra, dec)))

        print()
        print("Get change for several RAs and DECs for the same JD")
        ra = np.arange(10., 50., 10.)
        dec = np.arange(30., 70., 10.)
        res = pyasl.co_aberration(np.repeat(jd, ra.size), ra, dec)
        print(res[0], res[1])

        print()
        print("Get change for several RAs and DECs for different JDs")
        jds = np.arange(jd, jd + ra.size, 1)
        res = pyasl.co_aberration(jds, ra, dec)
        print("JD             delta(RA)   delta(DEC)")
        for i in range(ra.size):
            print("%12.5f   %8.5f   %8.5f" % (jds[i], res[0][i], res[1][i]))


class SanityOfAltitude(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_altitudeExample(self):
        """
        Sanity of example for `altitude`
        """
        from PyAstronomy import pyasl
        import datetime
        import numpy as np

        # Assume, a star is observed at an altitude of 50 degrees
        alt = 50.
        # Now one wants to know the real altitude of the star, i.e.,
        # the altitude corrected for atmospheric refraction.
        print()
        print("Get apparent (real) altitude of a star with observed altitude of " +
              str(alt) + " degrees")
        print("  ->  Apparent altitude = ",
              alt - pyasl.co_refract_forward(alt))

        print()
        print("You are not observing from sea level, but from an altitude of 5000 meter.")
        print(("Apparent altitude = %9.5f, estimated pressure [mbar] = %9.5f, " +
               "estimated temperature [K] = %9.5f") %
              pyasl.co_refract(alt, observer_alt=5000, convert_to_observed=False))

        print()
        print("Convert apparent (real) altitude into observed altitude.")
        print("Apparent altitude = " + str(alt) + " degrees", end=' ')
        print(" -> Observed altitude = " + str(pyasl.co_refract(alt, full_output=False,
                                                                convert_to_observed=True)[0]))

        print()
        print("The same object observed from different observer altitudes")
        apparentAltitudes = np.repeat(30.0, 10)
        obsalts = np.linspace(0., 5000., len(apparentAltitudes))
        r = pyasl.co_refract(
            apparentAltitudes, observer_alt=obsalts, convert_to_observed=True)
        for i in range(len(r[0])):
            print("Observed altitude [deg] = %g, pressure [mbar] = %g, temperature [K] = %g"
                  % (r[0][i], r[1][i], r[2][i]))


class SanityOfHorizontalCoordinates(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_hadec2altazExample(self):
        """
        Sanity of example for `hadec2altaz`
        """
        from PyAstronomy import pyasl
        import datetime
        import numpy as np

        # Hour angle 0. means transiting the local meridian.
        ha = 0.
        # Declination of object
        dec = 30.
        # Latitude of the observer (here Hamburger Sternwarte)
        lat = +53.48
        print("Get altitude and azimuth of object in DEGREES")
        print(pyasl.hadec2altaz(ha, dec, lat))

        # List of coordinates
        ha = np.arange(0., 20., 5.)
        dec = np.arange(30., 50., 5.)
        lat = np.zeros(dec.size) + 53.48
        print()
        print("Get altitude and azimuth for a list of objects from same observer latitude")
        altaz = pyasl.hadec2altaz(ha, dec, lat)
        print("alt: ", altaz[0])
        print("az: ", altaz[1])

    def sanity_eq2horExample(self):
        """
        Sanity of example for `eq2hor`
        """
        from PyAstronomy import pyasl
        import matplotlib.pylab as plt
        import datetime
        import numpy as np

        # Convert calendar date to JD
        # use the datetime package
        jd = datetime.datetime(2013, 4, 16)
        jd = pyasl.jdcnv(jd)
        # Specific RA and DEC
        ra = 10.
        dec = 30.
        print()
        print("Get horizontal coordinates (alt, az, ha) from JD, RA,")
        print("  and DEC for the Hamburger Sternwarte")
        print(pyasl.eq2hor(jd, ra, dec, observatory="HS"))

        print()
        print("From a list of Julian dates ...")
        jds = np.arange(jd, jd + 1, .2)
        ras = np.zeros(jds.size) + ra
        decs = np.zeros(jds.size) + dec
        alt, az, ha = pyasl.eq2hor(
            jds, ras, decs, lon=-70.4042, lat=-24.6272, alt=2635.)

        for i in range(alt.size):
            print("JD = %g : alt = % g,  az = % g,  ha = % g" %
                  (jds[i], alt[i], az[i], ha[i]))

        print()
        print("For one object and different times at the VLT...")
        jds = np.arange(jd - .25, jd + .25, .01)
        ras = np.zeros(jds.size) + 130.
        decs = np.zeros(jds.size) - 30.
        res = pyasl.eq2hor(jds, ras, decs, lon=-70.4042,
                           lat=-24.6272, alt=2635.)

        plt.plot(jds, res[0])
        plt.xlabel("Julian date")
        plt.ylabel("Altitude [deg]")
#     plt.show()


class SanityOfObservatory(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example(self):
        """
        Sanity of observatory example.
        """
        from PyAstronomy import pyasl

        # List all available observatory data
        pyasl.listObservatories()

        print()
        print("Data for Kitt Peak National Observatory")
        print(pyasl.observatory("kpno"))
        print("(longitude and latitude in degrees, altitude in meters, and")
        print("time zone in hours West of Greenwich")


class SanityOfDampingConstConversion(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example(self):
        """
        Checking example for converting damping constant into line width.
        """
        from PyAstronomy import pyasl

        # Einstein coefficient relevant for hydrogen LyA
        gLya = 6.258085e8

        print("Width of H LyA line at 1215.67 A = %e cm" %
              pyasl.convertDampingConstant(gLya, 1215.67))


class SanityOfCardinalPoint(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_checkNumbers(self):
        """
        Check some numbers for cardinal point
        """
        self.assertEqual(getCardinalPoint(0.0), "N")
        self.assertEqual(getCardinalPoint(360.0), "N")
        self.assertEqual(getCardinalPoint(135.0), "E")
        self.assertEqual(getCardinalPoint(45.0), "N")
        self.assertEqual(getCardinalPoint(225.0), "S")
        self.assertEqual(getCardinalPoint(315.0), "W")
        self.assertEqual(getCardinalPoint(217.0), "S")

    def sanity_checkExample(self):
        """
        Check cardinal point example
        """
        from PyAstronomy import pyasl
        import numpy as np

        # Get the cardinal point for 10 azimuth angles
        azimuths = np.random.random(10) * 360.
        for azimuth in azimuths:
            cp = pyasl.getCardinalPoint(azimuth)
            print("Azimuth: {0:6.2f} deg, Cardinal point: {1:1s}".format(
                azimuth, cp))


class SanityOfTwilight(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_checkNumbers(self):
        """
        Check twilightName example
        """
        from PyAstronomy import pyasl
        import numpy as np

        for alt in np.linspace(-20., 5., 15):
            print("Altitude = {0:6.2f}, Twilight is called: ".format(alt),
                  pyasl.twilightName(alt))


class SanityOfAngDist(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_checkExample(self):
        """
        Check angular distance example
        """
        from PyAstronomy import pyasl

        print("Angular distance between the poles (deg):")
        print(pyasl.getAngDist(98.0, -90.0, 100., +90.0))

        print("Angular distance between Vega and Altair (deg)")
        print(pyasl.getAngDist(279.23473479, +
                               38.78368896, 297.69582730, +08.86832120))


class SanityOfTransit(unittest.TestCase, SaniBase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_transitDurationExample(self):
        """
        Transit duration example
        """
        from PyAstronomy import pyasl
        from PyAstronomy import constants as pc

        # Earth radius expressed in Jovian radii
        reJ = pc.REarth / pc.RJ
        print("Earth radius in Jovian units: ", reJ)

        # Estimate the duration of Earth's transit
        td = pyasl.transitDuration(1.0, reJ, 1.0, 90.0, 365.0)
        print("The transit of Earth lasts about: %5.3f days" % td)
        
        # Estimate the duration of in- and egress
        ti = pyasl.ingressDuration(1.0, reJ, 1.0, 90.0, 365.0)
        print("The in- and egress last: %6.4f days or %4.2f hours" % (ti, ti*24))
    
    def sanity_transitDuration_RsExample(self):
        """
        Check example of transitDuration_Rs
        """
        from PyAstronomy import pyasl
        
        # Semi-major axis in units of stellar radius
        sma = 8.8
        # Radius ratio (Rp/Rs)
        rprs = 0.16
        
        # Estimate the duration of Earth's transit
        td = pyasl.transitDuration_Rs(sma, rprs, 85.7, 2.21858)
        print("The transit of HD 189733 b lasts about: %5.2f hours" % (td*24.))
    
        # Estimate the duration of in- and egress
        ti = pyasl.ingressDuration_Rs(sma, rprs, 85.7, 2.21858)
        print("The in- and egress of HD 189733 b lasts: %5.2f hours" % (ti*24.))
    
    def sanity_transitDuration(self):
        """
        """
        from PyAstronomy import pyasl
        from PyAstronomy import constants as PC
        
        td1 = pyasl.transitDuration(0.05, 0.1, 1.5, 88.0, 10.)
        
        pc = PC.PyAConstants()
        rs = 1.5*pc.RSun
        
        trs = pyasl.transitDuration_Rs(0.05*pc.AU/rs, 0.1*pc.RJ/rs, 88., 10.)
        
        print("td1: ", td1, ", trs: ", trs)
        self.assertAlmostEqual(td1, trs, delta=1e-10, msg="transitDuration and transitDuration_Rs in conflict!")
    
    def sanity_ingressDuration(self):
        """
        Check sanity and consistency of ingressDuration
        """
        from PyAstronomy import pyasl
        import numpy as np
        from PyAstronomy import constants as pc
        
        rp = pc.RSun / pc.RJ * 0.4
        rprs = 0.4
        
        for sma in [0.2, 0.6, 0.8]:
            for P in [2., 60., 300.]:
                t = pyasl.ingressDuration(sma, rp, 1.0, 90, P)
                # Linear approximation
                app = rp*pc.RJ*P/(sma*pc.AU*np.pi)
                print(sma, P, (t-app)/P)
                self.assertAlmostEqual((t-app)/P, 0.0, delta=1e-6, msg="Ingress duration for 90 deg case does not match approximation.")
                t2 = pyasl.ingressDuration_Rs(sma*pc.AU/pc.RSun, rprs, 90, P)
                print(t-t2)
                self.assertAlmostEqual(t, t2, delta=1e-12, msg="ingressDuration and ingressDuration_Rs differ.") 
        
        print()
        for sma in [0.05, 0.1, 0.15]:
            for P in [2., 60., 300.]:
                t = pyasl.ingressDuration(sma, rp, 1.0, 89., P)
                app = rp*pc.RJ*P/(sma*pc.AU*np.pi)
                print(sma, P, t, app, (t-app)/P)
                self.assertAlmostEqual((t-app)/P, 0.0, delta=2e-3, msg="Ingress duration for 89 deg case does not match approximation.")
        
        
            
    def sanity_transitDurationExact(self):
        """
        Check approximation of transit duration.
        """
        from PyAstronomy import pyasl
        import numpy as np
        td1 = pyasl.transitDuration(0.05, 0.1, 1.0, 90.0, 365.0, exact=False)
        td2 = pyasl.transitDuration(0.05, 0.1, 1.0, 90.0, 365.0, exact=True)
        self.assertEqual(td1, td2, msg="Exact and approximate transit duration differ by " + str(td1-td2) + " at inclination of 90 deg.")
        td1 = pyasl.transitDuration(0.05, 0.1, 1.0, 88.0, 365.0, exact=False)
        td2 = pyasl.transitDuration(0.05, 0.1, 1.0, 88.0, 365.0, exact=True)
        self.assertAlmostEqual(365./np.pi*np.arcsin(np.sin(td1*np.pi/365.)/np.sin(88./180.*np.pi)), td2, delta=1e-12, msg="Relation between exact and approximate transit duration broken.")

    def sanity_inTransitExample_1(self):
        """
        inTransit---individual point in time
        """
        from PyAstronomy import pyasl

        # Time of interest
        time = 2476357.756234
        # Define some (arbitrary) transit parameters
        T0 = 2475123.01245
        period = 3.4789112
        duration = 2.2 / 24.0

        # Check whether the time is in-transit
        print("Time is within transit? ", end=' ')
        if not pyasl.isInTransit(time, T0, period, duration / 2.0):
            print("No")
        else:
            print("Yes")

    def sanity_inTransitExample_2(self):
        """
        inTransit---Series of points in time
        """
        from PyAstronomy import pyasl
        import numpy as np

        # Times of interest
        times = 2476357.756234 + np.linspace(0.0, 5.0, 300)
        # Define some (arbitrary) transit parameters
        T0 = 2475123.01245
        period = 3.4789112
        duration = 2.2 / 24.0

        # Check whether the time is in-transit
        print("Indices if time points within transit: ", end=' ')
        print(pyasl.isInTransit(times, T0, period, duration / 2.0))

        print()
        print("For each time point, a flag indicating whether it")
        print("is in- or off-transit:")
        print(pyasl.isInTransit(times, T0, period,
                                duration / 2.0, boolOutput=True))

    def sanity_transit_T1_T4_ell_example(self):
        """
        Checking example for transit_T1_T4_ell example
        """
        from PyAstronomy import pyasl
        
        # SMA in stellar radii
        sma = 5.67
        # Rp/Rs
        rprs = 0.15
        # Orbital inclination
        inc = 89.2
        # Orbital period (time units are arbitrary but must be consistent)
        p = 2.0
        # Eccentricity
        e = 0.63
        # Argument of periastron (planetary orbit)
        w = 155.
        # Time of periastron passage
        tau = 2412345.346 
        
        # Contact times for primary transit
        pts = pyasl.transit_T1_T4_ell(sma, rprs, inc, p, tau, e, w, transit="p")
        # Contact times for secondary transit
        sts = pyasl.transit_T1_T4_ell(sma, rprs, inc, p, tau, e, w, transit="s")
        
        print("Transit times at arbitrary epoch (N*period may added)")
        print("Primary transit T1-4: ", pts)
        print("Secondary trabnsit T1-4: ", sts)
        print()
        print("Duration of primary and secondary transit: %5.3f, %5.3f " % (pts[3]-pts[0], sts[3]-sts[0]))

    def sanity_ell_vs_circ_duration(self):
        """
        Checking sanity of T1-4 elliptical contact points vs. circular solution
        """
        from PyAstronomy import pyasl
        import numpy as np

        np.random.seed(23847746)
        
        for i in range(20):
        
            sma = 3 + 10*np.random.random()
            rprs = 0.01+0.2*np.random.random()
            inc = 87. + 3*np.random.random()
            p = 1 + 10*np.random.random()
        
            a = pyasl.transitDuration_Rs(sma, rprs, inc, p, exact=True)
            #print(a)
            d = pyasl.ingressDuration_Rs(sma, rprs, inc, p)
            #print(d)
        
            b = pyasl.transit_T1_T4_ell(sma, rprs, inc, p, 0.0, 0.0, -90, transit="p")
            t21 = b[1] - b[0]
            t43 = b[3] - b[2]
            
            bs = pyasl.transit_T1_T4_ell(sma, rprs, inc, p, 0.0, 0.0, -90, transit="s")
            db = [ (b[j]-bs[j]) % p - p/2 for j in range(4)]
            
            self.assertAlmostEqual(a - (b[3]-b[0]), 0, delta=1e-5, msg="Problem with transit duration at iteration %d" % i)
            self.assertAlmostEqual(d - (b[1]-b[0]), 0, delta=1e-5, msg="Problem with ingress duration at iteration %d" % i)
            
            x = np.max(np.abs(db))
            self.assertAlmostEqual(x, 0, delta=1e-5, msg="Problem with offset between primary and secondary transit at iteration %d" % i)

        
        np.random.seed(23847746+1)
        
        for i in range(10):
        
            sma = 3 + 10*np.random.random()
            rprs = 0.01+0.2*np.random.random()
            inc = 87. + 3*np.random.random()
            p = 1 + 10*np.random.random()
            e = np.random.random() * 0.6
            w = 360 * np.random.random()
            tau = 1e3 * np.random.random()
            
            print(i,sma, rprs, inc, p, tau, e, w)
        
            b = pyasl.transit_T1_T4_ell(sma, rprs, inc, p, tau, e, w, transit="p")
            
            b = [np.NaN if b[j] is None else b[j] for j in range(4)]
            #print(b)
            
            bim = sma * np.cos(inc/180.*np.pi)
            da = 1/np.pi * p / sma * np.sqrt((1+rprs)**2-bim**2) * np.sqrt(1-e**2) / (1+e*np.sin((w+180)/180.*np.pi))
         
            self.assertAlmostEqual((b[3]-b[0]-da)/(b[3]-b[0]), 0.0, delta=0.1, msg="Problem with transit duration (elliptical) at iteration %d" % i)

            k = e*np.sin((w+180)/180.*np.pi)
            tt = (1+k)/(1-k)
            
            bs = pyasl.transit_T1_T4_ell(sma, rprs, inc, p, tau, e, w, transit="s")
            bs = [np.NaN if bs[j] is None else bs[j] for j in range(4)]
            
            self.assertAlmostEqual((tt-(bs[3]-bs[0])/(b[3]-b[0]))/tt, 0.0, delta=0.1, msg="Problem with transit duration ratio (elliptical) at iteration %d" % i)
            
        

class SanityOfAirmass(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_airmassSphericalExample(self):
        """
        Example airmass spherical, homogeneous atmosphere.
        """
        from PyAstronomy import pyasl

        obsAlt = 2400.0

        for za in range(0, 90, 10):
            print("Zenith angle: %2d, airmass = %7.2f" %
                  (za, pyasl.airmassSpherical(za, obsAlt)))

    def sanity_airmassPPExample(self):
        """
        Example for plane-parallel airmass
        """
        from PyAstronomy import pyasl

        print("Airmass for plane-parallel atmosphere")
        for za in range(0, 70, 10):
            print("Zenith angle: %2d deg, airmass = %7.2f" %
                  (za, pyasl.airmassPP(za)))

    def sanity_compAirmasses(self):
        """
        Compare different airmass calculations
        """
        from PyAstronomy import pyasl
        for za in range(0, 70, 10):
            ampp = pyasl.airmassPP(za)
            amsp = pyasl.airmassSpherical(za, 0.0)
            self.assertAlmostEqual(ampp / amsp, 1.0, delta=0.01)


class SanityOfTransitTimes(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example(self):
        """
        Checking example of transitTimes
        """
        from PyAstronomy import pyasl
        import datetime as dt

        # Get the data for WASP-7 from NEXA data base
        nexa = pyasl.NasaExoplanetArchive()
        dat = nexa.selectByPlanetName("Wasp-7 b")

        # Convert July 14th, 2018, 10pm to a Julian date
        d = dt.datetime(2018, 7, 14, 22)
        jd = pyasl.jdcnv(d)

        # Calculate transit data for transits within 100 days starting
        # form July 14th, 2018.
        dat = pyasl.transitTimes(jd, jd + 100., dat, nexaInput=True,
                                 observatory="esoparanal", obsOffset=1. / 24.,
                                 minAltitude=10.0)
        # Plot the result
#    pyasl.transitVisibilityPlot(dat, markTransit=True)

    def sanity_ETD_compare(self):
        """
        Check transit of Wasp-7 on Oct. 03, 2018, with ETD result.
        """
        from PyAstronomy import pyasl
        import datetime as dt

        # Get the data for WASP-7 from NEXA data base
        nexa = pyasl.NasaExoplanetArchive()
        dat = nexa.selectByPlanetName("Wasp-7 b")

        # These are the numbers from Hellier 2008. Meanwhile, NEXA uses
        # Southworth 2011. Keeping the 'old' numbers ensures that results remain
        # comparable to the numbers provided by ETC used for cross checking. 
        dat["pl_tranmid"] = 2453985.0149
        dat["pl_orbper"] = 4.954658
        dat["pl_radj"] = 0.915
        dat["st_rad"] = 1.236
        dat["pl_orbincl"] = 89.6
        dat["pl_orbsmax"] = 0.0618

        # Convert July 14th, 2018, 10pm to a Julian date
        d = dt.datetime(2018, 7, 14, 22)
        jd = pyasl.jdcnv(d)

        # Calculate transit data for transits within 100 days starting
        # form July 14th, 2018.
        dat = pyasl.transitTimes(jd, jd + 100., dat, nexaInput=True,
                                 observatory="esoparanal", obsOffset=0.0,
                                 minAltitude=10.0)

        # Choose first transit
        d = dat[2]
        print(d)
        mt = d["Obs cal"][1]
        self.assertEqual(mt[0], 2018, "Year does not match")
        self.assertEqual(mt[1], 10, "Month does not match")
        self.assertEqual(mt[2], 3, "Day does not match")
        self.assertAlmostEqual(mt[3], 3. + 53. / 60.,
                               delta=5. / 60., msg="Hour does not match")
        self.assertAlmostEqual(
            d["Star alt"][0], 67., msg="Stellar altitude does not match", delta=3.0)
        self.assertAlmostEqual(
            d["Star alt"][1], 47., msg="Stellar altitude does not match", delta=3.0)
        self.assertAlmostEqual(
            d["Star alt"][2], 26., msg="Stellar altitude does not match", delta=3.0)
        # Transit duration is specified in Nexa (therefore, 15 min uncertainty here)
        self.assertAlmostEqual((d["Transit jd"][2] - d["Transit jd"][0]) * 24.0,
                               3.567, delta=15. / 60., msg="Duration does not match")


class SanityOfDecimalYear(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example(self):
        """
        Check example of decimal year.
        """
        import datetime as dt
        from PyAstronomy import pyasl

        # Convert July 2nd, 1998, 12:30:59 into decimal
        # representation
        d = dt.datetime(1998, 7, 2, 12, 30, 59)

        # Obtain float representation of decimal year
        decy = pyasl.decimalYear(d)
        print("Decimal representation: ", decy)

        # Convert back into gregorian date first
        print("The decimal year %10.5f correspond to " % decy +
              pyasl.decimalYearGregorianDate(decy, "yyyy-mm-dd hh:mm:ss"))
        print(" ... or equivalently (y, m, d, h, m, s, ms): ",
              pyasl.decimalYearGregorianDate(decy, "tuple"))

    def sanity_decimalYear(self):
        """
        Decimal year: Check that conversion is reasonable.
        """
        import datetime as dt
        from PyAstronomy import pyasl
        d = dt.datetime(1998, 7, 2, 12, 30, 59)
        # Obtain float representation of decimal year
        decy = pyasl.decimalYear(d)
        t = pyasl.decimalYearGregorianDate(decy, "tuple")

        self.assertEqual(t[0], 1998, msg="Year does not match")
        self.assertEqual(t[1], 7, msg="Month does not match")
        self.assertEqual(t[2], 2, msg="Day does not match")
        self.assertEqual(t[3], 12, msg="Hour does not match")
        self.assertEqual(t[4], 30, msg="Minute does not match")
        self.assertEqual(t[5], 59, msg="Second does not match")

        t = pyasl.decimalYearGregorianDate(2008.8196, "tuple")


class SanityOfBroad(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example(self):
        """
        Checking sanity of 'instrBroadGaussFast' example.
        """
        from PyAstronomy import pyasl
        import matplotlib.pylab as plt
        import numpy as np

        # Set up an input spectrum
        x = np.linspace(5000.0, 5100.0, 1000)
        y = np.ones(x.size)

        # Introduce some delta-peaked lines
        y[165] = 0.7
        y[187] = 0.3
        y[505] = 0.1
        y[610] = 0.1
        y[615] = 0.7

        # Apply Gaussian instrumental broadening, setting the resolution to 10000.
        r, fwhm = pyasl.instrBroadGaussFast(x, y, 10000,
                                            edgeHandling="firstlast", fullout=True)

        # Apply Gaussian instrumental broadening, setting the resolution to 10000.
        # Limit the extent of the Gaussian broadening kernel to five standard
        # deviations.
        r2, fwhm = pyasl.instrBroadGaussFast(x, y, 10000,
                                             edgeHandling="firstlast", fullout=True, maxsig=5.0)

        print("FWHM used for the Gaussian kernel: ", fwhm, " A")

        # Plot the output
        plt.plot(x, r, 'r--p', label="Broadened curve (full)")
        plt.plot(x, r2, 'k:', label="Broadened curve (5 stds)")
        plt.plot(x, y, 'b-', label="Input")
        plt.legend(loc=4)
#     plt.show()

    def sanity_tests(self):
        """
        Checking sanity of 'instrBroadGaussFast'
        """
        from PyAstronomy import pyasl
        import numpy as np

        # Set up an input spectrum
        x = np.linspace(5000.0, 5100.0, 1000)
        y = np.zeros(x.size)

        # Introduce some delta-peaked lines
        y[165] = 0.7
        y[187] = 0.3
        y[505] = 0.1
        y[610] = 0.1
        y[615] = 0.7

        # Apply Gaussian instrumental broadening.
        resolution = 12700.
        r, fwhm = pyasl.instrBroadGaussFast(
            x, y, resolution, edgeHandling=None, fullout=True)

        self.assertAlmostEqual(fwhm, np.mean(
            x) / resolution, 6, "FWHM does not match")

        s1 = y.sum()
        s2 = r.sum()
        self.assertAlmostEqual(s1, s2, 6, "EW in spectrum did change")

    def sanity_convolutionGaussian(self):
        """
        Check the sanity of broadGaussFast (properties of result in case of two Gaussians)
        """
        import numpy as np
        from PyAstronomy import pyasl
        from PyAstronomy import funcFit as fuf

        x = np.arange(-8.0, 8.0, 0.001)
        gf = fuf.GaussFit1d()
        gf["A"] = 0.8
        gf["sig"] = 0.471
        y = gf.evaluate(x)

        y2 = pyasl.broadGaussFast(
            x, y, 1.5 - gf["sig"], edgeHandling="firstlast")

        print(np.mean(x * y) / np.sum(y), np.mean(x * y2) / np.sum(y2))
        print(np.sqrt(np.sum(x**2 * y) / np.sum(y)))
        print(np.sqrt(np.sum(x**2 * y2) / np.sum(y2)))
        print(np.sum(y), np.sum(y2))

        mean2 = np.mean(x * y2) / np.sum(y2)
        self.assertAlmostEqual(mean2, 0.0, delta=1e-9,
                               msg="Barycenter of convoluted Gaussian (" + str(mean2) + ") deviates from 0.0.")
        std2 = np.sqrt(np.sum(x**2 * y2) / np.sum(y2))
        std2_nom = np.sqrt(gf["sig"]**2 + (1.5 - gf["sig"])**2)
        self.assertAlmostEqual(std2, std2_nom, delta=1e-5,
                               msg="Std of convoluted Gaussian (" + str(std2) + ") deviates from nominal value.")
        self.assertAlmostEqual(np.sum(y), np.sum(y2), delta=1e-4,
                               msg="Normalization of convoluted Gaussian (" + str(np.sum(y2)) + ") is incorrect.")

    def sanity_convolutionGaussianMaxsig(self):
        """
        Check the sanity of broadGaussFast including maxsig
        """
        import numpy as np
        from PyAstronomy import pyasl
        from PyAstronomy import funcFit as fuf

        x = np.arange(-10.0, 10.0, 0.001)
        gf = fuf.GaussFit1d()
        gf["A"] = 0.8
        gf["sig"] = 0.471
        y = gf.evaluate(x)

        y3 = pyasl.broadGaussFast(
            x, y, 1.5 - gf["sig"], edgeHandling="firstlast", maxsig=5.0)
        y2 = pyasl.broadGaussFast(
            x, y, 1.5 - gf["sig"], edgeHandling="firstlast", maxsig=None)

        self.assertAlmostEqual(np.max(np.abs(y2 - y3)), 0.0, delta=2e-7,
                               msg="Deviation in profiles using maxsig (broadGaussFast)")

        print(np.mean(x * y) / np.sum(y), np.mean(x * y2) / np.sum(y2))
        print(np.sqrt(np.sum(x**2 * y) / np.sum(y)))
        print(np.sqrt(np.sum(x**2 * y2) / np.sum(y2)))
        print(np.sum(y), np.sum(y2))

        mean2 = np.mean(x * y2) / np.sum(y2)
        self.assertAlmostEqual(mean2, 0.0, delta=1e-9,
                               msg="Barycenter of convoluted Gaussian (" + str(mean2) + ") deviates from 0.0.")
        std2 = np.sqrt(np.sum(x**2 * y2) / np.sum(y2))
        std2_nom = np.sqrt(gf["sig"]**2 + (1.5 - gf["sig"])**2)
        self.assertAlmostEqual(std2, std2_nom, delta=1e-5,
                               msg="Std of convoluted Gaussian (" + str(std2) + ") deviates from nominal value.")
        self.assertAlmostEqual(np.sum(y), np.sum(y2), delta=1e-4,
                               msg="Normalization of convoluted Gaussian (" + str(np.sum(y2)) + ") is incorrect.")

    def sanity_thermalBroadening(self):
        """
        Check sanity of thermal broadening width.
        """
        from PyAstronomy import pyasl
        import numpy as np
        fwhm = pyasl.thermalBroadeningWidth(4200.0, 7500.)
        sig = pyasl.thermalBroadeningWidth(4200.0, 7500., fwhm=False)
        self.assertAlmostEqual(fwhm, 0.26, delta=1e-2,
                               msg="FWHM of thermal broadening (" + str(fwhm) + ") is incorrect")
        self.assertAlmostEqual(fwhm, sig * 2 * np.sqrt(2. * np.log(2.0)), delta=1e-9,
                               msg="FWHM and std differ.")

    def sanity_thermalExample(self):
        """
        Checking example of thermal broadening
        """
        from PyAstronomy import pyasl
        
        w0 = 6564.0
        T = 9567.0
        
        linefwhm = pyasl.thermalBroadeningWidth(w0, T)
        tbroad = pyasl.tempFromthermalBroadeningWidth(w0, linefwhm, awidth=0.0)
        
        print("Line width [FWHM]: %5.2f" % linefwhm)
        print("Thermal broadening temperature: %6.1f" % tbroad)
        
    def sanity_thermalBroadeningAndReverse(self):
        """
        Checking thermal broadening and its reverse
        """
        from PyAstronomy import pyasl

        for w in [1000, 5000, 10000]:
            for T in [100., 1000., 10000.]:
        
                linefwhm = pyasl.thermalBroadeningWidth(w, T)
                tbroad = pyasl.tempFromthermalBroadeningWidth(w, linefwhm, awidth=0.0)
                self.assertAlmostEqual(T, tbroad, delta=1e-8, msg="Input and output temperature differ in thermal broadening.")
                tbroad2 = pyasl.tempFromthermalBroadeningWidth(w, linefwhm, awidth=linefwhm/2.)
                self.assertAlmostEqual(T*0.75, tbroad2, delta=1e-8, msg="Input and output temperature differ in thermal broadening width awidth=0.5*linewidth.")



class SanityOfCrosscor(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_example(self):
        """
        Checking example of crosscorr
        """
        from PyAstronomy import pyasl
        import numpy as np
        import matplotlib.pylab as plt

        # Create the template
        tw = np.linspace(5000, 5010, 1000)
        tf = np.exp(-(tw - 5004.0)**2 / (2. * 0.1**2))

        # Create data, which are not that well sampled
        dw = np.linspace(5000, 5010, 200)
        df = np.exp(-(dw - 5004.17)**2 / (2. * 0.1**2))

        # Plot template and data
        plt.title("Template (blue) and data (red)")
        plt.plot(tw, tf, 'b.-')
        plt.plot(dw, df, 'r.-')
#     plt.show()

        # Carry out the cross-correlation.
        # The RV-range is -30 - +30 km/s in steps of 0.6 km/s.
        # The first and last 20 points of the data are skipped.
        rv, cc = pyasl.crosscorrRV(
            dw, df, tw, tf, -30., 30., 30. / 50., skipedge=20)

        # Find the index of maximum cross-correlation function
        maxind = np.argmax(cc)

        print("Cross-correlation function is maximized at dRV = ",
              rv[maxind], " km/s")
        if rv[maxind] > 0.0:
            print("  A red-shift with respect to the template")
        else:
            print("  A blue-shift with respect to the template")

        plt.plot(rv, cc, 'bp-')
        plt.plot(rv[maxind], cc[maxind], 'ro')
#     plt.show()

    def sanity_ShiftedGaussian(self):
        """
        Checking the shift of a single Gaussian.
        """
        from PyAstronomy import pyasl
        import numpy as np

        # Create the template
        tw = np.linspace(5000, 5010, 1000)
        tf = np.exp(-(tw - 5004.0)**2 / (2. * 0.1**2))

        # Create data, which are not that well sampled
        dw = np.linspace(5000, 5010, 200)
        df = np.exp(-(dw - 5004.0)**2 / (2. * 0.1**2))

        rv1, cc1 = pyasl.crosscorrRV(
            dw, df, tw, tf, -30., 30., 60. / 100., skipedge=20, mode="doppler")
        rv2, cc2 = pyasl.crosscorrRV(
            dw, df, tw, tf, -30., 30., 60. / 100., skipedge=20, mode="lin")

        m1 = np.argmax(cc1)
        m2 = np.argmax(cc2)

        self.assertAlmostEqual(rv1[m1], rv2[m2], delta=1e-10)

        accu = 60. / 100.
        for dwl in np.linspace(-0.2, 0.2, 20):
            # Create the template
            tw = np.linspace(5000, 5010, 1000)
            tf = np.exp(-(tw - 5004.0)**2 / (2. * 0.1**2))

            # Create data, which are not that well sampled
            dw = np.linspace(5000, 5010, 200)
            df = np.exp(-(dw - (5004.0 + dwl))**2 / (2. * 0.1**2))

            rv1, cc1 = pyasl.crosscorrRV(
                dw, df, tw, tf, -30., 30., 60. / 100., skipedge=20, mode="doppler")
            rv2, cc2 = pyasl.crosscorrRV(
                dw, df, tw, tf, -30., 30., 60. / 100., skipedge=20, mode="lin")

            m1 = np.argmax(cc1)
            m2 = np.argmax(cc2)

            rv = dwl / np.mean(tw) * 299792.458

            self.assertAlmostEqual(rv1[m1], rv2[m2], delta=accu)
            self.assertAlmostEqual(rv1[m1], rv, delta=accu)
            self.assertAlmostEqual(rv2[m1], rv, delta=accu)

    def sanity_randomDIstribution(self):
        """
        Checking the shift with random numbers.
        """
        from PyAstronomy import pyasl
        import numpy as np

        # Create the template
        tw = np.linspace(10000, 10010, 100)
        tf = np.random.normal(0.0, 1.0, len(tw))

        dw = tw
        df = np.roll(tf, 2)

        rv1, cc1 = pyasl.crosscorrRV(
            dw, df, tw, tf, -30., 30., 60. / 1000., skipedge=20, mode="doppler")

        m1 = np.argmax(cc1)
        rv = (tw[2] - tw[0]) / np.mean(tw) * 299792.458

        accu = 60. / 1000.
        self.assertAlmostEqual(rv1[m1], rv, delta=accu)


class SanityOfQuadExtreme(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_quadextreme(self):
        """
        Checking quadExtreme: find max and min of parabola.
        """
        import numpy as np
        from PyAstronomy import pyasl

        x = np.arange(100.0)

        sign = {"min": +1, "max": -1}

        for mode in ["min", "max"]:

            for delta in np.arange(30, 70, np.pi / 2.):
                y = sign[mode] * (x - delta)**2
                p, indi = pyasl.quadExtreme(x, y, mode=mode)
                self.assertAlmostEqual(p - delta, 0.0, delta=1e-14)
                p, indi = pyasl.quadExtreme(x, y, mode=mode, dp=(3, 7))
                self.assertAlmostEqual(p - delta, 0.0, delta=1e-14)

    def sanity_quadextremeFullPoint(self):
        """
        Checking quadExtreme: find max and min of parabola (ordinate values included).
        """
        import numpy as np
        from PyAstronomy import pyasl

        x = np.arange(100.0)

        sign = {"min": +1, "max": -1}

        for mode in ["min", "max"]:

            for delta in np.arange(30, 70, np.pi / 2.):
                c = (np.random.random() - 0.5) * 500.0
                y = sign[mode] * (x - delta)**2 + c
                p, indi = pyasl.quadExtreme(x, y, mode=mode, fullPoint=True)
                self.assertAlmostEqual(p[0] - delta, 0.0, delta=1e-8)
                self.assertAlmostEqual(p[1], c, delta=1e-8)
                p, indi = pyasl.quadExtreme(
                    x, y, mode=mode, dp=(3, 7), fullPoint=True)
                self.assertAlmostEqual(p[0] - delta, 0.0, delta=1e-8)
                self.assertAlmostEqual(p[1], c, delta=1e-8)

    def sanity_quadextremeExample(self):
        """
        Checking example for quadExtreme.
        """
        import numpy as np
        from PyAstronomy import pyasl
        import matplotlib.pylab as plt

        # Create some data (a Gaussian)
        x = np.arange(100.0)
        y = np.exp(-(x - 50.2714)**2 / (2. * 5.**2))

        # Find the maximum
        epos, mi = pyasl.quadExtreme(x, y, mode="max")
        print("Maximum found at index: ", mi, ", value at maximum: ", y[mi])
        print("Maximum found by parabolic fit: ", epos)
        print()

        # Find the maximum, use a wider range for the
        # parabolic fit.
        print("Using 5 points to each side of the maximum")
        epos, mi = pyasl.quadExtreme(x, y, mode="max", dp=(5, 5))
        print("Maximum found at index: ", mi, ", value at maximum: ", y[mi])
        print("Maximum found by parabolic fit: ", epos)
        print()

        # Do as above, but get the full output
        print("Using 2 points to each side of the maximum")
        epos, mi, xb, yb, p = pyasl.quadExtreme(
            x, y, mode="max", dp=(2, 2), fullOutput=True)
        # Evaluate polynomial at a number of points.
        # Note that, internally, the x-value of the extreme point has
        # been subtracted before the fit. Therefore, we need to re-shift
        # it in the plot.
        newx = np.linspace(min(xb), max(xb), 100)
        model = np.polyval(p, newx)

        # Plot the "data"
        plt.plot(x, y, 'bp')
        # Mark the points used in the fitting (shifted, because xb is shifted)
        plt.plot(xb + x[mi], yb, 'rp')
        # Overplot the model (shifted, because xb is shifted)
        plt.plot(newx + x[mi], model, 'r--')
#     plt.show()


class SanityOfPositionAngle(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_Example(self):
        """
        Checking position angle example
        """
        from PyAstronomy import pyasl

        # Position of Mizar: 200.98141867 +54.92535197
        # Position of Alcor: 201.30640764 +54.98795966

        # Calculate position angle
        r = pyasl.positionAngle(
            200.98141867, +54.92535197, 201.30640764, +54.98795966)

        print("Position angle of Alcor (from Mizar): %4.2f deg" % r)

    def sanity_positionAngle(self):
        """
        Checking sanity of positionAngle
        """
        # Check North
        self.assertEqual(positionAngle(10.0, 10.0, 10.0, 11.0), 0.0)
        # Check South
        self.assertEqual(positionAngle(10.0, 10.0, 10.0, 9.0), 180.0)
        # Check East
        self.assertEqual(positionAngle(10.0, 0.0, 10.1, 0.0), 90.0)
        # Check West
        self.assertEqual(positionAngle(10.0, 0.0, 9.9, 0.0), 270.0)


class SanityOfSMW_RHK(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_ExampleConversion(self):
        """
        Check the S-index -> RHK conversion example
        """
        from PyAstronomy import pyasl

        ss = pyasl.SMW_RHK()

        bv = 0.8
        teff = 5100.0
        s = 0.4

        print("Convert S-index to RHK assuming a giant")
        ss.SMWtoRHK(s, teff, bv, lc="g", verbose=True)

        print()
        print()
        print("Convert S-index to RHK assuming a main-sequence star")
        ss.SMWtoRHK(s, teff, bv, lc="ms", verbose=True)

    def sanity_ExampleShowCCF(self):
        """
        Check example showing the Ccf conversion factor.
        """
        from PyAstronomy import pyasl
        import numpy as np
        import matplotlib.pylab as plt

        ss = pyasl.SMW_RHK()

        bv = np.arange(0.4, 0.9, 0.05)
        ccfn = bv * 0.0
        ccfr = bv * 0.0
        ccfrg = bv * 0.0

        for i in range(len(bv)):
            ccfn[i] = ss.log10ccfNoyes(bv[i])
            ccfr[i] = ss.log10ccfRutten(bv[i])
            ccfrg[i] = ss.log10ccfRutten(bv[i], lc="g")

#     plt.plot(bv, ccfn, 'b.-', label="Noyes")
#     plt.plot(bv, ccfr, 'r.-', label="Rutten (ms)")
#     plt.plot(bv, ccfrg, 'g.-', label="Rutten (g)")
#     plt.xlabel("B - V [mag]")
#     plt.ylabel("Ccf")
#     plt.legend()
#     plt.show()


class SanityOfAtomicNo(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_Example(self):
        """
        Check example of atomic number (AtomicNo)
        """
        from PyAstronomy import pyasl

        an = pyasl.AtomicNo()

        # Show table with all atomic numbers, elemental
        # symbols, and the names of the elements
        an.showAll()

        print(an.getElSymbol(26))

        print(an.getAtomicNo("He"))

        print(an.getElementName(25))

    def sanity_AtomicNo(self):
        """
        Check a number of associations
        """
        from PyAstronomy import pyasl

        an = pyasl.AtomicNo()

        self.assertEqual("Fe", an.getElSymbol(26))
        self.assertEqual("Ca", an.getElSymbol(20))
        self.assertEqual("H", an.getElSymbol(1))

        self.assertEqual(an.getAtomicNo("He"), 2)
        self.assertEqual(an.getAtomicNo("O"), 8)

        self.assertEqual(an.getElementName(3), "Lithium")

        print(an.getElementName(25))


class SanityOfFitsSpec(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        import os
        os.remove("test1.fits")
        os.remove("test2.fits")

    def sanity_write1dFitsSpecExample(self):
        """
        Sanity of write1dFitsSpec example
        """
        import numpy as np
        from PyAstronomy import pyasl

        # Generate a "spectrum"
        wvl = np.arange(5000., 5010., 0.01)
        flux = np.random.normal(1.0, 0.01, wvl.size)

        # Write spectrum providing wavelength array
        pyasl.write1dFitsSpec("test1.fits", flux, wvl=wvl, clobber=True)

        # Write spectrum specifying wavelength-related header keywords
        # manually
        wp = {"CRVAL1": 5000., "CDELT1": 0.01, "CRPIX1": 1}
        pyasl.write1dFitsSpec("test2.fits", flux, waveParams=wp, clobber=True)


class SanityOftgal_uvw(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def sanity_example(self):
        """
        Checking sanity of gal_uvw example
        """
        from PyAstronomy import pyasl
        
        ra, dec = 017.4294147291650, +61.5473037149441
        pmra, pmdec = 627.89, 77.84
        d = 144
        rv = -321.4
        
        print("U, V, W [kms/s]: % 5.3f, % 5.3f, % 5.3f " % \
            pyasl.gal_uvw(ra, dec, pmra, pmdec, d, rv))
        print("U, V, W (LSR corrected) [kms/s]:  % 5.3f, % 5.3f, % 5.3f" % \
            pyasl.gal_uvw(ra, dec, pmra, pmdec, d, rv, lsr=(8.5, 13.38, 6.49) ) )

    def sanity_galuvw(self):
        """
        Checking sanity of gal_uvw 
        """
        from PyAstronomy import pyasl
        
        ra, dec = 017.4294147291650, +61.5473037149441
        pmra, pmdec = 627.89, 77.84
        d = 144
        rv = -321.4

        u0, v0, w0 = -154, -493, 97
        
        u, v, w = pyasl.gal_uvw(ra, dec, pmra, pmdec, d, rv, lsr=(8.5, 13.38, 6.49))
        self.assertAlmostEqual(u-u0, 0, delta=1.5, msg="Problem with U")
        self.assertAlmostEqual(v-v0, 0, delta=1.5, msg="Problem with V")
        self.assertAlmostEqual(w-w0, 0, delta=1.5, msg="Problem with W")
        
        ra, dec = 200., -50.
        pmra, pmdec = 235., -571.
        d = 37.6
        rv = 112.
        u, v, w = pyasl.gal_uvw(ra, dec, pmra, pmdec, d, rv)
        # from http://www.astro.ucla.edu/~drodrigu/UVWCalc.html#data
        u0, v0, w0 = 104.66,   -86.50,   -78.86
        self.assertAlmostEqual(u-u0, 0, delta=0.1, msg="Problem with U (no 2)")
        self.assertAlmostEqual(v-v0, 0, delta=0.1, msg="Problem with V (no 2)")
        self.assertAlmostEqual(w-w0, 0, delta=0.1, msg="Problem with W (no 2)")
        

class SanityOfProperMotion(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def sanity_pm(self):
        """
        Test sanity of proper motion
        """
        from PyAstronomy import pyasl
        import numpy as np
        # Year 2018 (GJ 1215)
        ra, dec = 259.4319131981014, +11.6678903453170
        pmra, pmdec = -352.897, -445.558
        
        # 2100 by SIMBAD
        ra2100, dec2100 = 259.4237057511850, +11.6577413569938
        # 2050 by SIMBAD
        ra2050, dec2050 = 259.4287102205105, +11.6639297921180

        ra2, dec2 = pyasl.applyProperMotion(ra, dec, pmra, pmdec, dt=2100-2018, fixes=1)
        ra5, dec5 = pyasl.applyProperMotion(ra, dec, pmra, pmdec, dt=2050-2018, fixes=1)
        
        self.assertAlmostEqual(ra2, ra2100, delta=1e-6, msg="Ra2100 does not match")
        self.assertAlmostEqual(ra5, ra2050, delta=1e-6, msg="Ra2050 does not match")
        self.assertAlmostEqual(dec2, dec2100, delta=1e-6, msg="Dec2100 does not match")
        self.assertAlmostEqual(dec5, dec2050, delta=1e-6, msg="Dec2050 does not match")
        
        ra, dec = 0, 0
        pmra, pmdec = 1, 0
        ra2, dec2 = pyasl.applyProperMotion(ra, dec, pmra, pmdec, dt=5, fixes=1)
        self.assertAlmostEqual(ra2, 5/3600/1000, delta=1e-6, msg="Ra 1,0 does not match")
        self.assertAlmostEqual(dec2, 0, delta=1e-6, msg="Dec 1,0 does not match")
        
        np.random.seed(997177)
        for i in range(40):
        
            ra, dec = np.random.random()+360, np.random.random()*178-79.
            pmra, pmdec = np.random.random()*1000-500, np.random.random()*1000-500
            ra2, dec2 = pyasl.applyProperMotion(ra, dec, pmra, pmdec, dt=43, fixes=1)
            ra2, dec2 = pyasl.applyProperMotion(ra2, dec2, pmra, pmdec, dt=-43, fixes=1)
            self.assertAlmostEqual(ra, ra2, delta=1e-6, msg="Forth-Back test failed (ra, no. %d)" % i)
            self.assertAlmostEqual(dec, dec2, delta=1e-6, msg="Forth-Back test failed (dec, no. %d)" % i)

        # Year 2018 (GJ 1215)
        ra, dec = 259.4319131981014, +11.6678903453170
        pmra, pmdec = -352.897, -445.558
        ra2, dec2 = pyasl.applyProperMotion(ra, dec, pmra, pmdec, dt=2100-2018, fixes=1)
        for i in range(1, 6, 1):
            ra3, dec3 = pyasl.applyProperMotion(ra, dec, pmra, pmdec, dt=2100-2018, fixes=i)
            self.assertAlmostEqual(ra2, ra3, delta=1e-6, msg="Fixes test failed (ra, no. %d)" % i)
            self.assertAlmostEqual(dec2, dec3, delta=1e-6, msg="Fixes test failed (dec, no. %d)" % i)

    def sanity_example(self):
        """
        Testing proper motion example
        """
        from PyAstronomy import pyasl
        
        # Year 2018 coordinates and proper motion (mas/yr)
        # of GJ 1215
        ra, dec = 259.4319131981014, +11.6678903453170
        pmra, pmdec = -352.897, -445.558
        
        # Calculate 2050 position
        ra5, dec5 = pyasl.applyProperMotion(ra, dec, pmra, pmdec, dt=2050-2018, fixes=1)
        
        print("Position 2018: %10.6f  % 10.6f" % (ra, dec))
        print("Position 2050: %10.6f  % 10.6f" % (ra5, dec5))


class SanityOfBetaSigma(unittest.TestCase):
 
    def sanity_test1(self):
        """
        Beta Sigma sanity test 1
        """
        import numpy as np
        from PyAstronomy import pyasl
        
        np.random.seed(192837465)
        x = np.arange(1000)
        y = np.random.normal(0, 1, len(x))
 
        na = pyasl.BSArbSamp()
        ne = pyasl.BSEqSamp()
 
        e, se = ne.betaSigma(y, 0, 1, ibs=True)
        a, sa = na.betaSigma(x, y, 0, 1, ibs=True)
 
        print("e, se: ", e, se)
        print("a, sa: ", a, sa)
 
        self.assertAlmostEqual(e, a, delta=1e-10, msg="SanityOfBetaSigma, test1: e and a do not match!")
        self.assertAlmostEqual(se, sa, delta=1e-10,
                               msg="SanityOfBetaSigma, test1: se and sa do not match!")
        self.assertAlmostEqual(e, 1., delta=2 * se,
                               msg="SanityOfBetaSigma, test1: e does not match 1 within 2 std.")
 
    def sanity_test2(self):
        """
        Beta Sigma sanity test 2
        """
        import numpy as np
        from PyAstronomy import pyasl
        
        np.random.seed(192837466)
        x = np.arange(1000)
        y = np.random.normal(0, 1, len(x))
 
        na = pyasl.BSArbSamp()
        ne = pyasl.BSEqSamp()
 
        stds = []
 
        for N in smo.range(10):
            for j in smo.range(1, 10, 1):
 
                e, se = ne.betaSigma(y, N, j, ibs=True)
                a, sa = na.betaSigma(x, y, N, j, ibs=True)
 
                print("N, j, n(betaSample): ", N, j, len(
                    ne.betaSample), len(na.betaSample))
                print("    e, se: ", e, se)
                print("    a, sa: ", a, sa)
 
                stds.append(e)
 
                self.assertAlmostEqual(
                    e, a, delta=1e-7, msg="SanityOfBetaSigma, test2: e and a do not match! e - a = " + str(e - a))
                self.assertAlmostEqual(
                    se, sa, delta=1e-7, msg="SanityOfBetaSigma, test2: se and sa do not match!")
                self.assertAlmostEqual(
                    e, 1., delta=4 * se, msg="SanityOfBetaSigma, test2: e does not match 1 within 4 std.")
 
    def sanity_test_betaSampleIndices(self):
        """
        Beta sigma sample indices for independent sample construction
        """
 
        import numpy as np
        from PyAstronomy import pyasl
        
        x = np.arange(9)
        y = np.random.normal(0, 1, len(x))
 
        ne = pyasl.BSEqSamp()
 
        _ = ne.betaSigma(x, 1, 1, ibs=True)
 
        self.assertEqual(True, np.any(ne.karr[0] == np.array([0,  3,  6])))
        self.assertEqual(True, np.any(ne.karr[1] == np.array([0,  3,  6]) + 1))
        self.assertEqual(True, np.any(ne.karr[2] == np.array([0,  3,  6]) + 2))
 
        _ = ne.betaSigma(x, 0, 3, ibs=True)
 
        self.assertEqual(True, np.any(ne.karr[0] == np.array([0, 1, 2])))
        self.assertEqual(True, np.any(ne.karr[1] == np.array([0, 1, 2]) + 3))
 
        x = np.arange(1000)
        y = np.random.normal(0, 1, len(x))
        for N in smo.range(10):
            for j in smo.range(10):
                _ = ne.betaSigma(x, N, j + 1, ibs=True)
                ka = np.empty(shape=(0,))
                for k in ne.karr:
                    ka = np.concatenate((ka, k))
 
                u = np.unique(ka)
                self.assertEqual(len(ka), len(u))
                print(len(ka), len(x), N, j)
                self.assertTrue(len(ka) >= len(x) - (N + 2) * (j + 1))
 
    def test3(self):
        """
        Beat sigma test3
        """
 
        import numpy as np
        from PyAstronomy import pyasl
        
        np.random.seed(773311)
        x = np.arange(1000)
        y = np.random.normal(0, 1, len(x))
 
        ne = pyasl.BSEqSamp()
 
        std, e = ne.stdUnbiased(y)
        madstd = ne.estimateStdMAD(y)
 
        self.assertAlmostEqual(std, madstd, delta=e,
                               msg="SanityOfBetaSigma, test3: std and madstd do not match")
 
 
    def sanity_example_quickstart(self):
        """
        Beta-Sigma Quickstart example
        """
        from PyAstronomy import pyasl
        import numpy as np
        
        # Set up artificial data set
        x = np.linspace(0., 1., 200)
        y = 10.2 + x + 4 * x**2
        # Add some independent, Gaussian noise
        gstd = 0.1
        y += np.random.normal(0., gstd, len(x))
        
        # Estimate noise using robust estimate
        beq = pyasl.BSEqSamp()
        # Define order of approximation (use larger values such as 2,3, or 4 for
        # faster varying or less well sampled data sets; also 0 is a valid order)
        N = 1
        # Define 'jump parameter' (use larger values such as 2,3, or 4 if correlation
        # between adjacent data point is suspected)
        j = 1
        # Estimate noise assuming equidistant sampling (often a good approximation even
        # if data are not strictly equidistant) and robust estimation (often advantageous
        # in working with real data)
        nstd, nstdstd = beq.betaSigma(y, N, j, returnMAD=True)
        
        print("Input noise std: ", gstd)
        print("Estimated noise std = %5.3f +/- %5.3f" % (nstd, nstdstd))

    def sanity_example_different_N(self):
        """
        BS example: Noise estimates with different orders of approximations (N)
        """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy import pyasl
        
        def g(t):
            """
            Function determining the behavior of the data.
            """
            return 1.3 - 0.003*t + 1.3*np.sin(t/5.) * np.exp(-t/100.)
        
        
        # Number of data points
        nd = 200
        
        ti = np.arange(float(nd))
        gi = g(ti)
        
        mdiff = np.max(np.abs(gi[1:] - gi[0:-1]))
        print("""Maximum absolute difference between consecutive
            values of g(t): """, mdiff)
        
        # Standard deviation of noise
        istd = 0.02
        
        print("Input standard deviation: ", istd)
        print("Number of 'data points': ", nd)
        print()
        
        # Add Gaussian noise to data
        yi = gi + np.random.normal(0.0, istd, nd)
        
        # Create class instance for equidistant sampling
        bseq = pyasl.BSEqSamp()
        
        # Specify jump parameter (j) for construction of beta sample
        j = 1
        
        # Order of approximation to use
        Ns = [0,1,2,3]
        
        # Use to store noise estimates
        smads, dsmads = [], []
        
        # Loop over orders of approximation between 0 and 3
        for N in Ns:
        
            # Get estimates of standard deviation based on robust (MAD-based) estimator
            smad, dsmad = bseq.betaSigma(yi, N, j, returnMAD=True)
            print("Order of approximation (N): ", N)
        
            print("    Size of beta sample: ", len(bseq.betaSample))
            print("    Robust estimate of noise std: %6.3f +/- %6.3f" % (smad, dsmad))
            # Save result
            smads.append(smad)
            dsmads.append(dsmad)
        
#         # Plot g(t) and the synthetic data
#         plt.subplot(2,1,1)
#         plt.title("Data (top) and noise estimates (bottom)")
#         plt.plot(ti, gi, 'b.-', label="$g(t_i)$")
#         plt.errorbar(ti, yi, yerr=np.ones(nd)*istd, fmt='r+', label="$y_i$")
#         plt.legend()
#         plt.subplot(2,1,2)
#         plt.title("N=0 is insufficient")
#         plt.errorbar(Ns, smads, yerr=dsmads, fmt='k+', label="Noise estimates")
#         plt.plot([min(Ns)-0.5, max(Ns)+0.5], [istd]*2, 'k--', label="Input value")
#         plt.legend()
#         plt.xlabel("Order of approximation (N)")
#         plt.ylabel("Noise STD")
#         plt.tight_layout()
#         plt.show()

    def sanity_example_outliers(self):
        """
        BS example: The effect of outliers: A case for robust estimation
        """
        import numpy as np
        from PyAstronomy import pyasl
        import matplotlib.pylab as plt
        
        # Number of data points
        nd = 200
        
        # Input standard deviation
        istd = 0.2
        
        # Create some synthetic data (Gaussian noise) with
        # input standard deviation.
        y = np.random.normal(1.0, istd, nd)
        
        # Introduce some outliers
        # Number of outliers
        no = 15
        for _ in range(no):
            # Choose an index (could occur more than once)
            index = np.random.randint(0, high=nd)
            # Add point from normal distribution with
            # STD of 50
            y[index] = np.random.normal(1.0, 50.0)
        
        # Create class instance for equidistant sampling
        bseq = pyasl.BSEqSamp()
        
        # Specify order of approximation (N) and jump parameter (j) for
        # construction of beta sample
        N = 0
        j = 1
        
        print("Order of approximation: ", N)
        print("Jump parameter: ", j)
        print()
        
        # Get estimates of standard deviation based on MV estimator and ...
        smv, dsmv = bseq.betaSigma(y, N, j, returnMAD=False)
        # ... robust (MAD-based) estimator
        smad, dsmad = bseq.betaSigma(y, N, j, returnMAD=True)
        
        print("Input standard deviation: ", istd)
        print("Number of 'data points': ", nd)
        print("Size of beta sample: ", len(bseq.betaSample))
        print()
        print("Minimum-variance estimate: %6.3f +/- %6.3f" % (smv, dsmv))
        print("Robust estimate: %6.3f +/- %6.3f" % (smad, dsmad))
        
#         plt.subplot(2,1,1)
#         plt.title("Synthetic data")
#         plt.plot(y, 'bp')
#         plt.subplot(2,1,2)
#         plt.title("Histogram of $\\beta$ sample")
#         plt.hist(bseq.betaSample, 30)
#         plt.show()
 
    def sanity_example_eqarbsamp(self):
        """
        Equidistant vs. arbitrary sampling
        """
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy import pyasl
        
        def g(t):
            """
            Function determining the behavior of the data.
            """
            return 1.3 - 10.0*t
        
        
        # Number of data points
        nd = 30
        
        # Creating non-equidistant samping axis (ti)
        te = np.arange(float(nd))
        ti = (te**3) / float(nd**2)
        
        # Get values of g(t)
        gi = g(ti)
        
        # Standard deviation of noise
        istd = 0.3
        
        # Add Gaussian noise to data
        yi = gi + np.random.normal(0.0, istd, nd)
        
        print("Input standard deviation: ", istd)
        print("Number of 'data points': ", nd)
        print()
        
        # Create class instance for equidistant sampling
        bseq = pyasl.BSEqSamp()
        # Create class instance for arbitrary sampling
        bsar = pyasl.BSArbSamp()
        
        # Get estimates assung equidistant and arbitrary sampling
        # using N = 1 and j = 1. From the definition of g(t), N = 1
        # will be sufficient for the case of arbitrary sampling, but
        # not necessarily for (assumed) equidistant sampling.
        smv_es, dsmv_es  = bseq.betaSigma(yi, 1, 1)
        smv_as, dsmv_as  = bsar.betaSigma(ti, yi, 1, 1)
        
        print("Estimates for N=1 and j=1")
        print("    Equidistant sampling: %5.3f +/- %5.3f" % (smv_es, dsmv_es))
        print("    Arbitrary sampling: %5.3f +/- %5.3f" % (smv_as, dsmv_as))
        print()
        
        # Get estimates for N=2 and 3 assuming equidistant sampling
        smv_es2, dsmv_es2  = bseq.betaSigma(yi, 2, 1)
        smv_es3, dsmv_es3  = bseq.betaSigma(yi, 3, 1)
        
        print("Estimates for N=2 and 3 based on equidistant sampling")
        print("    N = 2: %5.3f +/- %5.3f" % (smv_es2, dsmv_es2))
        print("    N = 3: %5.3f +/- %5.3f" % (smv_es3, dsmv_es3))
        
#         plt.subplot(2,1,1)
#         plt.title("Data with true sampling")
#         plt.plot(ti, gi, 'b-')
#         plt.errorbar(ti, yi, yerr=np.ones(nd)*istd, fmt='b+')
#         plt.subplot(2,1,2)
#         plt.title("Same data assuming equidistant sampling")
#         plt.plot(te, gi, 'r-')
#         plt.errorbar(te, yi, yerr=np.ones(nd)*istd, fmt='r+')
#         plt.tight_layout()
#         plt.show()


