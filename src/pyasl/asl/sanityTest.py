from airtovac import airtovac, vactoair
import numpy
import unittest
from PyAstronomy.pyaC import SaniBase
from aitoffLegacy import inverseAitoff, aitoff
from astroTimeLegacy import daycnv, juldate
from keplerOrbit import KeplerEllipse, MarkleyKESolver
import datetime
from dopplerShift import dopplerShift
from fluxConversion import flux2photons, photons2flux
from rotBroad import rotBroad, fastRotBroad
from cardinalPoint import getCardinalPoint

class SanityOfPyasl(unittest.TestCase, SaniBase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_airtovac(self):
    """
      Checking `airtovac` and back...
    """
    wvl_air = numpy.random.random(100) * 10000.0 + 1000.0
    wvl_vac = airtovac(wvl_air)
    wvl_air2 = vactoair(wvl_vac)
    self.assertTrue(self.mrd(wvl_air, wvl_air2) < 1e-6)
  
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
    wvl = np.linspace(5000.,5010.,1000.)
    flux = g.evaluate(wvl)
    # ... and add some noise
    flux += np.random.normal(0., 0.02, len(wvl))
    
    # Assume that this spectrum is one observed
    # or modeled in vacuum conditions. We want
    # to convert it into a spectrum in air, but
    # retain the old wavelength axis.
    airflux, vind = pyasl.specAirVacConvert(wvl, flux, \
                    direction="vactoair")
    
    # Plot the result
    plt.plot(wvl, flux, "b.-")
    plt.plot(wvl, airflux, "r.-")
#    plt.show()
  
  def sanity_aitoff(self):
    """
      Checking Aitoff projection.
    """
    for i in xrange(100):
      l = numpy.random.random() * 360.0
      b = numpy.random.random() * 180.0 - 90.0
      x, y = aitoff(l, b)
      l2, b2 = inverseAitoff(x, y)
      self.assertTrue(self.mrd(l2, l) < 1e-6)
      self.assertTrue(self.mrd(b2, b) < 1e-6)
  
  def sanity_daycnv(self):
    """
      Checking `daycnv`.
    """
    y, m, d, h = tuple(daycnv(2451596.812500))
    self.assertEqual(y, 2000)
    self.assertEqual(m, 2)
    self.assertEqual(d, 22)
    self.assertEqual(h, 7.5)
  
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
    for i in xrange(len(time)):
      pos1 = ke.xyzPos(time[i])
      self.assertEqual( numpy.sum(numpy.abs((pos[i,::]-pos1))), 0.0)
  
  def sanity_MarkleyKESolver_precision(self):
    """
      Checking precision of Markley solver for Kepler's equation.
    """
    ks = MarkleyKESolver()
    p = ks.precisionTest()
    self.assertLess(p, 1e-14)


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
    print "Eccentric anomaly: ", ks.getE(M, e)
  
  def sanity_example2(self):
    """
      Checking Kepler orbit calculation example.
    """
    import numpy
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    
    # Instantiate a Keplerian elliptical orbit with
    # semi-major axis of 1.3 length units,
    # period of 2 time units, eccentricity of 0.5, and
    # longitude of ascending node of 70 degrees.
    ke = pyasl.KeplerEllipse(1.3, 2., e=0.5, Omega=70.)
    
    # Get a time axis
    t = numpy.linspace(0, 6.5, 200)
    
    # Calculate the orbit position at the given points
    # in a Cartesian coordinate system.
    pos = ke.xyzPos(t)
    print "Shape of output array: ", pos.shape
    
    # x, y, and z coordinates for 50th time point
    print "x, y, z for 50th point: ", pos[50, ::]
    
    # Calculate orbit radius as a function of the
    radius = ke.radius(t)
    
    # Plot x and y coordinates of the orbit
    plt.subplot(2,1,1)
    plt.plot(pos[::,0], pos[::,1], 'bp')
    # Plot orbit radius as a function of time
    plt.subplot(2,1,2)
    plt.plot(t, radius, 'bp')
    # plt.show()


class SanityOfBinnin(unittest.TestCase, SaniBase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_example_binningx0dt(self):
    """
      Checking `binningx0dt` example.
    """
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy.pyasl import binningx0dt
    
    # Generate some data
    x = np.arange(999)
    y = np.sin(x/100.)
    y += np.random.normal(0,0.1,len(x))
    
    # Bin using fixed number of bins and start at x0 = -10.
    # Use beginning of bin as starting value.
    r1, dt1 = binningx0dt(x, y, nbins=50, x0=-10, useBinCenter=False)
    # Use fixed bin width. Specify another (wrong) error estimate and
    # use bin center.
    r2, dt2 = binningx0dt(x, y, yerr=np.ones(len(x))*0.2, dt=dt1, \
                          x0=-10, useBinCenter=True, removeNoError=True)
    
    print "dt1, dt2: ", dt1, dt2
    print "Input data points in last bin: ", r2[-1,3]
    
    # Use the reducedBy flag to indicate the binning. In this case, x0
    # will be set to the lowest x value in the data, and the number of
    # bins will be calculated as: int(round(len(x)/float(reduceBy))).
    # Here, we will, thus, obtain 100 bins.
    r3, dt3 = binningx0dt(x, y, \
                          useBinCenter=True, removeNoError=True, reduceBy=10)
    
    print "dt3: ", dt3
    print "Number of bins in theird version: ", len(r3[::,0])
    
    
    # Plot the output
    plt.plot(x,y)
    plt.errorbar(r1[::,0], r1[::,1], yerr=r1[::,2], fmt='kp--')
    plt.errorbar(r2[::,0], r2[::,1], yerr=r2[::,2], fmt='rp--')
    plt.errorbar(r3[::,0], r3[::,1], yerr=r3[::,2], fmt='gp--')
#    plt.show()


class SanityOfPhotonConversion(unittest.TestCase, SaniBase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_fluxConversion(self):
    """
      Checking `flux2photons` and `photons2flux`.
    """
    wvl = numpy.linspace(300., 20000., 100.)
    flux = 1./wvl
    
    photons = flux2photons(wvl, flux)
    flux2 = photons2flux(wvl, photons)
    self.assertLess(max(numpy.abs(flux-flux2)), 1e-12)

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
    print "%g erg/s at %g A correspond to %g photons/s" \
            % (flux, wvl, photons)
    
    # Converting back
    flux2 = pyasl.photons2flux(wvl, photons)
    
    print "%g photons/s at %g A correspond to %g erg/s" \
            % (photons, wvl, flux2)



class SanityOfFolding(unittest.TestCase, SaniBase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_foldAtExample(self):
    """
      Checking `flotAt` example.
    """
    from PyAstronomy.pyasl import foldAt
    import matplotlib.pylab as plt
    import numpy as np
    
    # Generate some data ...
    time = np.random.random(1000) * 100.
    flux = 0.05 * np.sin(time*(2.*np.pi/21.5) + 15)
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
    flux = numpy.ones(len(wvl)) + (wvl/wvl.min())*0.05
    flux -= numpy.exp( -(wvl-6050.)**2/(2.*0.5**2) )*0.05
    
    nflux, wlprime = dopplerShift(wvl, flux, 20., edgeHandling="firstlast")
    nflux2, wlprime = dopplerShift(wvl, nflux, -20., edgeHandling="firstlast")

    self.assertLess(max(numpy.abs(flux-nflux2)), 1e-5)

  def sanity_checkExample(self):
    """
      Doppler shift -- Checking whether example works
    """
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    import numpy as np
    
    # Create a "spectrum" with 0.01 A binning ...
    wvl = np.linspace(6000., 6100., 10000)
    # ... a gradiant in the continuum ...
    flux = np.ones(len(wvl)) + (wvl/wvl.min())*0.05
    # ... and a Gaussian absoption line
    flux -= np.exp( -(wvl-6050.)**2/(2.*0.5**2) )*0.05
    
    # Shift that spectrum redward by 20 km/s using
    # "firstlast" as edge handling method.
    nflux1, wlprime1 = pyasl.dopplerShift(wvl, flux, 20., edgeHandling="firstlast")
    
    # Shift the red-shifted spectrum blueward by 20 km/s, i.e.,
    # back on the initial spectrum.
    nflux2, wlprime = pyasl.dopplerShift(wvl, nflux1, -20., \
                            edgeHandling="fillValue", fillValue=1.0)
    
    # Check the maximum difference in the central part
    indi = np.arange(len(flux)-200) + 100
    print "Maximal difference (without outer 100 bins): ", \
                    max(np.abs(flux[indi]-nflux2[indi]))
    
    # Plot the outcome
    plt.title("Initial (blue), shifted (red), and back-shifted (green) spectrum")
    plt.plot(wvl, flux, 'b.-')
    plt.plot(wvl, nflux1, 'r.-')
    plt.plot(wvl, nflux2, 'g.-')
#    plt.show()


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
    print "Initional EW [A]: ", 4. - sci.trapz(flux, wvl)
    print "After broadening without LD: ", 4. - sci.trapz(rflux, wvl)
    print "After broadening with LD: ", 4. - sci.trapz(lflux, wvl)
    
    # Plot the results
    plt.title("Rotational broadening")
    plt.xlabel("Wavelength [A]")
    plt.ylabel("Normalized flux")
    plt.plot(wvl, flux, 'b-')
    plt.plot(wvl, rflux, 'r-')
    plt.plot(wvl, lflux, 'g-')
#    plt.show()

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
    
    for broad in np.linspace(0.1, 100., 1.):
      for eps in np.linspace(0., 1., 0.2):
        nflux = rotBroad(wvl, flux, epsilon=eps, vsini=broad)
        ew = 4. - sci.trapz(nflux, wvl)
        self.assertAlmostEqual(ew, refEW, 5, "The EW in rotbroad changes for: vsini = " + \
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
    self.assertAlmostEqual(numpy.max(numpy.abs( bfast[indi] - bslow[indi] )), 0.0, delta=1e-5)
    
    flux[150] = 1.0
    A = -0.05
    s = 0.1
    mu = 5004.1635788
    flux += A/numpy.sqrt(2.*numpy.pi*s**2) * \
      numpy.exp(-(wvl-mu)**2/(2.*s**2))
    
    bfast = fastRotBroad(wvl, flux, 0.81, 11.37)
    bslow = rotBroad(wvl, flux, 0.81, 11.37)
   
    indi = numpy.where(numpy.logical_and(wvl > 5000., wvl < 5010.))[0]
    self.assertAlmostEqual(numpy.max(numpy.abs( bfast[indi] - bslow[indi] )), 0.0, delta=1e-5)

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
    
    print "Earth's velocity at JD: ", jd
    print "Heliocentric velocity [km/s]: ", heli
    print "Barycentric velocity [km/s] : ", bary
    
    # Coordinates of Sirius
    ra  = 101.28715535
    dec = -16.71611587
    
    vh, vb = pyasl.baryCorr(jd, ra, dec, deq=2000.0)
    print "Barycentric velocity of Earth toward Sirius: ", vb

  def sanity_baryCorrAltair(self):
    """
      Checking barycentric correction for Altair (IDL example)
    """
    from PyAstronomy import pyasl
    vh, vb = pyasl.baryCorr(2449398.5000000000, 297.69488, 8.8676389, deq=2000)
    self.assertAlmostEqual(vb, 10.583091, delta=1e-5)
    
    

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
    nflux = pyasl.rotBroad(wvl, template, 0.5, 23.45, edgeHandling="firstlast")
    
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
    nflux = pyasl.rotBroad(wvl, template, 0.5, 23.45, edgeHandling="firstlast")
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
      chi.append( ((nflux[mind] - m)**2/0.005**2).sum() / len(mind) )
    
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
    flux = pyasl.planck(T=5778., lam=wvl*1e-10)
    
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
    y = np.random.normal(0.0, 1.0/SNR, N) + 1.0
    y += np.sin(x/500.0*2.*np.pi)*0.1
    
    # Estimate the signal to noise ratio. Check whether the
    # estimate fits the input...
    # Use a chunk length of 20 data points, a polynomial of degree
    # one, and produce a "control plot".
    snrEsti = pyasl.estimateSNR(x, y, 20, deg=1, controlPlot=False)
    print "Estimate of the SNR: ", snrEsti["SNR-Estimate"]
    
    # Use a chunks with a length of 27, a polynomial of degree
    # two, and produce a "control plot".
    snrEsti = pyasl.estimateSNR(x, y, 27, deg=2, controlPlot=False, xlenMode="excerpt")
    print "Estimate of the SNR: ", snrEsti["SNR-Estimate"]

    
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
    x = np.array(map(lambda x: float(x),
            "-0.25 0.68 0.94 1.15 1.20 1.26 1.26 1.34 1.38 1.43 1.49 1.49 \
              1.55 1.56 1.58 1.65 1.69 1.70 1.76 1.77 1.81 1.91 1.94 1.96 \
              1.99 2.06 2.09 2.10 2.14 2.15 2.23 2.24 2.26 2.35 2.37 2.40 \
              2.47 2.54 2.62 2.64 2.90 2.92 2.92 2.93 3.21 3.26 3.30 3.59 \
              3.68 4.30 4.64 5.34 5.42 6.01".split()))
    
    # Apply the generalized ESD
    r = pyasl.generalizedESD(x, 10, 0.05, fullOutput=True)
    
    print "Number of outliers: ", r[0]
    print "Indices of outliers: ", r[1]
    print "        R      Lambda"
    for i in range(len(r[2])):
      print "%2d  %8.5f  %8.5f" % ((i+1), r[2][i], r[3][i])
    
    # Plot the "data"
    plt.plot(x, 'b.')
    # and mark the outliers.
    for i in range(r[0]):
      plt.plot(r[1][i], x[r[1][i]], 'rp')
#    plt.show()

  def sanity_pointDistGESDExample(self):
    """
      Sanity of pointDistGESD example
    """
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import pyasl
    
    # Get some data
    x = np.random.normal(0.,0.1,50)
    
    # Introduce outliers
    x[27] = 1.0
    x[43] = -0.66
    
    # Run distance based outlier detection
    r = pyasl.pointDistGESD(x, 5)
    
    print "Number of outliers detected: ", r[0]
    print "Indices of these outliers: ", r[1]
    
    plt.plot(x, 'b.')
    for i in range(len(r[1])):
      plt.plot(r[1][i], x[r[1][i]], 'rp')
#    plt.show()


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
    print "Absolute bolometric magnitude of the Sun: ", absMagSun
    print "  Absolute luminosity [erg/s]: ", pyasl.absMagToPower(absMagSun)


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
    print "JD = " + str(jd)
    pos = pyasl.sunpos(jd, full_output=True)
    print "Coordinates of the Sun (ra, dec): %g, %g" % (pos[1], pos[2])
    print "Solar elongation = %g and obliquity = %g" % (pos[3], pos[4])
    
    # Get the Sun's RA and DEC values for a period of time.
    startjd = datetime.datetime(2013, 4, 16)
    endjd = datetime.datetime(2013, 6, 16)
    # Convert into Julian dates
    startjd = pyasl.jdcnv(startjd)
    endjd = pyasl.jdcnv(endjd)
    print
    pos = pyasl.sunpos(startjd, end_jd=endjd, jd_steps=10, plot=False, full_output=True)
    
    for i in range(len(pos[0])):
      print "At JD = %g: ra = %g, dec = %g" % (pos[0][i], pos[1][i], pos[2][i])


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
    print "Nutation for the date."
    res = pyasl.nutate(jd)
    print "JD = " + str(jd) + ", Longitude = " + str(res[0]) + \
          ", Obliquity = " + str(res[1])
    
    # Get nutation for an array of JDs.
    startjd = datetime.datetime(2013, 4, 16)
    endjd = datetime.datetime(2013, 6, 16)
    startjd = pyasl.jdcnv(startjd)
    endjd = pyasl.jdcnv(endjd)
    jds = np.arange(startjd, endjd, .5)
    print
    print "Plot the results"
    res = pyasl.nutate(jds, plot=False)
    
    print "Longitude: ", res[0]
    print "Obliquity: ", res[1]
  
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
    print "Get change in RA and DEC due to Earth's nutation for JD = " \
          + str(jd)
    print pyasl.co_nutate(jd, ra, dec)
    
    print
    print "Get change for several RAs and DECs for the same JD"
    ra = np.arange(0.,160.,20.)
    dec = np.arange(-80.,80.,20.)
    res = pyasl.co_nutate(np.repeat(jd, ra.size), ra, dec)
    print res[0], res[1]
    
    print
    print "Get change for several RAs and DECs for different JDs"
    jds = np.arange(jd,jd+ra.size,1)
    res = pyasl.co_nutate(jds, ra, dec)
    print "JD             delta(RA)   delta(DEC)"
    for i in range(ra.size):
      print "%12.5f   %8.5f   %8.5f" % (jds[i], res[0][i], res[1][i])


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
    print "Get change in RA and DEC due to annual aberration"
    print "  for JD = " + str(jd) + ":", \
          np.ravel(pyasl.co_aberration(jd, ra, dec))
    
    print
    print "Get change for several RAs and DECs for the same JD"
    ra = np.arange(10.,50.,10.)
    dec = np.arange(30.,70.,10.)
    res = pyasl.co_aberration(np.repeat(jd, ra.size), ra, dec)
    print res[0], res[1]
    
    print
    print "Get change for several RAs and DECs for different JDs"
    jds = np.arange(jd,jd+ra.size,1)
    res = pyasl.co_aberration(jds, ra, dec)
    print "JD             delta(RA)   delta(DEC)"
    for i in range(ra.size):
      print "%12.5f   %8.5f   %8.5f" %(jds[i], res[0][i], res[1][i])


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
    print
    print "Get apparent (real) altitude of a star with observed altitude of " + \
          str(alt) + " degrees"
    print "  ->  Apparent altitude = ", alt - pyasl.co_refract_forward(alt)
    
    print
    print "You are not observing from sea level, but from an altitude of 5000 meter."
    print ("Apparent altitude = %9.5f, estimated pressure [mbar] = %9.5f, " + \
          "estimated temperature [K] = %9.5f") % \
          pyasl.co_refract(alt, observer_alt=5000, convert_to_observed=False)
    
    print
    print "Convert apparent (real) altitude into observed altitude."
    print "Apparent altitude = " + str(alt) + " degrees",
    print " -> Observed altitude = " + str(pyasl.co_refract(alt, full_output=False,\
                                           convert_to_observed=True)[0])
    
    print
    print "The same object observed from different observer altitudes"
    apparentAltitudes = np.repeat(30.0, 10)
    obsalts = np.linspace(0.,5000.,len(apparentAltitudes))
    r = pyasl.co_refract(apparentAltitudes, observer_alt=obsalts, convert_to_observed=True)
    for i in range(len(r[0])):
      print "Observed altitude [deg] = %g, pressure [mbar] = %g, temperature [K] = %g" \
            % (r[0][i], r[1][i], r[2][i])


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
    print "Get altitude and azimuth of object in DEGREES"
    print pyasl.hadec2altaz(ha, dec, lat)
    
    # List of coordinates
    ha = np.arange(0.,20.,5.)
    dec = np.arange(30.,50.,5.)
    lat = np.zeros(dec.size)+53.48
    print
    print "Get altitude and azimuth for a list of objects from same observer latitude"
    altaz = pyasl.hadec2altaz(ha, dec, lat)
    print "alt: ", altaz[0]
    print "az: ", altaz[1]
  
  def sanity_eq2horExample(self):
    """
      Sanity of example for `eq2hor`
    """
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    import datetime
    import numpy as np
    
    # Convert calender date to JD
    # use the datetime package
    jd = datetime.datetime(2013, 4, 16)
    jd = pyasl.jdcnv(jd)
    # Specifiy RA and DEC
    ra = 10.
    dec = 30.
    print
    print "Get horizontal coordinates (alt, az, ha) from JD, RA,"
    print "  and DEC for the Hamburger Sternwarte"
    print pyasl.eq2hor(jd, ra, dec, observatory="HS")
    
    print
    print "From a list of Julian dates ..."
    jds = np.arange(jd,jd+1,.2)
    ras = np.zeros(jds.size) + ra
    decs = np.zeros(jds.size) + dec
    alt, az, ha = pyasl.eq2hor(jds, ras, decs, lon=-70.4042, lat=-24.6272, alt=2635.)
    
    for i in range(alt.size):
      print "JD = %g : alt = % g,  az = % g,  ha = % g" % (jds[i], alt[i], az[i], ha[i])
    
    
    print
    print "For one object and different times at the VLT..."
    jds = np.arange(jd-.25,jd+.25,.01)
    ras = np.zeros(jds.size) + 130.
    decs = np.zeros(jds.size) - 30.
    res = pyasl.eq2hor(jds, ras, decs, lon=-70.4042, lat=-24.6272, alt=2635.)
    
    plt.plot(jds, res[0])
    plt.xlabel("Julian date")
    plt.ylabel("Altitude [deg]")
#    plt.show()


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
    
    print
    print "Data for Kitt Peak National Observatory"
    print pyasl.observatory("kpno")
    print "(longitude and latitude in degrees, altitude in meters, and"
    print "time zone in hours West of Greenwich"


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
    
    print "Width of H LyA line at 1215.67 A = %e cm" % \
          pyasl.convertDampingConstant(gLya, 1215.67)


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
      print "Azimuth: {0:6.2f} deg, Cardinal point: {1:1s}".format(azimuth, cp)


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
      print "Altitude = {0:6.2f}, Twilight is called: ".format(alt), \
            pyasl.twilightName(alt)


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
    
    print "Angular distance between the poles (deg):"
    print pyasl.getAngDist(98.0, -90.0, 100., +90.0)
    
    print "Angular distance between Vega and Altair (deg)"
    print pyasl.getAngDist(279.23473479, +38.78368896,297.69582730, +08.86832120)


class SanityOfTransit:
  
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
    reJ = pc.REarth/pc.RJ
    print "Earth radius in Jovian units: ", reJ
    
    # Estimate the duration of Earth's transit
    td = pyasl.transitDuration(1.0, reJ, 1.0, 90.0, 365.0)
    print "The transit of Earth lasts about: %5.3f days" % td
  
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
    duration = 2.2/24.0
    
    # Check whether the time is in-transit
    print "Time is within transit? ",
    if not pyasl.isInTransit(time, T0, period, duration/2.0):
      print "No"
    else:
      print "Yes"
  
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
    duration = 2.2/24.0
    
    # Check whether the time is in-transit
    print "Indices if time points within transit: ",
    print pyasl.isInTransit(times, T0, period, duration/2.0)
    
    print
    print "For each time point, a flag indicating whether it"
    print "is in- or off-transit:"
    print pyasl.isInTransit(times, T0, period, duration/2.0, boolOutput=True)


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
    
    for za in range(0,90,10):
      print "Zenith angle: %2d, airmass = %7.2f" % \
        (za, pyasl.airmassSpherical(za, obsAlt))
  
  def sanity_airmassPPExample(self):
    """
      Example for plane-parallel airmass
    """
    from PyAstronomy import pyasl

    print "Airmass for plane-parallel atmosphere"
    for za in range(0,70,10):
      print "Zenith angle: %2d deg, airmass = %7.2f" % \
        (za, pyasl.airmassPP(za))
  
  def sanity_compAirmasses(self):
    """
      Compare different airmass calculations
    """
    from PyAstronomy import pyasl
    for za in range(0,70,10):
      ampp = pyasl.airmassPP(za)
      amsp = pyasl.airmassSpherical(za, 0.0)
      self.assertAlmostEqual(ampp/amsp, 1.0, delta=0.01)
