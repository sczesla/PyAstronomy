from __future__ import print_function, division
import unittest

class PyPeriodExSanity(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_FourierSpec(self):
    import numpy
    import matplotlib.pylab as plt
    # Import pyTiming
    from PyAstronomy.pyTiming import pyPeriod
    
    # Create some evenly sampled artificial data (Poisson noise)
    time = numpy.arange(1000.)/10.
    flux = numpy.random.poisson(1, len(time))
    
    # Build the TimeSeries instance
    lc = pyPeriod.TimeSeries(time, flux)
    
    # Compute the Leahy-normalized Fourier transform,
    # plot the time series, and check that the mean
    # power level is 2 as expected.
    fft = pyPeriod.Fourier(lc)
    fig, ax = fft.plot()
    print('Mean power level:',numpy.mean(fft.power))
#     plt.show()
  
  def sanity_ErrorWeightLS(self):
    import numpy
    import matplotlib.pylab as plt
    from PyAstronomy.pyTiming import pyPeriod
    
    # Create unevenly saplted data with frequency=0.1,
    # measurement error and Gaussian noise
    time = numpy.arange(1000.) + numpy.random.normal(0., 0.1, 1000)
    flux = 0.15 * numpy.sin(2. * numpy.pi * time / 10.)
    # Add some noise
    flux += numpy.random.normal(0, 1, time.size) * 0.5
    error = numpy.ones(time.size) * 0.5
    
    # Plot the light curve in top panel
    plt.subplot(3,1,1)
    plt.errorbar(time, flux, yerr=error)
    
    # Build the TimeSeries instance
    lc = pyPeriod.TimeSeries(time, flux, error)
    
    # Compute and plot fast Lomb-Scargle periodogram,
    # which does not take errors into account.
    ls = pyPeriod.LombScargle(lc, ofac=1, hifac=1)
    # Plot the Lomb-Scargle periodogram in middle panel
    plt.subplot(3,1,2)
    plt.plot(ls.freq, ls.power, 'r-')
    
    # Compute the full error-weighted Lomb-Periodogram
    # in 'Cumming' normalization and calculate the
    # significance of the maximum peak.
    clp = pyPeriod.Gls(lc, ofac=10, hifac=1, norm="Cumming")
    maxPower = numpy.max(clp.power)
    print("GLS maximum power: ", maxPower)
    print("GLS statistics of maximum power peak: ", clp.stats(maxPower))
    
    # Plot the generalized Lomb-Scargle periodogram in
    # bottom panel.
    plt.subplot(3,1,3)
    plt.plot(clp.freq, clp.power)
    # Show the results
#     plt.show()

