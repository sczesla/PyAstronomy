from __future__ import print_function, division
import unittest


class TestPyPeriodExSanity(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testsanity_FourierSpec(self):
        import numpy
        import matplotlib.pylab as plt

        # Import pyTiming
        from PyAstronomy.pyTiming import pyPeriod

        # Create some evenly sampled artificial data (Poisson noise)
        time = numpy.arange(1000.0) / 10.0
        flux = numpy.random.poisson(1, len(time))

        # Build the TimeSeries instance
        lc = pyPeriod.TimeSeries(time, flux)

        # Compute the Leahy-normalized Fourier transform,
        # plot the time series, and check that the mean
        # power level is 2 as expected.
        fft = pyPeriod.Fourier(lc)
        fig, ax = fft.plot()
        print("Mean power level:", numpy.mean(fft.power))

    #     plt.show()

    def testsanity_ErrorWeightLS(self):
        import numpy
        import matplotlib.pylab as plt
        from PyAstronomy.pyTiming import pyPeriod

        # Create unevenly saplted data with frequency=0.1,
        # measurement error and Gaussian noise
        time = numpy.arange(1000.0) + numpy.random.normal(0.0, 0.1, 1000)
        flux = 0.15 * numpy.sin(2.0 * numpy.pi * time / 10.0)
        # Add some noise
        flux += numpy.random.normal(0, 1, time.size) * 0.5
        error = numpy.ones(time.size) * 0.5

        # Plot the light curve in top panel
        plt.subplot(3, 1, 1)
        plt.errorbar(time, flux, yerr=error)

        # Build the TimeSeries instance
        lc = pyPeriod.TimeSeries(time, flux, error)

        # Compute and plot fast Lomb-Scargle periodogram,
        # which does not take errors into account.
        ls = pyPeriod.LombScargle(lc, ofac=1, hifac=1)
        # Plot the Lomb-Scargle periodogram in middle panel
        plt.subplot(3, 1, 2)
        plt.plot(ls.freq, ls.power, "r-")

        # Compute the full error-weighted Lomb-Periodogram
        # in 'Cumming' normalization and calculate the
        # significance of the maximum peak.
        clp = pyPeriod.Gls(lc, ofac=10, hifac=1, norm="Cumming")
        maxPower = numpy.max(clp.power)
        print("GLS maximum power: ", maxPower)
        print("GLS statistics of maximum power peak: ", clp.stats(maxPower))

        # Plot the generalized Lomb-Scargle periodogram in
        # bottom panel.
        plt.subplot(3, 1, 3)
        plt.plot(clp.freq, clp.power)
        # Show the results

    #     plt.show()

    def testsanity_GlsCodeExample(self):
        """
        Gls in-code example
        """
        import numpy as np
        from PyAstronomy.pyTiming.pyPeriod import Gls

        time = np.random.uniform(54000.0, 56000.0, 1000)
        flux = 0.15 * np.sin(2.0 * np.pi * time / 10.0)

        error = 0.5 * np.ones(time.size)
        flux += np.random.normal(0, error)

        # Compute the full error-weighted Lomb-Periodogram
        # in 'ZK' normalization and calculate the significance
        # of the maximum peak.
        gls = Gls((time, flux, error), verbose=True)

        maxPower = gls.pmax
        print("GLS maximum power: ", maxPower)
        print("GLS statistics of maximum power peak: ", gls.stats(maxPower))


#    gls.plot(block=True)


class TestPyPeriodGLSSanity(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testsanity_example1(self):
        """
        Checking: Simple periodogram with default options
        """
        #     from __future__ import print_function, division
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy.pyTiming import pyPeriod

        # Create some evenly sampled data including a periodic term.
        # Number of data points, input frequency, amplitude of input signal,
        # and standard deviation of noise.
        N = 1000
        f = 0.1
        A = 0.15
        sig = 0.2

        time = np.arange(float(N))
        flux = A * np.sin(2.0 * np.pi * time * f)
        # Adding the noise
        flux += np.random.normal(0, sig, time.size)

        # Compute the GLS periodogram with default options.
        clp = pyPeriod.Gls((time, flux))

        # Print helpful information to screen
        clp.info()

        # and plot power vs. frequency.

    #     plt.xlabel("Frequency")
    #     plt.ylabel("Power")
    #     plt.plot(clp.freq, clp.power, 'b.-')
    #     plt.show()

    def testsanity_example2(self):
        """
        Checking: Adding an error column and FAP levels
        """
        #     from __future__ import print_function, division
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy.pyTiming import pyPeriod

        # Create some evenly sampled data including a periodic term.
        # Number of data points, input frequency, amplitude of input signal,
        # and standard deviation of noise.
        N = 1000
        f = 0.1
        A = 0.05
        sig = 0.2

        time = np.arange(float(N))
        flux = A * np.sin(2.0 * np.pi * time * f)
        # Adding the noise
        flux += np.random.normal(0, sig, time.size)
        # Adding an error column
        err = np.ones(N) * sig

        # Compute the GLS periodogram with default options.
        # Choose Zechmeister-Kuerster normalization explicitly
        clp = pyPeriod.Gls((time, flux, err), norm="ZK")

        # Print helpful information to screen
        clp.info()

        # Define FAP levels of 10%, 5%, and 1%
        fapLevels = np.array([0.1, 0.05, 0.01])
        # Obtain the associated power thresholds
        plevels = clp.powerLevel(fapLevels)

    #     # and plot power vs. frequency.
    #     plt.xlabel("Frequency")
    #     plt.ylabel("Power")
    #     plt.plot(clp.freq, clp.power, 'b.-')
    #     # Add the FAP levels to the plot
    #     for i in range(len(fapLevels)):
    #         plt.plot([min(clp.freq), max(clp.freq)], [plevels[i]]*2, '--',
    #         label="FAP = %4.2f%%" % fapLevels[i])
    #     plt.legend()
    #
    #     plt.show()

    def testsanity_example3(self):
        """
        Checking: Getting highest-power (best-fit) sine curve
        """
        #     from __future__ import print_function, division
        import numpy as np
        import matplotlib.pylab as plt
        from PyAstronomy.pyTiming import pyPeriod

        # Create some evenly sampled data including a periodic term.
        # Number of data points, input frequency, amplitude of input signal,
        # and standard deviation of noise.
        N = 1000
        f = 0.171
        A = 0.5
        sig = 0.2

        time = np.arange(float(N))
        flux = A * np.sin(2.0 * np.pi * time * f)
        # Adding the noise
        flux += np.random.normal(0, sig, time.size)
        # Adding an error column
        err = np.ones(N) * sig

        # Compute the GLS periodogram with default options.
        # Choose Zechmeister-Kuerster normalization explicitly
        clp = pyPeriod.Gls((time, flux, err), norm="ZK")

        # Get index associated with highest power
        ifmax = np.argmax(clp.power)
        # and highest power and associated frequency
        pmax = clp.power[ifmax]
        fmax = clp.freq[ifmax]
        # Convert frequency into period
        hpp = 1.0 / fmax
        print("Highest-power period: ", hpp)

        # Calculate sine wave associated with 'best-fit' frequency
        bestSine = clp.sinmod(time)


#     plt.subplot(2,1,1)
#     plt.title("Data and sine asscoiated with highest-power frequency")
#     plt.plot(time, flux, 'b.')
#     plt.plot(time, bestSine, 'r--')
#     plt.subplot(2,1,2)
#     plt.title("Folded data")
#     plt.plot(time/hpp-np.floor(time/hpp), flux, 'b.')
#     plt.show()
