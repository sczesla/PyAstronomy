.. _GLS:

The Generalized Lomb-Scargle Periodogram (GLS)
===============================================

.. currentmodule:: PyAstronomy.pyTiming.pyPeriod

The GLS class provides an implementation of the Generalized Lomb-Scargle periodogram
as described by Zechmeister & Kuerster 2009 (A&A 496, 577) and earlier by Ferraz-Mello 1981 (AJ 86, 691).
Compared to the periodogram presented, e.g., by Lomb 1976, the GLS takes
into account measurement errors and a constant term in the fit of the trigonometric function.

The periodogram calculation can be controlled by a number of keyword argument.
Amongst others, the boundaries of the frequency range
and its sampling can be adjusted or specified explicitly, various conventions for the
periodogram normalization can be adopted, and the classical Lomb-Scargle periodogram
can be obtained.

In the following, we provide a number of examples, presenting some of the functionality.  

Example: Simple periodogram with default options
----------------------------------------------------------

This example shows a basic application of the GLS class, using
only default options and an artificial input data set. No errors
are specified.

::
    
    from __future__ import print_function, division
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
    flux = A * np.sin(2. * np.pi*time*f)
    # Adding the noise
    flux += np.random.normal(0, sig, time.size)
    
    # Compute the GLS periodogram with default options.
    clp = pyPeriod.Gls((time, flux))
    
    # Print helpful information to screen
    clp.info()
    
    # and plot power vs. frequency.
    plt.xlabel("Frequency")
    plt.ylabel("Power")
    plt.plot(clp.freq, clp.power, 'b.-')
    plt.show()


Example: Adding an error column and FAP levels
----------------------------------------------------------

In this example, an error column will be added. Additionally, we will
calculate and show a number of False-Alarm-Probability levels to help
distinguish between significant and spurious peaks.

::
    
    from __future__ import print_function, division
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
    flux = A * np.sin(2. * np.pi*time*f)
    # Adding the noise
    flux += np.random.normal(0, sig, time.size)
    # Adding an error column
    err = np.ones(N)*sig
    
    # Compute the GLS periodogram with default options.
    # Choose Zechmeister-Kuerster normalization explicitly
    clp = pyPeriod.Gls((time, flux, err), norm="ZK")
    
    # Print helpful information to screen
    clp.info()
    
    # Define FAP levels of 10%, 5%, and 1%
    fapLevels = np.array([0.1, 0.05, 0.01])
    # Obtain the associated power thresholds
    plevels = clp.powerLevel(fapLevels)
    
    # and plot power vs. frequency.
    plt.xlabel("Frequency")
    plt.ylabel("Power")
    plt.plot(clp.freq, clp.power, 'b.-')
    # Add the FAP levels to the plot
    for i in range(len(fapLevels)):
        plt.plot([min(clp.freq), max(clp.freq)], [plevels[i]]*2, '--',
                 label="FAP = %4.1f%%" % (fapLevels[i]*100))
    plt.legend()
    
    plt.show()




Example: Getting highest-power (best-fit) sine curve
----------------------------------------------------------

The best-fit sine curve is extracted from the periodogram class. 

::
    
    from __future__ import print_function, division
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
    flux = A * np.sin(2. * np.pi*time*f)
    # Adding the noise
    flux += np.random.normal(0, sig, time.size)
    # Adding an error column
    err = np.ones(N)*sig
    
    # Compute the GLS periodogram with default options.
    # Choose Zechmeister-Kuerster normalization explicitly
    clp = pyPeriod.Gls((time, flux, err), norm="ZK")
    
    # Get index associated with highest power
    ifmax = np.argmax(clp.power)
    # and highest power and associated frequency
    pmax = clp.power[ifmax]
    fmax = clp.freq[ifmax]
    # Convert frequency into period
    hpp = 1./fmax
    print("Highest-power period: ", hpp)
    
    # Calculate sine wave associated with 'best-fit' frequency
    bestSine = clp.sinmod(time)
    
    plt.subplot(2, 1, 1)
    plt.title("Data and sine asscoiated with highest-power frequency")
    plt.plot(time, flux, 'b.')
    plt.plot(time, bestSine, 'r--')
    plt.subplot(2, 1, 2)
    plt.title("Folded data")
    plt.plot(time/hpp-np.floor(time/hpp), flux, 'b.')
    plt.show()



API documentation
-----------------------

.. autoclass:: Gls
   :members:
   :inherited-members:
   
