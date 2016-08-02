Exponentially correlated Gaussian random numbers
==================================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: expCorrRN

Example: Generating correlated numbers
---------------------------------------------

::
  
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
    
    plt.subplot(2,1,1)
    plt.plot(range(200), c1, 'bp-')
    plt.subplot(2,1,2)
    plt.plot(range(200), c2, 'bp-')
    plt.plot(range(200), g, 'g.')
    plt.show()


Example: Exploring the autocorrelation function
---------------------------------------------------

::
  
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import pyasl
    
    # Generate n exponentially correlated Gaussian
    # random numbers with a decay time, tau
    n = 500
    tau = 5.
    c1 = pyasl.expCorrRN(n, tau)
    
    # Obtain autocorrelation function
    ac = np.correlate(c1, c1, mode="full")[n-1:]
    
    # Plot correlated random numbers and autocorrelation
    # function along with exponential model.
    x = np.arange(float(n))
    plt.subplot(2,1,1)
    plt.plot(x, c1, 'bp-')
    plt.subplot(2,1,2)
    plt.plot(x, ac, 'b.')
    plt.plot(x, np.exp(-x/tau)*ac.max(), 'r--')
    plt.show()