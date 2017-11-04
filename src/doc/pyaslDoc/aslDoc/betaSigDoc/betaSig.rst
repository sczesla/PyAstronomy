Beta Sigma Noise Estimation
===================================

The basic assumptions

    - Noise contribution in individual measurements/bins is independent
    - The noise distribution is Gaussian with standard deviation $\\sigma_0$


The effect of outliers: A case for robust estimation
------------------------------------------------------------

Outliers are a common nuisance in working with real data. They
can strongly affect the :math:`\beta\sigma` noise estimate even in otherwise very benign data sets.
**However, the effect of outliers on the noise estimate is much less significant, if robust estimation is used.**
In many cases, the gain in precision from using a more efficient estimator is insignificant, because the
main uncertainties arise from deviations from the assumption of Gaussian noise. 

::
    
    from __future__ import print_function, division
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
    
    plt.subplot(2,1,1)
    plt.title("Synthetic data")
    plt.plot(y, 'bp')
    plt.subplot(2,1,2)
    plt.title("Histogram of $\\beta$ sample")
    plt.hist(bseq.betaSample, 30)
    plt.show()
