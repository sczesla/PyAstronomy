Examples of Beta Sigma Noise Estimation
=====================================================

In the following, a number of examples demonstrate aspects of the
:math:`\beta\sigma` procedure. Please see `the corresponding paper <http://cdsads.u-strasbg.fr/abs/2018A%26A...609A..39C>`_ 
for further discussion.

    - :ref:`bsOrders`
    - :ref:`bsRobust`
    - :ref:`bsEquisamp`

.. currentmodule:: PyAstronomy.pyasl



.. _bsOrders:

Noise estimates with different orders of approximations (N)
-----------------------------------------------------------------

The required order of approximation (N) depends on the level of variation in
the data (or, equivalently, the sampling cadence in relation to the variability
time scale of the signal). Comparing estimates obtained from various orders of
approximation provides a useful cross-check for the plausibility of the
result. In the example below, it is clearly seen that the zeroth order yields a
too high estimate for the amplitude of noise because the signal varies too quickly.

::
    
    from __future__ import print_function, division
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
    Ns = [0, 1, 2, 3]
    
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
    
    # Plot g(t) and the synthetic data
    plt.subplot(2, 1, 1)
    plt.title("Data (top) and noise estimates (bottom)")
    plt.plot(ti, gi, 'b.-', label="$g(t_i)$")
    plt.errorbar(ti, yi, yerr=np.ones(nd)*istd, fmt='r+', label="$y_i$")
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.title("N=0 is insufficient")
    plt.errorbar(Ns, smads, yerr=dsmads, fmt='k+', label="Noise estimates")
    plt.plot([min(Ns)-0.5, max(Ns)+0.5], [istd]*2, 'k--', label="Input value")
    plt.legend()
    plt.xlabel("Order of approximation (N)")
    plt.ylabel("Noise STD")
    plt.tight_layout()
    plt.show()

.. _bsRobust:

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
    
    plt.subplot(2, 1, 1)
    plt.title("Synthetic data")
    plt.plot(y, 'bp')
    plt.subplot(2, 1, 2)
    plt.title("Histogram of $\\beta$ sample")
    plt.hist(bseq.betaSample, 30)
    plt.show()

  
.. _bsEquisamp:
  
Equidistant vs. arbitrary sampling
--------------------------------------------

The classes :py:func:`BSEqSamp` and :py:func:`BSArbSamp` treat the cases of equidistant and arbitrary
sampling of the signal explicitly. Equidistant sampling is technically more simple to treat and is
often a good approximation *even if* the actual sampling is not equidistant. This is true when there
is a "not too complicated" transformation relating the actual and an equidistant sampling axes.
The following example demonstrates a case, where equidistant sampling can be assumed but a higher
order of approximation is required.
    
::   
    
    from __future__ import print_function, division
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
    smv_es, dsmv_es = bseq.betaSigma(yi, 1, 1)
    smv_as, dsmv_as = bsar.betaSigma(ti, yi, 1, 1)
    
    print("Estimates for N=1 and j=1")
    print("    Equidistant sampling: %5.3f +/- %5.3f" % (smv_es, dsmv_es))
    print("    Arbitrary sampling: %5.3f +/- %5.3f" % (smv_as, dsmv_as))
    print()
    
    # Get estimates for N=2 and 3 assuming equidistant sampling
    smv_es2, dsmv_es2 = bseq.betaSigma(yi, 2, 1)
    smv_es3, dsmv_es3 = bseq.betaSigma(yi, 3, 1)
    
    print("Estimates for N=2 and 3 based on equidistant sampling")
    print("    N = 2: %5.3f +/- %5.3f" % (smv_es2, dsmv_es2))
    print("    N = 3: %5.3f +/- %5.3f" % (smv_es3, dsmv_es3))
    
    plt.subplot(2, 1, 1)
    plt.title("Data with true sampling")
    plt.plot(ti, gi, 'b-')
    plt.errorbar(ti, yi, yerr=np.ones(nd)*istd, fmt='b+')
    plt.subplot(2, 1, 2)
    plt.title("Same data assuming equidistant sampling")
    plt.plot(te, gi, 'r-')
    plt.errorbar(te, yi, yerr=np.ones(nd)*istd, fmt='r+')
    plt.tight_layout()
    plt.show()
