.. _tutMCMCSampler:

.. _matplotlib: http://matplotlib.sourceforge.net/
.. _pymc: https://github.com/pymc-devs/pymc
.. _SciPy: www.scipy.org/
.. _numpy: numpy.scipy.org/
.. _XSPEC: http://heasarc.nasa.gov/xanadu/xspec/
.. _emcee: http://dan.iel.fm/emcee/current/

MCMC sampling with funcFit tutorial
======================================

Currently, funcFit supports MCMC sampling either via the pymc_ or the emcee_ package.
To do this, the model objects provides the `fitMCMC` method (pymc) and the `fitEMCEE`
method (emcee). Both are basically independent and can be used separately.

Although the method names start with *fit*, we emphasize that
fitting, in the sense of optimization, is not the exact purpose of the analysis
carried out here, but sampling from the posterior.

.. note:: The Markov chains produced by fitMCMC and fitEMCEE can be analyzed using the TraceAnalysis class.
          See this tutorial: :ref:`traceAnalysisTut`.

.. contents:: Sections of this tutorial

Sampling using pymc (fitMCMC)
--------------------------------------------

The *fitMCMC* method provided by funcFit is not an MCMC sampler itself, but
it is a wrapper around functionality provided by a third party package,
namely, PyMC_.

pymc is a powerful Python package providing a wealth of functionality concerning Bayesian
analysis. *fitMCMC* provides an easy to use interface to pymc sampling, which
allows to carry out a basic Bayesian data analysis quickly. 



.. note:: To run these examples, pymc must be installed (check the output of funcFit.status() shown at the beginning
          of this tutorial to see whether this is the case on your system).
          
.. warning:: Unfortunately, pymc is only supported in version 2.x. In particular, version 3.x is not supported.

Using the Markov-Chain Monte-Carlo (MCMC) sampler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following example demonstrates
how the funcFit interface can be used to carry out a Bayesian analysis using pymc.
For a deeper understanding of the working, adaptability, and logic implemented by
pymc, we refer the reader to their web page (PyMC_).

::
    
    from __future__ import print_function, division
    # Import some required modules
    from numpy import arange, sqrt, exp, pi, random, ones
    import matplotlib.pylab as plt
    import pymc
    # ... and now the funcFit package
    from PyAstronomy import funcFit as fuf
    
    # Creating a Gaussian with some noise
    # Choose some parameters...
    gPar = {"A":-5.0, "sig":10.0, "mu":10.0, "off":1.0, "lin":0.0}
    # Calculate profile
    x = arange(100) - 50.0
    y = gPar["off"] + gPar["A"] / sqrt(2*pi*gPar["sig"]**2) \
    * exp(-(x-gPar["mu"])**2/(2*gPar["sig"]**2))
    # Add some noise
    y += random.normal(0.0, 0.01, x.size)
    
    # Now let us come to the fitting
    # First, we create the Gauss1d fit object
    gf = fuf.GaussFit1d()
    # See what parameters are available
    print("List of available parameters: ", gf.availableParameters())
    # Set guess values for the parameters
    gf["A"] = -10.0
    gf["sig"] = 15.77
    gf["off"] = 0.87
    gf["mu"] = 7.5
    # Let us see whether the assignment worked
    print("Parameters and guess values: ", gf.parameters())
    
    # Which parameters shall be variable during the fit?
    # 'Thaw' those (the order is irrelevant)
    gf.thaw(["A", "sig", "off", "mu"])
    
    # Now start a simplex fit
    gf.fit(x,y,yerr=ones(x.size)*0.01)
    
    # Obtain the best-fit values derived by the simplex fit.
    # They are to be used as start values for the MCMC sampling.
    # Note that 'A' is missing - we will introduce this later.
    X0 = {"sig":gf["sig"], "off":gf["off"], "mu":gf["mu"]}
    
    # Now we specify the limits within which the individual parameters
    # can be varied (for those parameters listed in the 'X0' dictionary).
    Lims = {"sig":[-20.,20.], "off":[0.,2.], "mu":[5.,15.]}
    
    # For the parameters contained in 'X0', define the step widths, which
    # are to be used by the MCMC sampler. The steps are specified using
    # the same scale/units as the actual parameters.
    steps = {"A":0.01, "sig":0.1, "off":0.1, "mu":0.1}
    
    # In this example, we wish to define our ``own'' PyMC variable for the parameter
    # 'A'. This can be useful, if nonstandard behavior is desired. Note that this
    # is an optional parameter and you could simply include the parameter 'A' into
    # The framework of X0, Lims, and steps.
    ppa = {}
    ppa["A"] = pymc.Uniform("A", value=gf["A"], lower=-20., \
                            upper=10.0, doc="Amplitude")
    
    # Start the sampling. The resulting Marchov-Chain will be written
    # to the file 'mcmcExample.tmp'. In default configuration, pickle
    # is used to write that file.
    # To save the chain to a compressed 'hdf5'
    # file, you have to specify the dbArgs keyword; e.g., use:
    #   dbArgs = {"db":"hdf5", "dbname":"mcmcExample.hdf5"}
    gf.fitMCMC(x, y, X0, Lims, steps, yerr=ones(x.size)*0.01, \
                pymcPars=ppa, iter=2500, burn=0, thin=1, \
                dbfile="mcmcExample.tmp")
    
    # Reload the database (here, this is actually not required, but it is
    # if the Marchov chain is to be analyzed later).
    db = pymc.database.pickle.load('mcmcExample.tmp')
    # Plot the trace of the amplitude, 'A'.
    plt.hist(db.trace("A", 0)[:])
    plt.show()


Some points shall be emphasized in this example:
 * For MCMC sampling the exact same fit object is used as for "normal" fitting.
 * If the *yerr* keyword is specified in the call to fitMCMC, a Gaussian distribution is
   assumed for the data points. Otherwise a Poisson distribution is assumed.
 * We used the normal simplex fit to obtain starting values for the Markov chain. You may also
   use, e.g., burn-in.
 * In the example, we demonstrated how a uniformly distributed PyMC variable is created. Normally,
   the fitMCMC method does this for you.
 * The result, i.e., the Markov chain, is saved to the file mcmcExample.tmp and is reloaded to
   obtain the trace of the amplitude.

.. note:: A convenient analysis of the resulting traces can be carried out using the *TraceAnalysis*
          class (see :ref:`traceAnalysisClass`)

We emphasize that PyMC is a powerful and highly adaptable package, which can do a lot more. A more
detailed introduction is, however, beyond the scope of this tutorial.



Simplifying initialization (MCMCautoParameters)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It can become cumbersome to define the starting values, steps,
and ranges for uniform priors as done in the above example.
Using the "auto" methods defined in the fitting class, you can take a short
cut. However, be warned:

.. warning:: There is NO guarantee that the auto functions produce reasonable
             results. You need to check that.
             
::
    
    from __future__ import print_function, division
    from PyAstronomy import funcFit as fuf
    import numpy as np
    import matplotlib.pylab as plt
    
    x = np.linspace(0,30,1000)
    gauss = fuf.GaussFit1d()
    gauss["A"] = 1
    gauss["mu"] = 23.
    gauss["sig"] = 0.5
    # Generate some "data" to fit
    yerr = np.random.normal(0., 0.05, len(x))
    y = gauss.evaluate(x) + yerr
    # Thaw the parameters A, mu, and sig
    gauss.thaw(["A","mu","sig"])
    
    # Define the ranges, which are used to construct the
    # uniform priors and step sizes.
    # Note that for "sig", we give only a single value.
    # In this case, the limits for the uniform prior will
    # be constructed as [m0-1.5, m0+1.5], where m0 is the
    # starting value interpreted as the current value of
    # mu (23. in this case).
    ranges = {"A":[0,10],"mu":3, "sig":[0.1,1.0]}
    # Generate default input for X0, lims, and steps
    X0, lims, steps = gauss.MCMCautoParameters(ranges)
    
    # Show what happened...
    print()
    print("Auto-generated input parameters:")
    print("X0: ", X0)
    print("lims: ", lims)
    print("steps: ", steps)
    print()
    # Call the usual sampler
    gauss.fitMCMC(x, y, X0, lims, steps, yerr=yerr, iter=1000)
    
    # and plot the results
    plt.plot(x, y, 'k+')
    plt.plot(x, gauss.evaluate(x), 'r--')
    plt.show()


You may even shorten the short-cut by using the `autoFitMCMC` method.
**However, note that the same warning remains valid here.**

::

  from PyAstronomy import funcFit as fuf
  import numpy as np
  import matplotlib.pylab as plt
  
  x = np.linspace(0,30,1000)
  gauss = fuf.GaussFit1d()
  gauss["A"] = 1
  gauss["mu"] = 23.
  gauss["sig"] = 0.5
  # Generate some "data" to fit
  yerr = np.random.normal(0., 0.05, len(x))
  y = gauss.evaluate(x) + yerr
  
  # Define the ranges, which are used to construct the
  # uniform priors and step sizes.
  # Note that for "sig", we give only a single value.
  # In this case, the limits for the uniform prior will
  # be constructed as [m0-1.5, m0+1.5], where m0 is the
  # starting value interpreted as the current value of
  # mu (23. in this case).
  ranges = {"A":[0,10],"mu":3, "sig":[0.1,1.0]}
  
  # Call the auto-sampler
  # Note that we set picky to False here. In this case, the
  # parameters specified in ranges will be thawed automatically.
  # All parameters not mentioned there, will be frozen.
  gauss.autoFitMCMC(x, y, ranges, yerr=yerr, picky=False, iter=1000)
  
  # and plot the results
  plt.plot(x, y, 'k+')
  plt.plot(x, gauss.evaluate(x), 'r--')
  plt.show()


Sampling with Gaussian and uniform priors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The use of prior information is inherent in Bayesian analyses. The
following example demonstrates how prior information can explicitly
be included in the sampling. We note, however, that some kind of prior
is implicitly assumed for all parameters; in this case, a uniform one.

::

  from PyAstronomy import funcFit as fuf
  import numpy as np
  import matplotlib.pylab as plt
  import pymc
  
  # Create a Gauss-fit object
  gf = fuf.GaussFit1d()
  
  # Choose some parameters
  gf["A"] = -0.65
  gf["mu"] = 1.0
  gf["lin"] = 0.0
  gf["off"] = 1.1
  gf["sig"] = 0.2
  
  # Simulate data with noise
  x = np.linspace(0., 2., 100)
  y = gf.evaluate(x)
  y += np.random.normal(0, 0.05, len(x))
  
  gf.thaw(["A", "off", "mu", "sig"])
  
  # Set up a normal prior for the offset parameter
  # Note!---The name (first parameter) must correspond to that
  #         of the parameter.
  # The expectation value us set to 0.9 while the width is given
  # as 0.01 (tau = 1/sigma**2). The starting value is specified
  # as 1.0.
  offPar = pymc.Normal("off", mu=0.9, tau=(1./0.01)**2, value=1.0)
  # Use a uniform prior for mu.
  muPar = pymc.Uniform("mu", lower=0.95, upper=0.97, value=0.96)
  
  # Collect the "extra"-variables in a dictionary using
  # their names as keys
  pymcPars = {"mu":muPar, "off":offPar}
  
  # Specify starting values, X0, and limits, lims, for
  # those parameter distributions not given specifically.
  X0 = {"A":gf["A"], "sig":gf["sig"]}
  lims = {"A":[-1.0,0.0], "sig":[0., 1.0]}
  # Still, the steps dictionary has to contain all
  # parameter distributions.
  steps = {"A":0.02, "sig":0.02, "mu":0.01, "off":0.01}
  
  # Carry out the MCMC sampling
  gf.fitMCMC(x, y, X0, lims, steps, yerr=np.ones(len(x))*0.05, \
             pymcPars=pymcPars, burn=1000, iter=3000)
  
  # Setting parameters to mean values
  for p in gf.freeParameters():
    gf[p] = gf.MCMC.trace(p)[:].mean()
  
  # Show the "data" and model in the upper panel
  plt.subplot(2,1,1)
  plt.title("Data and model")
  plt.errorbar(x, y, yerr=np.ones(len(x))*0.05, fmt="bp")
  # Plot lowest deviance solution
  plt.plot(x, gf.evaluate(x), 'r--')
  
  # Show the residuals in the lower panel
  plt.subplot(2,1,2)
  plt.title("Residuals")
  plt.errorbar(x, y-gf.evaluate(x), yerr=np.ones(len(x))*0.05, fmt="bp")
  plt.plot([min(x), max(x)], [0.0,0.0], 'r-')
  
  plt.show()

Clearly, the plot shows that the solution does not fit well in a Chi-square sense,
because the prior information has a significant influence on the outcome.
Whether this should be considered reasonable or not is not a question the sampler could answer.



Sampling using emcee (fitEMCEE)
---------------------------------------------

The emcee_ package relies on sampling with an ensemble of chains, the so-called walkers. In this
way, the sampling is automatically adapted to the scale of the problem, which simplifies
obtained reasonable acceptance rates.


Basic sampling with emcee
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following example shows a basic application of fitEMCEE. By default,
the resulting Markov chain is saved to a file called 'chain.emcee'.

::

    # Import numpy and matplotlib
    from numpy import arange, sqrt, exp, pi, random, ones
    import matplotlib.pylab as plt
    # ... and now the funcFit package
    from PyAstronomy import funcFit as fuf
    
    # Before we can start fitting, we need something to fit.
    # So let us create some data...
    
    # Choose some signal-to-noise ratio
    snr = 25.0
    
    # Creating a Gaussian with some noise
    # Choose some parameters...
    gf = fuf.GaussFit1d()
    gf.assignValues({"A":-5.0, "sig":2.5, "mu":10.0, "off":1.0, "lin":0.0})
    # Calculate profile
    x = arange(100) - 50.0
    y = gf.evaluate(x)
    # Add some noise
    y += random.normal(0.0, 1.0/snr, x.size)
    
    # Define the free parameters
    gf.thaw(["A", "sig", "mu", "off"])
    
    # Start a fit (quite dispensable here)
    gf.fit(x, y, yerr=ones(x.size)/snr)
    
    # Say, we want 200 burn-in iterations and, thereafter,
    # 1000 further iterations (per walker).
    sampleArgs = {"iters":1000, "burn":200}
    
    # Start the sampling (ps could be used to continue the sampling)
    ps = gf.fitEMCEE(x, y, yerr=ones(x.size)/snr, sampleArgs=sampleArgs)
    
    # Plot the distributions of the chains
    # NOTE: the order of the parameters in the chain object is the same
    #       as the order of the parameters returned by freeParamNames()
    for i, p in enumerate(gf.freeParamNames()):
      plt.subplot(len(gf.freeParamNames()), 1, i+1)
      plt.hist(gf.emceeSampler.flatchain[::,i], label=p)
      plt.legend()
    plt.show()
    
  
Sampling with prior information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, we use a very simple constant model to, first, compare
the result from sampling with classical error estimation; note that in this
simple case Bayesian credibility intervals are, indeed, numerically identical
to classical confidence intervals, which is, however, not generally the case.
Second, we introduce (strong) prior information and repeat the sampling.

A number of ready-to-use priors are implemented here: :ref:`emceePriors`
    
::    

    from __future__ import print_function, division
    # Import numpy and matplotlib
    from numpy import arange, sqrt, exp, pi, random, ones
    import matplotlib.pylab as plt
    # ... and now the funcFit package
    from PyAstronomy import funcFit as fuf
    import numpy as np
    
    # Before we can start fitting, we need something to fit.
    # So let us create some data...
    
    # Choose some signal-to-noise ratio
    snr = 25.0
    
    # Choosing an arbitrary constant and ...
    c = 10.0
    # ... an equally arbitrary number of data points
    npoint = 10
    
    # Define 'data'
    x = arange(npoint)
    y = np.ones(len(x)) * c
    # Add some noise
    y += random.normal(0.0, 1.0/snr, x.size)
    
    # A funcFit object representing a constant
    pf = fuf.PolyFit1d(0)
    pf["c0"] = c
    
    # The only parameter shall be free
    pf.thaw("c0")
    
    # Say, we want 200 burn-in iterations and, thereafter,
    # 2500 further iterations (per walker).
    sampleArgs = {"iters":2500, "burn":200}
    
    # Start the sampling (ps could be used to continue the sampling)
    ps = pf.fitEMCEE(x, y, yerr=ones(x.size)/snr, sampleArgs=sampleArgs)
    print()
    
    # Plot the distributions of the chains
    # NOTE: the order of the parameters in the chain object is the same
    #       as the order of the parameters returned by freeParamNames()
    h = plt.hist(pf.emceeSampler.flatchain[::,0], label="c0", normed=True)
    # Construct "data points" in the middle of the bins
    xhist = (h[1][1:] + h[1][0:-1]) / 2.0
    yhist = h[0]
    
    # Fit the histogram using a Gaussian
    gf = fuf.GaussFit1d()
    gf.assignValues({"A":1.0, "mu":c, "sig":1.0/snr/np.sqrt(npoint)})
    # First fitting only "mu" is simply quite stable 
    gf.thaw("mu")
    gf.fit(xhist, yhist)
    gf.thaw(["A", "sig"])
    gf.fit(xhist, yhist)
    
    print()
    print("  --- Sampling results ---")
    print("Posterior estimate of constant: ", np.mean(pf.emceeSampler.flatchain[::,0]))
    print("Nominal error of the mean: ", 1.0/snr/np.sqrt(npoint))
    print("Estimate from Markov chain: ", np.std(pf.emceeSampler.flatchain[::,0]), end=' ')
    print(" and from Gaussian fit to distribution: ", gf["sig"])
    
    # Evaluate best-fit model ...
    xmodel = np.linspace(c - 10.0/snr, c + 10.0/snr, 250)
    ymodel = gf.evaluate(xmodel)
    # ... and plot
    plt.plot(xhist, yhist, 'rp')
    plt.plot(xmodel, ymodel, 'r--')
    plt.legend()
    plt.show()
    
    
    # Defining a prior on c0. Prior knowledge tells us that its value
    # is around 7. Let us choose the standard deviation of the prior so
    # that the estimate will lie in the middle between 7 and 10. Here we
    # exploit symmetry and make the prior information as strong as the
    # information contained in the likelihood function.
    priors = {"c0":fuf.FuFPrior("gaussian", sig=1.0/snr/np.sqrt(npoint), mu=7.0)}
    
    # Start the sampling (ps could be used to continue the sampling)
    ps = pf.fitEMCEE(x, y, yerr=ones(x.size)/snr, sampleArgs=sampleArgs, priors=priors)
    
    print()
    print("  --- Sampling results with strong prior information ---")
    print("Posterior estimate of constant: ", np.mean(pf.emceeSampler.flatchain[::,0]), end=' ')
    print(" +/-", np.std(pf.emceeSampler.flatchain[::,0]))
    
    plt.hist(pf.emceeSampler.flatchain[::,0], label="c0", normed=True)
    plt.show()
    
