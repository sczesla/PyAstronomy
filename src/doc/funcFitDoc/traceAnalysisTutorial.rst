.. _traceAnalysisTut:

Analyze Markov-Chains using `TraceAnalysis`
=============================================

.. currentmodule:: PyAstronomy.funcFit

After having using `pymc`, e.g., via PyAstronomy's
wrapper `fitMCMC` and having created a Markov Chain,
the information needs to be extracted from that Chain.
:py:class:`TraceAnalysis` is a helper class designed to help you
doing this.

This collection of examples demonstrates how to create
a Markov Chain and carry out a basic analysis.

Create a Markov-Chain to work with
-------------------------------------

The following example demonstrates how to **create
the Markov-Chain**, which will later be analyzed. 

::
    
    from __future__ import print_function, division
    # Import some required modules
    import numpy as np
    import matplotlib.pylab as plt
    # ... and now the funcFit package
    from PyAstronomy import funcFit as fuf
    
    # Starting from with Voigt profile
    vp = fuf.Voigt1d()
    # Set some values to create a model
    vp["A"] = -0.4
    vp["al"] = 0.7
    vp["mu"] = 5500.
    vp["ad"] = 0.3
    vp["off"] = 1.0
    
    x = np.linspace(5490., 5510., 200)
    # Create our data with some noise
    yerr = np.ones(len(x))*0.01
    y = vp.evaluate(x) + np.random.normal(0.0, 0.01, len(x))
    
    # Say, we have a guess of the parameters, which is, however,
    # not entirely correct
    vp["A"] = -0.376
    vp["al"] = 0.9
    vp["mu"] = 5499.7
    vp["ad"] = 0.4
    vp["off"] = 1.0
    
    # Plot the data and our guess
    plt.errorbar(x, y, yerr=yerr, fmt='b.-')
    plt.plot(x, vp.evaluate(x), 'r--')
    plt.show()
    
    # Thaw the parameters, which we wish to vary
    # during the sampling
    vp.thaw(["A", "al", "mu", "ad"])
    
    # Use current parameters as starting point for the sampling
    X0 = vp.freeParameters()
    print("Starting point for sampling: ", X0)
    
    
    ##############################################################
    ### Using pymc for sampling, for emcee see below
    ##############################################################
    
    # Now we specify the limits within which the individual parameters
    # can be varied. Actually, you specify the limits of uniform priors
    # here.
    lims = {"A":[-1.0,0.0], "al":[0.0,3.], "ad":[0.0,3.0], "mu":[5495., 5505.]}
    
    # Provide a guess for the proposal step widths.
    # Try to guess the scale of the problem in the individual
    # parameters.
    steps = {"A":0.02, "al":0.01, "ad":0.01, "mu":0.05}
    
    # Start the sampling. The resulting Marchov-Chain will be written
    # to the file 'mcmcTA.tmp'. In default configuration, pickle
    # is used to write that file.
    # To save the chain to a compressed 'hdf5'
    # file, you have to specify the dbArgs keyword; e.g., use:
    #   dbArgs = {"db":"hdf5", "dbname":"mcmcExample.hdf5"}
    vp.fitMCMC(x, y, X0, lims, steps, yerr=yerr, \
               iter=2500, burn=0, thin=1, \
               dbfile="mcmcTA.tmp")
    
    ##############################################################
    ### Uncomment to use emcee for sampling
    ##############################################################           
               
    #priors = {"al":fuf.FuFPrior("limuniform", lower=0.0, upper=100.), \
              #"ad":fuf.FuFPrior("limuniform", lower=0.0, upper=100.)}
    
    ## Note that the filename should end in .emcee. Substitute this filename
    ## in the following examples.
    #vp.fitEMCEE(x, y, yerr=yerr, sampleArgs={"iters":2500}, \
                #dbfile="mcmcTA.emcee", priors=priors)

Investigate convergence behavior using the "deviance"
------------------------------------------------------

This example demonstrates how to investigate whether the
chain has converged using the **deviance** and setting
an appropriate **burn-in**.

::

  from PyAstronomy import funcFit as fuf
  
  # Create an instance of TraceAnalysis
  # telling it which file to use
  ta = fuf.TraceAnalysis("mcmcTA.tmp")
  
  # Have a look at the deviance to check if and when
  # the chains reached equilibrium.
  ta.plotTrace("deviance")
  ta.show()
  
  # Say, we are sure that after 500 iterations, the chain
  # reached equilibrium. We use this as the burn-in phase
  ta.setBurn(500)
  
  # Have a second look at the deviance, this time considering
  # the burn-in. Note that the first 500 iterations are not
  # removed from the chain. They are just not considered any
  # more.
  ta.plotTrace("deviance")
  ta.show()


Accessing the traces
----------------------

Here, we show how to **access the traces** and calculate
**mean, median, standard deviation, and Highest Probability
Density (HPD) credibility intervals**.

::

    from __future__ import print_function, division
    from PyAstronomy import funcFit as fuf
    
    # Create an instance of TraceAnalysis
    # telling it which file to use
    ta = fuf.TraceAnalysis("mcmcTA.tmp")
    
    # Use the burn-in from the previous example
    ta.setBurn(500)
    
    # See which model parameters have been sampled
    print("Available parameters: ", ta.availableParameters())
    
    # Access the traces of these parameters
    print("Trace for A: ", ta["A"])
    
    # Calculate mean, median, standard deviation, and
    # credibility interval for the available parameters
    for p in ta.availableParameters():
      hpd = ta.hpd(p, cred=0.95)
      print("Parameter %5s, mean = % g, median = % g, std = % g, 95%% HPD = % g - % g" \
            % (p, ta.mean(p), ta.median(p), ta.std(p), hpd[0], hpd[1]))
    


Analyze parameter correlations
-------------------------------

**Parameter correlations** are important in many
analysis. This example shows how such analyses can
be studied using the Markov-Chain.  

::

    from __future__ import print_function, division
    from PyAstronomy import funcFit as fuf
    
    # Create an instance of TraceAnalysis
    # telling it which file to use
    ta = fuf.TraceAnalysis("mcmcTA.tmp")
    
    # Use the burn-in from the previous example
    ta.setBurn(500)
    
    # Have a look at the parameter correlations
    ta.correlationTable()
    
    # Calculate Pearson's and Spearman's r-coefficients
    print("Pearson: ", ta.pearsonr("ad", "al"))
    print("Spearman: ", ta.spearmanr("ad", "al"))
    
    # Show a plot of the correlation
    # Note that the plotCorrEnh method can also
    # be used, which is useful in the case of long
    # chains.
    ta.plotCorr(parsList=["ad", "al"])
    # ta.plotCorrEnh(parsList=["ad", "al"])
    ta.show()
   

Obtain lowest-deviance, mean, and median set of parameters
------------------------------------------------------------

This example demonstrates how to obtain certain parameter
sets and use them with a model.

::
    
    from __future__ import print_function, division
    from PyAstronomy import funcFit as fuf
    import matplotlib.pylab as plt
    import numpy as np
    
    # Create an instance of TraceAnalysis
    # telling it which file to use
    ta = fuf.TraceAnalysis("mcmcTA.tmp")
    
    # Use the burn-in from the previous example
    ta.setBurn(500)
    
    # Find sets of parameters
    # First, the lowest deviance set
    lds, index = ta.parameterSet(prescription="lowestDev")
    print("Lowest deviance set: ", lds)
    print("  at chain index: ", index)
    means = ta.parameterSet(prescription="mean")
    print("Set of mean values: ", means)
    medians = ta.parameterSet(prescription="median")
    print("Set of median values: ", means)
    
    # Create Voigt model and plot the models belonging
    # to the lowest deviance, mean, and median parameter
    # set.
    vp = fuf.Voigt1d()
    # Generate the model wavelength axis
    x = np.linspace(5490., 5510., 200)
    # Calculate and plot the models
    vp.assignValues(lds)
    plt.plot(x, vp.evaluate(x), 'b.-')
    vp.assignValues(means)
    plt.plot(x, vp.evaluate(x), 'r.-')
    vp.assignValues(medians)
    plt.plot(x, vp.evaluate(x), 'g.-')
    plt.show()


Other useful visualization utilities
---------------------------------------

This example shows a couple of other useful
routines to visualize and analyze the
Markov-Chains. 

::
  
    from PyAstronomy import funcFit as fuf
    
    # Create an instance of TraceAnalysis
    # telling it which file to use
    ta = fuf.TraceAnalysis("mcmcTA.tmp")
    
    # Use the burn-in from the previous example
    ta.setBurn(500)
    
    # Investigate a trace
    ta.plotTrace("mu")
    ta.show()
    # and its distribution.
    ta.plotHist("mu")
    ta.show()
    # Combine trace and distribution
    ta.plotTraceHist("mu")
    ta.show()
    # Plot correlations
    ta.plotCorr(parsList=["mu", "ad", "al"])
    ta.show()
    