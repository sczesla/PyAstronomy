Line list based Gaussian spectral model
========================================

.. p23ready

.. currentmodule:: PyAstronomy.modelSuite

.. autoclass:: LLGauss
   :members:

Example: Evaluation and fitting
-------------------------------   

::
  
    from PyAstronomy import modelSuite as ms
    import numpy as np
    import matplotlib.pylab as plt
    
    # Create our line list with 4 line
    lineList = np.zeros((4,3))
    # Assign wavelengths (in A)
    lineList[0,0] = 5002.37
    lineList[1,0] = 5005.9
    lineList[2,0] = 5007.52
    lineList[3,0] = 5007.64
    # Assign EWs (in A)
    lineList[0,1] = 0.01
    lineList[1,1] = 0.05
    lineList[2,1] = 0.009
    lineList[3,1] = 0.12
    # Assign depths (0-1)
    lineList[0,2] = 0.97
    lineList[1,2] = 0.9
    lineList[2,2] = 0.99
    lineList[3,2] = 0.35
    
    wvl = np.arange(5000., 5010., 0.01)
    
    # Get an instance of the LLGauss class
    llg = ms.LLGauss(lineList)
    # Have a look at the model parameters
    llg.parameterSummary()
    # Evaluate the model
    m1 = llg.evaluate(wvl)
    # Now apply rotational broadening [km/s]
    # with limb-darkening of 0.6
    llg["vsini"] = 61.0
    llg["eps"] = 0.6
    # and evaluate again
    mvsini = llg.evaluate(wvl)
    # Next, apply a Doppler shift [km/s]
    llg["vrad"] = -32.7
    # and evaluate
    mvrad = llg.evaluate(wvl)
    # Plot the results
    plt.subplot(2,1,1)
    plt.plot(wvl, m1, 'b.-')
    plt.plot(wvl, mvsini, 'g.-')
    plt.plot(wvl, mvrad, 'y.-')
    
    # Now use the model for fitting
    # We need "data" ...
    data = llg.evaluate(wvl)
    # ... with noise
    data += np.random.normal(0.0, 0.01, len(data))
    # Lets modify the strengths of the Gaussians
    # and get it back.
    for i in range(llg.numberOfLines()):
      llg["A"+str(i+1)] += np.random.normal(0.0, 0.1)
    # Use all line strengths for fitting
    llg.thawLineStrengths()
    # and fit
    llg.fit(wvl, data)
    # Plot the result
    plt.subplot(2,1,2)
    plt.errorbar(wvl, data, yerr=np.ones(len(wvl))*0.01, fmt='bp')
    plt.plot(wvl, llg.evaluate(wvl), 'r--')
    plt.show()
