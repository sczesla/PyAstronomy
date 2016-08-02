Examples
==============

.. p23ready

The following examples should give a reasonably impression of how the transit model
class can be used. 

Calculating a model
-----------------------

In this first example, we demonstrate how to calculate a model transit light-curve. 

::
  
    from PyAstronomy.modelSuite import palTrans
    import matplotlib.pylab as plt
    import numpy as np
    
    # Create a PalLC instance
    plc = palTrans.PalLC()
    
    # Set parameter values
    plc["p"] = 0.1    # Planet radius / Stellar radius
    plc["per"] = 1.0  # Orbital period
    plc["a"] = 2.0    # Large semi major axis [R_S]
    plc["i"] = 90.0   # Orbital inclination [deg]
    # Specify limb darkening
    # (quadratic limb-darkening law)
    plc["linLimb"] = 0.4
    plc["quadLimb"] = 0.2
    # Specify T0 (central time of transit)
    plc["T0"] = -0.1
    # Specify binary contribution
    plc["b"] = 0.0
    
    # Check the parameters
    plc.parameterSummary()
    
    # Create a time axis
    time = np.arange(1000)/1000.0 - 0.5
    
    # Calculate the light curve using above set
    # model parameters
    lightcurve = plc.evaluate(time)
    
    # Plot the result
    plt.plot(time, lightcurve, 'bp')
    plt.show()


Fitting a transit
---------------------

The leading half of this example is very similar to the first one; we only
use a smaller number of point for the time axis. The second half demonstrates
how the fitting interface is invoked (more examples using simpler models
can are shown in the *funcFit tutorial*). 

.. warning:: If you are using *mpmath* for evaluation of the elliptical integrals, the
             calculations will be quite slow.

::
    
    from __future__ import print_function, division
    from PyAstronomy.modelSuite import palTrans
    import matplotlib.pylab as plt
    import numpy as np
    
    # Create a PalLC instance
    plc = palTrans.PalLC()
    
    # Set parameter values
    plc["p"] = 0.1    # Planet radius / Stellar radius
    plc["per"] = 1.0  # Orbital period
    plc["a"] = 7.5    # Large semi major axis [R_S]
    plc["i"] = 90.0   # Orbital inclination [deg]
    # Specify limb darkening
    # (quadratic limb-darkening law)
    plc["linLimb"] = 0.4
    plc["quadLimb"] = 0.2
    # Specify T0 (central time of transit)
    plc["T0"] = -0.1
    # Specify binary contribution
    plc["b"] = 0.0
    
    # Check the parameters
    print("Input parameters: ")
    plc.parameterSummary()
    
    # Create a time axis
    time = np.arange(100)/100.0 * 0.2 - 0.2
    
    # Calculate the light curve using above set
    # model parameters
    lc = plc.evaluate(time)
    
    # Save the result and add some noise
    flux = lc + np.random.normal(0.0, 0.002, time.size)
    
    # Now lets try to recover what we put in
    # Choose some "guess" parameters
    plc["p"] = 0.1     # Planet radius / Stellar radius
    plc["per"] = 1.0   # Orbital period
    plc["a"] = 7.5     # Large semi major axis [R_S]
    plc["i"] = 90.0    # Orbital inclination [deg]
    # Specify limb darkening
    # (quadratic limb-darkening law)
    plc["linLimb"] = 0.4
    plc["quadLimb"] = 0.2
    # Specify T0 (central time of transit)
    plc["T0"] = -0.08
    # Specify binary contribution
    plc["b"] = 0.0
    
    # Assume we want to fit "p", "a", "i", and "T0"
    plc.thaw(["T0", "i"])
    
    # Before we start fitting, check how the elliptical integrals
    # are evaluated (mpmath or Boost)
    print("Which elliptical integrals are used?: ", plc.whichEllInts())
    
    # Carry out the fit
    plc.fit(time, flux, yerr=np.ones(time.size)*0.002)
    
    print("Fit parameters: ")
    plc.parameterSummary()
    
    # Plot the result
    plt.plot(time, flux, 'bp')
    plt.plot(time, plc.model, 'r-')
    plt.show()

 
Obtain a model taking finite integration time into account
-------------------------------------------------------------

This example shows how to use the *PalLC_Rebin* class to take finite
integration times and the resulting light-curve distortion into account.
This example is very similar to the first one.

::
    
    from PyAstronomy.modelSuite import palTrans
    import matplotlib.pylab as plt
    import numpy as np
    
    # Create a PalLC_Rebin instance
    plc = palTrans.PalLC_Rebin()
    
    # Set parameter values
    plc["p"] = 0.1    # Planet radius / Stellar radius
    plc["per"] = 1.0  # Orbital period
    plc["a"] = 2.0    # Large semi major axis [R_S]
    plc["i"] = 90.0   # Orbital inclination [deg]
    # Specify limb darkening
    # (quadratic limb-darkening law)
    plc["linLimb"] = 0.8
    plc["quadLimb"] = 0.2
    # Specify T0 (central time of transit)
    plc["T0"] = -0.1
    # Specify binary contribution
    plc["b"] = 0.0
    
    # Check the parameters
    plc.parameterSummary()
    
    # Create a time axis
    time = np.arange(50)/50.0 - 0.51
    
    # Specify oversampling parameters.
    # Here use 10 points per observed bin.
    plc.setRebinArray_Ndt(time, 10, time[1]-time[0])
    
    # Calculate the light curve using above set
    # model parameters
    lc = plc.evaluate(time)
    
    # Plot the result (both the overbinned and final
    # model light-curves)
    plt.plot(plc.rebinTimes, plc.unbinnedModel, 'b.-')
    plt.plot(time, plc.binnedModel, 'rd--')
    plt.legend(["Overbinned LC", "Averaged LC"])
    plt.show()

