
Circular orbit---radial-velocity shift
==============================================

.. p23ready

The radial velocity is the radial component of the velocity of a source
relative to an observer and is usually inferred spectroscopically.
Light from an object being part of a binary system and orbiting the common
center of mass will be subject to the Doppler effect. Given a circular orbit,
its radial velocity shift will show sinusoidal variations
determined by the object's orbital elements.

.. currentmodule:: PyAstronomy.modelSuite.radVel
.. autoclass:: SinRadVel
   :members:




Example code - Fit RV curve
------------------------------

::
    
    from __future__ import print_function, division
    # Import some unrelated modules
    from numpy import arange, random, ones
    import matplotlib.pylab as plt
    # ... and now the radVel module
    from PyAstronomy.modelSuite import radVel as rv
    
    # Create Radial Velocity object
    r = rv.SinRadVel()
    # Set parameters
    r.assignValue({"P":1.8, "T0":0.25, "K": 0.5, "rv0": 10.0})
    # Choose some time axis and calculate model
    time = arange(100)/100.0 * 3.0 - 1.5
    y = r.evaluate(time)
    
    # Create some faked data by adding noise
    rvData = y + random.normal(0.0, 0.05, y.size)
    
    # Randomize starting parameters for fit
    for p, v in r.parameters().items():
      r[p] = v + (random.random() - 0.5) * v
    # Show starting values
    print("Starting values for fit:")
    r.parameterSummary()
    
    # Thaw all parameters
    r.thaw(list(r.parameters().keys()))
    # Start the fit
    r.fit(time, rvData, yerr=ones(y.size)*0.05)
    
    # Show fit results
    print("Fitted values:")
    r.parameterSummary()
    
    # Let's see what happened...
    plt.ylabel("Radial velocity [km/s]")
    plt.xlabel("Radial velocity [d]")
    plt.errorbar(time, rvData, yerr=ones(y.size)*0.05, fmt='b.')
    plt.plot(time, y, 'r-')
    plt.show()


