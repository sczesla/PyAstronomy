Examples
=============

.. p23ready

In the following, we provide examples demonstrating the use of the
Rossiter McLaughin model.

.. note:: For a more detailed introduction on the fitting see the funcFit tutorial.

Calculating models
----------------------

The following example demonstrates how to calculate models.

Circular orbit
~~~~~~~~~~~~~~~~~

::
    
    # Import some unrelated modules
    from numpy import arange, pi
    import matplotlib.pylab as plt
    # ... and the model suite
    from PyAstronomy import modelSuite as ms
    
    # Create Rossiter-McLaughlin object
    rmcl = ms.RmcL()
    # Set parameters
    rmcl.assignValue({"a": 6.7, "lambda": 7.2/180.0*pi, "epsilon": 0.5,
                      "P": 1.74, "T0": 0.2, "i": 87.8/180.*pi,
                      "Is": 90.0/180.0*pi, "Omega": 1.609e-5, "gamma": 0.2})
    # Choose some time axis and calculate model
    time = arange(100)/100.0 * 0.2 + 0.1
    rv = rmcl.evaluate(time)
    
    # Let's see what happened...
    plt.ylabel("Radial velocity [stellar-radii/s]")
    plt.xlabel("Time [d]")
    plt.plot(time, rv, '.')
    plt.show()


Elliptical orbit
~~~~~~~~~~~~~~~~~~~~~

::
    
    # Import some unrelated modules
    from numpy import arange, pi
    import matplotlib.pylab as plt
    # ... and the model suite
    from PyAstronomy import modelSuite as ms
    
    # Create Rossiter-McLaughlin object (circular orbit)
    rmcl = ms.RmcL()
    # and one for an elliptical orbit
    rmel = ms.RmcLell()
    
    # Assign parameter values
    rmcl.assignValue({"a": 6.7, "lambda": 7.2/180.0*pi, "epsilon": 0.5,
                      "P": 1.74, "T0": 0.2, "i": 87.8/180.*pi,
                      "Is": 90.0/180.0*pi, "Omega": 1.609e-5, "gamma": 0.2})
    rmel.assignValue({"a": 6.7, "lambda": 7.2/180.0*pi, "epsilon": 0.5,
                      "P": 1.74, "tau": 0.2, "i": 87.8/180.*pi, "w": -90/180.*pi,
                      "e": 0.05, "Is": 90.0/180.0*pi, "Omega": 1.609e-5, "gamma": 0.2})
    
    
    # Choose some time axis and calculate model
    time = arange(100)/100.0 * 0.2 + 0.1
    rvc = rmcl.evaluate(time)
    rve = rmel.evaluate(time)
    
    # Let's see what happened...
    plt.ylabel("Radial velocity [stellar-radii/s]")
    plt.xlabel("Time [d]")
    plt.plot(time, rvc, 'b.-', label="circular")
    plt.plot(time, rve, 'r.-', label="elliptical")
    plt.legend()
    plt.show()



Fitting a model
----------------

This example is an extension of the first. It demonstrates how a model fit can
be carried out. 

::
    
    # Import some unrelated modules
    from numpy import arange, pi, random
    import matplotlib.pylab as plt
    # ... and the model suite
    from PyAstronomy import modelSuite as ms
    
    # Create Rossiter-McLaughlin object
    rmcl = ms.RmcL()
    # Set parameters
    rmcl.assignValue({"a": 6.7, "lambda": 7.2/180.0*pi, "epsilon": 0.5,
                      "P": 1.74, "T0": 0.2, "i": 87.8/180.*pi,
                      "Is": 90.0/180.0*pi, "Omega": 1.609e-5, "gamma": 0.2})
    # Choose some time axis and calculate model
    time = arange(100)/100.0 * 0.2 + 0.1
    rv = rmcl.evaluate(time)
    
    # Add some noise.
    rv += random.normal(0.0, 0.05*rv.max(), rv.size)
    
    # Assign guess parameters
    rmcl.assignValue({"a": 6.0, "lambda": 7.2/180.0*pi, "epsilon": 0.5,
                      "P": 1.74, "T0": 0.17, "i": 87.8/180.*pi,
                      "Is": 90.0/180.0*pi, "Omega": 1.609e-5, "gamma": 0.2})
    
    # Thaw parameters and fit
    rmcl.thaw(["a", "T0"])
    rmcl.fit(time, rv)
    
    # Investigate the outcome
    rmcl.parameterSummary()
    
    # Let's see what happened...
    plt.ylabel("Radial velocity [stellar-radii/s]")
    plt.xlabel("Time [d]")
    plt.plot(time, rv, '.')
    plt.plot(time, rmcl.model, 'r--')
    plt.legend(["Observation", "Model"])
    plt.show()
