The Kepler-ellipse model
=========================

.. p23ready

.. currentmodule:: PyAstronomy.modelSuite

The *KeplerEllipseModel* class provides a fitting model for the Kepler orbit.
In particular, the model provides the mapping from time to Cartesian x, y, and
z coordinates. On creation of the model, one can choose which of the coordinate
axis are of interest to the user. The coordinates will then be returned as a
one-dimensional array as explained the `evaluate` method of the class.

.. autoclass:: KeplerEllipseModel
   :members:

Example: Evaluating and fitting model using MCMC
------------------------------------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy.modelSuite import KeplerEllipseModel
    import numpy as np
    import matplotlib.pylab as plt
    
    # Create a model class instance
    # In this case, we are only interested
    # in the x- and z-components of the orbit
    # solution.
    kem = KeplerEllipseModel(relevantAxes="xz")
    
    # Setting some guess parameters
    kem["a"] = 7.8
    kem["per"] = 12.3
    kem["e"] = 0.07
    kem["tau"] = 0.745
    kem["Omega"] = 143.
    kem["w"] = 0.2
    kem["i"] = 92.0
    
    # Evaluate the model
    time = np.linspace(0, kem["per"], 20)
    model = kem.evaluate(time)
    # Note that the model has twice the number of points
    # compared to the time axis. This is because it contains
    # the data for two axes
    print("Used " + str(len(time)) + " time points")
    print("-> length of model: ", len(model))
    
    # Isolating the model for the x-axis, i.e.,
    # every second data point starting from the
    # beginning.
    xmodel = model[0::2]
    # Isolating the model for the y-axis
    ymodel = model[1::2]
    
    # Use the model to obtain mock data
    # by introducing some scatter
    data = model + np.random.normal(0., 0.5, model.size)
    # Plot the resulting "data"
    plt.title("Kepler Ellipse Model --- Example")
    plt.errorbar(data[0::2], data[1::2], xerr=np.ones(20)*0.5, \
                  yerr=np.ones(20)*0.5, fmt="bp")
    
    # Use MCMC to sample from the posterior
    # Specify free parameters
    kem.thaw(["a", "per", "e", "tau", "Omega", "w", "i"])
    # Specify starting values
    X0 = {}
    steps = {}
    for p in kem.freeParameters():
      X0[p] = kem[p]
      steps[p] = kem[p] / 20.
    
    lims = {"a":[5.,10.], "per":[10.,15.], "e":[0.,1.], "tau":[0.5,1.], \
            "Omega":[0.,360.], "w":[-5.,5.], "i":[90., 95.]}
    
    kem.fitMCMC(time, data, X0, lims, steps, yerr=np.ones(len(data))*0.5, \
                iter=500, dbfile="kemExample.tmp")
    
    # Plot the lowest deviance model
    ldmodel = kem.evaluate(np.linspace(0, kem["per"], 200))
    plt.plot(ldmodel[0::2], ldmodel[1::2], 'r--')
    
    plt.show()
