The Kepler ellipse models
===========================

.. p23ready

.. currentmodule:: PyAstronomy.modelSuite

The Kepler ellipse models provide a fitting framework for the Kepler orbit.
In particular, the models provide the mapping from time to Cartesian x, y, and
z coordinates or velocities. Also a representation of the usual radial velocity curve
is available. Inevitably, the models are somewhat redundant.
The calculations are based on algorithms implemented
:ref:`in PyAstronomy's Astrolib <keplerorbitpyasl>`.

.. contents:: :local:

KeplerEllipseModel: Location or velocity
----------------------------------------------

The KeplerEllipseModel model allows to model positions and velocities along
a Keplerian orbit.
On creation of the model, one can choose which of the coordinate
axes are of interest to the user. The coordinates or velocities will then
be returned as a one-dimensional array as explained the `evaluate` method of
the class.

Example: Evaluating and fitting model using MCMC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
	plt.errorbar(data[0::2], data[1::2], xerr=np.ones(20)*0.5,
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
	
	lims = {"a": [5., 10.], "per": [10., 15.], "e": [0., 1.], "tau": [0.5, 1.],
	        "Omega": [0., 360.], "w": [-5., 5.], "i": [90., 95.]}
	
	# Generate functions serving as uniform priors with lower and upper limit
	def getprior(ll, ul):
	    """ lower (ll) and upper (ul) limit """
	    def prior(pardict, p):
	        if pardict[p] < ll:
	            return -np.inf
	        elif pardict[p] > ul:
	            return -np.inf
	        else:
	            return 0
	    return prior
	
	# Define priors for parameters according to lims
	priors = {par:getprior(l[0], l[1]) for par, l in lims.items()}
	
	kem.fitEMCEE(time, data, yerr=np.ones(len(data))*0.5, dbfile="kemExample.emcee", priors=priors)
	
	# Plot the lowest deviance model
	ldmodel = kem.evaluate(np.linspace(0, kem["per"], 200))
	plt.plot(ldmodel[0::2], ldmodel[1::2], 'r--')
	
	plt.show()


KeplerEllipseModel: API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
.. autoclass:: KeplerEllipseModel
   :members:
   

Radial velocity curves: KeplerRVModel
-----------------------------------------

Radial velocity curves are among the most important tools to study planetary and stellar
systems. The KeplerRVModel provides an implementation of a radial velocity curve model
for one or more planets (or stars). Note, however, that mutual interactions between the
constituents of the system are not taken into account.

Example: Best-fit RV model and error analysis using KeplerRVModel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example demonstrates how to fit a radial velocity curve and posterior-based error estimation
(and parameter estimation) using emcee.

::
    
    from __future__ import print_function
    
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy.modelSuite import KeplerRVModel
    from PyAstronomy import funcFit as fuf
    
    # Generate artificial data ...
    jd = np.arange(100)
    rv = 1.5 * np.sin(jd / 37.6 * 2.*np.pi)
    # ... with some error
    rverr = 0.5
    rv += np.random.normal(0, rverr, len(jd))
    rverr = np.ones(len(rv)) * rverr
    
    # Get RV model with one planet (mp) and a potential constant offset
    # in RV (deg = 0)
    krvm = KeplerRVModel(mp=1, deg=0)
    # To obtain some useful estimate of the minimum mass of the companion,
    # we must specify the mass of the star (in terms of solar masses)
    krvm["mstar"] = 0.5
    
    # Let us have a look at the available parameters.
    # Note that not all are meant for fitting in this model (MA and a)!
    # There is also not much use in fitting 'mstar'. It may, however, be
    # used in combination with a prior to take into account its uncertainty in
    # the estimates.
    krvm.parameterSummary(sorting="ps")
    
    # We specify some guess parameters.
    krvm["per1"] = 37.0
    krvm["K1"] = 1.0
    krvm["e1"] = 0.0
    krvm["tau1"] = 17.0
    krvm["w1"] = 180.
    
    # Let us fit all of these but period ...
    krvm.thaw(["K1", "tau1", "w1", "e1", "c0"])
    # ... and now also the period
    krvm.thaw(["per1"])
    krvm.fit(jd, rv, yerr=rverr)
    # and then get the best-fit model
    kmo = krvm.evaluate(jd)
    
    # What about chi-square and RMS?
    chi = np.sum((rv - krvm.model)**2 / rverr**2)
    # Reduced chi-square
    rchi = chi / (len(rv) - len(krvm.freeParameters()))
    print("chi-square and reduced chi-square: %6.3f, %6.3f" % (chi, rchi))
    rms = np.std(rv - krvm.model)
    print("RMS: ", rms)
    
    plt.title("RV data (blue) and model (red)")
    plt.errorbar(jd, rv, yerr=rverr, fmt='b+')
    plt.plot(jd, krvm.model, 'r-')
    plt.show()
    
    
    # Now let us do some posterior-based error analysis using MCMC
    
    # Say, we want 20 burn-in iterations and, thereafter,
    # 50 further iterations (per walker).
    sampleArgs = {"iters": 50, "burn": 100}
    
    # Specify a bounded uniform prior on the eccentricity. Note that restrictions are not
    # automatically converted into priors (they may not ne uniform). Potentially further prior,
    # e.g., on per1 may be required to prevent wandering into 'forbidden territory'.
    priors = {"e1": fuf.FuFPrior("limuniform", upper=1, lower=0)}
    
    # Start the sampling (ps could be used to continue the sampling)
    ps = krvm.fitEMCEE(jd, rv, yerr=rverr, sampleArgs=sampleArgs, scales={"e": 0.05}, dbfile="chain1.emcee",
                       priors=priors)
    
    
    # Have a look at the posterior
    ta = fuf.TraceAnalysis("chain1.emcee")
    
    # What about the deviance (-2 log(Likelihood))
    ta.plotTraceHist("deviance")
    ta.show()
    
    # Expectation value and highest probability density interval for eccentricity
    ta.plotTraceHist("e1")
    print("Expectation value for eccentricity: ", ta.mean("e1"))
    print("90% HPD for eccentricity: ", ta.hpd("e1", cred=0.9))
    ta.show()
    ta.show()




KeplerRVModel: API
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: KeplerRVModel
   :members:
