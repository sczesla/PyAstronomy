from __future__ import print_function, division
import unittest
import os

class ModSuiteSanity(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    # Clean up example output from KeplerEllipseModel example
    if os.path.isfile("kemExample.tmp"):
      os.remove("kemExample.tmp")
  
  def sanity_rmcl_model(self):
    # Import some unrelated modules
    from numpy import arange, pi
    import matplotlib.pylab as plt
    # ... and the model suite
    from PyAstronomy import modelSuite as ms
    
    # Create Rossiter-McLaughlin object
    rmcl = ms.RmcL()
    # Set parameters
    rmcl.assignValue({"a":6.7, "lambda":7.2/180.0*pi, "epsilon":0.5, \
                      "P":1.74, "T0":0.2, "i":87.8/180.*pi, \
                      "Is":90.0/180.0*pi, "Omega":1.609e-5, "gamma":0.2})
    # Choose some time axis and calculate model
    time = arange(100)/100.0 * 0.2 + 0.1
    rv = rmcl.evaluate(time)
    
    # Let's see what happened...
    plt.ylabel("Radial velocity [stellar-radii/s]")
    plt.xlabel("Time [d]")
    plt.plot(time, rv, '.')
#     plt.show()
  
  def sanity_rmcl_fit(self):
    # Import some unrelated modules
    from numpy import arange, pi, random
    import matplotlib.pylab as plt
    # ... and the model suite
    from PyAstronomy import modelSuite as ms
    
    # Create Rossiter-McLaughlin object
    rmcl = ms.RmcL()
    # Set parameters
    rmcl.assignValue({"a":6.7, "lambda":7.2/180.0*pi, "epsilon":0.5, \
                      "P":1.74, "T0":0.2, "i":87.8/180.*pi, \
                      "Is":90.0/180.0*pi, "Omega":1.609e-5, "gamma":0.2})
    # Choose some time axis and calculate model
    time = arange(100)/100.0 * 0.2 + 0.1
    rv = rmcl.evaluate(time)
    
    # Add some noise.
    rv += random.normal(0.0, 0.05*rv.max(), rv.size)
    
    # Assign guess parameters
    rmcl.assignValue({"a":6.0, "lambda":7.2/180.0*pi, "epsilon":0.5, \
                      "P":1.74, "T0":0.17, "i":87.8/180.*pi, \
                      "Is":90.0/180.0*pi, "Omega":1.609e-5, "gamma":0.2})
    
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
#     plt.show()

  def sanity_SinRadVel(self):
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
#     plt.show()

  def sanity_KeplerEllipseModel(self):
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
    
#     plt.show()


#  def sanity_atanProfile(self):
#    from PyAstronomy import modelSuite as ms
#    import numpy as np
#    import matplotlib.pylab as plt
#    
#    # Create an instance of the AtanProfile ...
#    ap = ms.AtanProfile()
#    # ... and define some starting values
#    ap["A"] = 1.0
#    ap["mu"] = 5.0
#    ap["scale"] = 0.4
#    ap["sig"] = 5.0
#    
#    # Plot profile on given x-axis
#    x = np.linspace(-5,15,100)
#    plt.plot(x, ap.evaluate(x), 'b.-')
#    
#    # Determine the locations of the inflection
#    # points
#    print "Inflection points: ", ap.inflectionPoints()
#    
#    # Create instance of damped profile and copy
#    # the values from the first profile
#    apd = ms.AtanProfileDamped()
#    for p, v in ap.parameters().iteritems():
#      apd[p] = v
#    
#    # Specify the additional damping parameter
#    apd["tau"] = 2.0
#    # and plot
#    plt.plot(x, apd.evaluate(x), 'r.-')
##    plt.show()
  
  def sanity_lineListGaussModel(self):
    """
      Checking example of line list Gauss model
    """
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
#    plt.show()

  def sanity_VoigtAstroPExample(self):
    """
      Sanity of VoigtAstroP example
    """
    from PyAstronomy import modelSuite as ms
    import numpy as np
    import matplotlib.pylab as plt
    
    # Obtain an object of type VoigtAstroP ...
    v = ms.VoigtAstroP()
    # ... and set some parameters
    v["b"] = 87.7
    v["f"] = 0.5
    v["w0"] = 1214.0
    # Damping constant [cm]
    v["gamma"] = 2e-9
    
    # Generate wavelength axis ...
    wvl = np.linspace(1212.,1216.,200)
    # ... and evaluate model
    m = v.evaluate(wvl)
    
    # Plot result
    plt.plot(wvl, m, 'b.-')
#    plt.show()

  def sanity_LyATransmission(self):
    """
      Checking sanity of LyATransmission example
    """
    from PyAstronomy import modelSuite as ms
    import numpy as np
    import matplotlib.pylab as plt
    
    la = ms.LyaTransmission()
    # Set some parameters
    la["N"] = 5e17
    la["b"] = 12.2
    la["Dfrac"] = 1.9e-5
    
    # Set up wavelength axis ...
    wvl = np.linspace(1214.,1217.,1000)
    # ... and evaluate model
    m = la.evaluate(wvl)
    
    # Plot the result
    plt.plot(wvl, m, 'b.-')
#    plt.show()

  def sanity_RotBroadProfileExample(self):
    """
      Example of rotational broadening.
    """
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import modelSuite as ms
    
    # Get an instance of the model ...
    x = ms.RotBroadProfile()
    # ... and define some starting value
    x["xmax"] = 60.0
    x["A"] = 1.0
    x["eps"] = 0.8
    x["off"] = 0.0
    
    # Define a radial velocity axis
    vv = np.linspace(-90.,90.,200)
    
    # Construct some "data" and ...
    data = x.evaluate(vv)
    # ... add noise
    data += np.random.normal(0.0, 1e-3, data.size)
    
    # Fit the model using A, xmax, and eps as free
    # parameters ...
    x.thaw(["A", "xmax", "eps"])
    x.fit(vv, data)
    # ... and show the resulting parameter values.
    x.parameterSummary()
    
    # Plot the data and the model
    plt.plot(vv, data, 'bp')
    plt.plot(vv, x.model, 'r--')
#    plt.show()
  
  def sanity_RotBroadProfile(self):
    """
      Checking RotBroadProfile
    """
    import numpy as np
    from PyAstronomy import modelSuite as ms
    import scipy.integrate as sci
    
    # Get an instance of the model ...
    x = ms.RotBroadProfile()
    vv = np.linspace(-90.,90.,200)
    
    for i in range(10):
      # ... and define some starting value
      x["xmax"] = np.random.random()*50.0 + 30.0
      x["A"] = np.random.random()*10.0 + 1.0
      x["eps"] = np.random.random()

      d = x.evaluate(vv)
      a = sci.trapz(d, vv)
      self.assertAlmostEqual(x["A"], a, delta=1.0/200., msg="Incorrect profile normalization (" + \
                             "%g vs %g)" % (x["A"], a))
  
    x["eps"] = 0.0
    x["xmax"] = 50.0
    x["A"] = 1.0
    vv = np.linspace(-x["xmax"], x["xmax"], 100)
    d = x.evaluate(vv)
    
    y = d - 2.0/(np.pi*x["xmax"])*np.sqrt(1.0 - (vv/x["xmax"])**2)
    self.assertFalse(np.any(np.abs(y) > 1e-6), msg="Incorrect profile for eps=0.0")