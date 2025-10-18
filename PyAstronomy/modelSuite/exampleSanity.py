from __future__ import print_function, division
import unittest
import os


class TestModSuiteSanity(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        # Clean up example output from KeplerEllipseModel example
        if os.path.isfile("kemExample.emcee"):
            os.remove("kemExample.emcee")

    def testsanity_rmcl_model(self):
        """ Checking sanity of RmcL calculation (example) """
        # Import some unrelated modules
        from numpy import arange, pi
        import matplotlib.pylab as plt
        # ... and the model suite
        from PyAstronomy import modelSuite as ms

        # Create Rossiter-McLaughlin object
        rmcl = ms.RmcL()
        # Set parameters
        rmcl.assignValue({"a": 6.7, "lambda": 7.2 / 180.0 * pi, "epsilon": 0.5,
                          "P": 1.74, "T0": 0.2, "i": 87.8 / 180. * pi,
                          "Is": 90.0 / 180.0 * pi, "Omega": 1.609e-5, "gamma": 0.2})
        # Choose some time axis and calculate model
        time = arange(100) / 100.0 * 0.2 + 0.1
        rv = rmcl.evaluate(time)

        # Let's see what happened...
        plt.ylabel("Radial velocity [stellar-radii/s]")
        plt.xlabel("Time [d]")
        plt.plot(time, rv, '.')
#     plt.show()

    def testsanity_rmcl_fit(self):
        """ Checking sanity of RmcL fit (example) """
        # Import some unrelated modules
        from numpy import arange, pi, random
        import matplotlib.pylab as plt
        # ... and the model suite
        from PyAstronomy import modelSuite as ms

        # Create Rossiter-McLaughlin object
        rmcl = ms.RmcL()
        # Set parameters
        rmcl.assignValue({"a": 6.7, "lambda": 7.2 / 180.0 * pi, "epsilon": 0.5,
                          "P": 1.74, "T0": 0.2, "i": 87.8 / 180. * pi,
                          "Is": 90.0 / 180.0 * pi, "Omega": 1.609e-5, "gamma": 0.2})
        # Choose some time axis and calculate model
        time = arange(100) / 100.0 * 0.2 + 0.1
        rv = rmcl.evaluate(time)

        # Add some noise.
        rv += random.normal(0.0, 0.05 * rv.max(), rv.size)

        # Assign guess parameters
        rmcl.assignValue({"a": 6.0, "lambda": 7.2 / 180.0 * pi, "epsilon": 0.5,
                          "P": 1.74, "T0": 0.17, "i": 87.8 / 180. * pi,
                          "Is": 90.0 / 180.0 * pi, "Omega": 1.609e-5, "gamma": 0.2})

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

    def testsanity_rmclell_calc(self):
        """
        Checking sanity of RmcLell (example)
        """
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
        rmcl.assignValue({"a":6.7, "lambda":7.2/180.0*pi, "epsilon":0.5, \
                            "P":1.74, "T0":0.2, "i":87.8/180.*pi, \
                            "Is":90.0/180.0*pi, "Omega":1.609e-5, "gamma":0.2})
        rmel.assignValue({"a":6.7, "lambda":7.2/180.0*pi, "epsilon":0.5, \
                            "P":1.74, "tau":0.2, "i":87.8/180.*pi, "w":-90/180.*pi, \
                            "e":0.05, "Is":90.0/180.0*pi, "Omega":1.609e-5, "gamma":0.2})
        
        
        # Choose some time axis and calculate model
        time = arange(100)/100.0 * 0.2 + 0.1
        rvc = rmcl.evaluate(time)
        rve = rmel.evaluate(time)
        
        # Let's see what happened...
#         plt.ylabel("Radial velocity [stellar-radii/s]")
#         plt.xlabel("Time [d]")
#         plt.plot(time, rvc, 'b.-', label="circular")
#         plt.plot(time, rve, 'r.-', label="elliptical")
#         plt.legend()
#         plt.show()

    def testsanity_rmcl_vs_rmclell(self):
        """ Cross-checking Rmcl and RmcLell """
        from numpy import arange, pi
        import numpy as np
        # ... and the model suite
        from PyAstronomy import modelSuite as ms
        
        # Create Rossiter-McLaughlin object
        rmcl = ms.RmcL()
        r2 = ms.RmcLell()
        
        np.random.seed(9234667)
        
        for i in range(10):
        
            a = np.random.random()*5 + 3
            l = np.random.random()*180 - 90
            inc = np.random.random()*2 + 88
            Omega = 1e-5 + np.random.random()*1e-5
        
            # Set parameters
            rmcl.assignValue({"a":a, "lambda":l/180.0*pi, "epsilon":0.5, \
                            "P":1.74, "T0":0.2, "i":inc/180.*pi, \
                            "Is":80.0/180.0*pi, "Omega":Omega, "gamma":0.2})
        
            # Set parameters
            r2.assignValue({"a":a, "lambda":l/180.0*pi, "epsilon":0.5, \
                            "P":1.74, "tau":0.2, "i":inc/180.*pi, \
                            "Is":80.0/180.0*pi, "Omega":Omega, "gamma":0.2,
                            "e":0.0, "w":-90/180*pi})
        
            # Choose some time axis and calculate model
            time = arange(20)/20.0 * 0.2 - 0.1 + rmcl["T0"]
            rv = rmcl.evaluate(time)
            rv2 = r2.evaluate(time)
            
            d = np.max(np.abs(rv-rv2))
            m = np.max(np.abs(rv))
            self.assertAlmostEqual(d/m, 0.0, delta=1e-8, msg="Elliptical and circular orbit solution for RmcL and RmcLell do not match. " + \
                                   str(r2.parameters()))

    def testsanity_SinRadVel(self):
        # Import some unrelated modules
        from numpy import arange, random, ones
        import matplotlib.pylab as plt
        # ... and now the radVel module
        from PyAstronomy.modelSuite import radVel as rv

        # Create Radial Velocity object
        r = rv.SinRadVel()
        # Set parameters
        r.assignValue({"P": 1.8, "T0": 0.25, "K": 0.5, "rv0": 10.0})
        # Choose some time axis and calculate model
        time = arange(100) / 100.0 * 3.0 - 1.5
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
        r.fit(time, rvData, yerr=ones(y.size) * 0.05)

        # Show fit results
        print("Fitted values:")
        r.parameterSummary()

        # Let's see what happened...
        plt.ylabel("Radial velocity [km/s]")
        plt.xlabel("Radial velocity [d]")
        plt.errorbar(time, rvData, yerr=ones(y.size) * 0.05, fmt='b.')
        plt.plot(time, y, 'r-')
#     plt.show()

    def testsanity_KeplerEllipseModel(self):
        """ Sanity of Kepler Ellipse Model example """
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
        
        #plt.show()



#  def testsanity_atanProfile(self):
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
# plt.show()

    def testsanity_lineListGaussModel(self):
        """
        Checking example of line list Gauss model
        """
        from PyAstronomy import modelSuite as ms
        import numpy as np
        import matplotlib.pylab as plt

        # Create our line list with 4 line
        lineList = np.zeros((4, 3))
        # Assign wavelengths (in A)
        lineList[0, 0] = 5002.37
        lineList[1, 0] = 5005.9
        lineList[2, 0] = 5007.52
        lineList[3, 0] = 5007.64
        # Assign EWs (in A)
        lineList[0, 1] = 0.01
        lineList[1, 1] = 0.05
        lineList[2, 1] = 0.009
        lineList[3, 1] = 0.12
        # Assign depths (0-1)
        lineList[0, 2] = 0.97
        lineList[1, 2] = 0.9
        lineList[2, 2] = 0.99
        lineList[3, 2] = 0.35

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
        plt.subplot(2, 1, 1)
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
            llg["A" + str(i + 1)] += np.random.normal(0.0, 0.1)
        # Use all line strengths for fitting
        llg.thawLineStrengths()
        # and fit
        llg.fit(wvl, data)
        # Plot the result
        plt.subplot(2, 1, 2)
        plt.errorbar(wvl, data, yerr=np.ones(len(wvl)) * 0.01, fmt='bp')
        plt.plot(wvl, llg.evaluate(wvl), 'r--')
#    plt.show()

    def testsanity_VoigtAstroPExample(self):
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
        wvl = np.linspace(1212., 1216., 200)
        # ... and evaluate model
        m = v.evaluate(wvl)

        # Plot result
        plt.plot(wvl, m, 'b.-')
#    plt.show()


    def testsanity_VoigtAstroP_R_Example(self):
        """
        Sanity of VoigtAstroP example with instrumental resolution
        """

        from PyAstronomy import modelSuite as ms
        import numpy as np
        import matplotlib.pylab as plt
        
        # Obtain an object of type VoigtAstroP ...
        v = ms.VoigtAstroP()
        # ... and set some parameters
        v["b"] = 40.7
        v["f"] = 0.5
        v["w0"] = 1214.0
        # Damping constant [cm]
        v["gamma"] = 2e-9
        
        # Generate wavelength axis ...
        wvl = np.linspace(1212.,1216.,200)
        # ... and evaluate model
        m = v.evaluate(wvl)
        
        # Add (Gaussian) instrumental broadening with resolution 5000
        v["R"] = 5000
        mr = v.evaluate(wvl)
        
        # Plot result
#         plt.plot(wvl, m, 'b.-', label="R = inf")
#         plt.plot(wvl, mr, 'r.-', label="R = 5000")
#         plt.legend()
#         plt.show()

    def testsanity_LyATransmission(self):
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
        wvl = np.linspace(1214., 1217., 1000)
        # ... and evaluate model
        m = la.evaluate(wvl)

        # Plot the result
        plt.plot(wvl, m, 'b.-')
#    plt.show()

    def testsanity_RotBroadProfileExample(self):
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
        vv = np.linspace(-90., 90., 200)

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

    def testsanity_RotBroadProfile(self):
        """
        Checking RotBroadProfile
        """
        import numpy as np
        from PyAstronomy import modelSuite as ms
        import scipy.integrate as sci

        # Get an instance of the model ...
        x = ms.RotBroadProfile()
        vv = np.linspace(-90., 90., 200)

        for i in range(10):
            # ... and define some starting value
            x["xmax"] = np.random.random() * 50.0 + 30.0
            x["A"] = np.random.random() * 10.0 + 1.0
            x["eps"] = np.random.random()

            d = x.evaluate(vv)
            a = sci.trapezoid(d, vv)
            self.assertAlmostEqual(x["A"], a, delta=1.0 / 200., msg="Incorrect profile normalization (" +
                                   "%g vs %g)" % (x["A"], a))

        x["eps"] = 0.0
        x["xmax"] = 50.0
        x["A"] = 1.0
        vv = np.linspace(-x["xmax"], x["xmax"], 100)
        d = x.evaluate(vv)

        y = d - 2.0 / (np.pi * x["xmax"]) * np.sqrt(1.0 - (vv / x["xmax"])**2)
        self.assertFalse(np.any(np.abs(y) > 1e-6),
                         msg="Incorrect profile for eps=0.0")


    def testsanity_KeplerRVModel_example(self):
        """
        Checking sanity of KeplerRVModel example
        """
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
        chi = np.sum( (rv - krvm.model)**2 / rverr**2 )
        # Reduced chi-square
        rchi = chi / (len(rv) - len(krvm.freeParameters()))
        print("chi-square and reduced chi-square: %6.3f, %6.3f" % (chi, rchi))
        rms = np.std(rv - krvm.model)
        print("RMS: ", rms)
        
        plt.title("RV data (blue) and model (red)")
        plt.errorbar(jd, rv, yerr=rverr, fmt='b+')
        plt.plot(jd, krvm.model, 'r-')
        #=======================================================================
        # plt.show()
        #=======================================================================
        
        
        
        # Now let us do some posterior-based error analysis using MCMC
        
        # Say, we want 20 burn-in iterations and, thereafter,
        # 50 further iterations (per walker).
        sampleArgs = {"iters":50, "burn":100}
        
        # Specify a bounded uniform prior on the eccentricity. Note that restrictions are not
        # automatically converted into priors (they may not ne uniform). Potentially further prior,
        # e.g., on per1 may be required to prevent wandering into 'forbidden territory'.
        priors = {"e1":fuf.FuFPrior("limuniform", upper=1, lower=0)}
        
        # Start the sampling (ps could be used to continue the sampling)
        ps = krvm.fitEMCEE(jd, rv, yerr=rverr, sampleArgs=sampleArgs, scales={"e":0.05}, dbfile="chain1.emcee", \
            priors=priors)
        
        
        # Have a look at the posterior
        ta = fuf.TraceAnalysis("chain1.emcee")
        
        # What about the deviance (-2 log(Likelihood))
        ta.plotTraceHist("deviance")
        #=======================================================================
        # ta.show()
        #=======================================================================
        
        # Expectation value and highest probability density interval for eccentricity
        ta.plotTraceHist("e1")
        print("Expectation value for eccentricity: ", ta.mean("e1"))
        print("90% HPD for eccentricity: ", ta.hpd("e1", cred=0.9))
        #=======================================================================
        # ta.show()
        #=======================================================================
        


class TestVoigtAstroPSanity(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def testsanity_normalization(self):
        """
        Normalization of AstroVoigtP
        """
        import numpy as np
        from PyAstronomy import modelSuite as ms
        from PyAstronomy import pyaC
        
        v = ms.VoigtAstroP()
        
        # Hypothetical wvls [A]
        w0s = [1000., 5000, 10000.]
        
        # (pi e**2)/(m_e c)
        const = (4.803e-10)**2*np.pi / (9.11e-28*29979245800.0)
        
        fs = [1, 100]
        
        for f in fs:
            for w0 in w0s:
        
                v["w0"] = w0
                v["b"] = 100.
                v["gamma"] = 1e-10
                v["f"] = f
        
                dw = 20.0
                w = np.linspace(w0-dw, w0+dw, 1000)
        
                m = v.evaluate(w)
        
                i = pyaC.ibtrapz(w/1e8, m*29979245800.0/(w/1e8)**2 , (w0-dw)/1e8, (w0+dw)/1e8)
                
                self.assertAlmostEqual(i/const, f, delta=1e-2, msg="Normalization of AstroVoigtP is broken: f, w0, i: % g, % g, % g" % (f, w0, i))

    def testsanity_instrumentalResolution(self):
        """
        Checking integrity of instrumental resolution in VoigtAstroP
        """
        import numpy as np
        from PyAstronomy import modelSuite as ms
        from PyAstronomy import pyasl
        
        v = ms.VoigtAstroP()
        
        w0 = 10830
        
        for R in [2500, 5000, 80000]:
        
            v["w0"] = w0
            v["b"] = 10.
            v["gamma"] = 1e-8
            v["f"] = 100.0
            v["R"] = 0
        
            dw = 40.0
            w = np.linspace(w0-dw, w0+dw, 4000)
        
            m = v.evaluate(w)
        
            v["R"] = R
            m2 = v.evaluate(w)
        
            
            fb = pyasl.instrBroadGaussFast(w, m, R, edgeHandling=None, fullout=False, maxsig=None)

            d = 1000
            
            self.assertAlmostEqual(np.max(np.abs(fb[d:-d] - m2[d:-d])/m2[d:-d]), 0.0, delta=1e-10,
                                   msg="VoigtAstroP instrumental broadening broken for R = " + str(R))

            
