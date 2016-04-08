from __future__ import print_function, division
import unittest
import os

class ExampleSanity(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass

  def sanity_firstExample(self):
    # Import numpy and matplotlib
    from numpy import arange, sqrt, exp, pi, random, ones
    import matplotlib.pylab as plt
    # ... and now the funcFit package
    from PyAstronomy import funcFit as fuf
    
    # Before we can start fitting, we need something to fit.
    # So let us create some data...
    
    # Creating a Gaussian with some noise
    # Choose some parameters...
    gPar = {"A":-5.0, "sig":10.0, "mu":10.0, "off":1.0, "lin":0.0}
    # Calculate profile
    x = arange(100) - 50.0
    y = gPar["off"] + gPar["A"] / sqrt(2*pi*gPar["sig"]**2) \
        * exp(-(x-gPar["mu"])**2/(2*gPar["sig"]**2))
    # Add some noise
    y += random.normal(0.0, 0.01, x.size)
    # Let us see what we have done...
    plt.plot(x, y, 'bp')
    
    # Now we can start exploiting the funcFit functionality to
    # fit a Gaussian to our data. In the following lines, we
    # create a fitting object representing a Gaussian and set guess parameters.
    
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
    print("Parameters and guess values: ")
    print("  A   : ", gf["A"])
    print("  sig : ", gf["sig"])
    print("  off : ", gf["off"])
    print("  mu  : ", gf["mu"])
    print("")
    
    # Now some of the strengths of funcFit are demonstrated; namely, the
    # ability to consider some parameters as free and others as fixed.
    # By default, all parameters of the GaussFit1d are frozen.
    
    # Show values and names of frozen parameters
    print("Names and values of FROZEN parameters: ", gf.frozenParameters())
    
    # Which parameters shall be variable during the fit?
    # 'Thaw' those (the order is irrelevant)
    gf.thaw(["A", "sig", "off", "mu"])
    
    # Let us assume that we know that the amplitude is negative, i.e.,
    # no lower boundary (None) and 0.0 as upper limit.
    gf.setRestriction({"A":[None,0.0]})
    
    # Now start the fit
    gf.fit(x, y, yerr=ones(x.size)*0.01)
    
    # Write the result to the screen and plot the best fit model
    gf.parameterSummary()
    plt.plot(x, gf.model, 'r--')
    
    # Show the data and the best fit model
#     plt.show()


  def sanity_CustomModel(self):
    # Import numpy and matplotlib
    from numpy import arange, random
    import matplotlib.pylab as plt
    # ... and now the funcFit package
    from PyAstronomy import funcFit as fuf
    
    class StraightLine(fuf.OneDFit):
      """
        Implements a straight line of the form y = "off" + x * "lin".
      """
    
      def __init__(self):
        fuf.OneDFit.__init__(self, ["off", "lin"])
    
      def evaluate(self, x):
        """
          Calculates and returns model according to the \
          current parameter values.
    
          Parameters:
            - `x` - Array specifying the positions at \
                    which to evaluate the model.
        """
        y = self["off"] + (self["lin"] * x)
        return y
      
    # Generate some data and add noise
    x = arange(100)
    y = 10.0 + 2.0 * x + random.normal(0.0, 5.0, 100)
    
    # Create fitting class instance and set initial guess
    # Note that all parameters are frozen by default
    lf = StraightLine()
    lf["off"] = 20.0
    lf["lin"] = 1.0
    # Thaw parameters
    lf.thaw(["off", "lin"])
    
    # Start fitting
    lf.fit(x, y)
    
    # Investigate the result
    lf.parameterSummary()
    plt.plot(x, y, 'bp')
    plt.plot(x, lf.model, 'r--')
#     plt.show()
    
  def sanity_Relations(self):
    # import numpy and matplotlib
    from numpy import arange, random
    import matplotlib.pylab as plt
    # ... and now the funcFit package
    from PyAstronomy import funcFit as fuf
    
    class StraightLine(fuf.OneDFit):
      """
        Implements a straight line of the form y = "off" + x * "lin".
      """
    
      def __init__(self):
        fuf.OneDFit.__init__(self, ["off", "lin"])
    
      def evaluate(self, x):
        """
          Calculates and returns model according to the current parameter values.
    
          Parameters:
            - x - Array specifying the positions at which to evaluate the model.
        """
        y = self["off"] + (self["lin"] * x)
        return y
    
    
    # Create a function, which defines the relation.
    
    def getLinearRelation(factor):
      def linOffRel(off):
        """
          Function used to relate parameters "lin" and "off".
        """
        return factor * off
      return linOffRel
    
    # Note, above we used a nested function (a closure) to define
    # the relation. This approach is very flexible. If we were already
    # sure about the value of ``factor'' (e.g., 10.0), we could
    # simply have used:
    #
    # def linOffRel(off):
    #   return 10.0 * off
    
    # Generate some data with noise
    x = arange(100)
    y = 100.0 + 2.0 * x + random.normal(0.0, 5.0, 100)
    
    # Create fitting class instance and set initial guess
    lf = StraightLine()
    lf["off"] = 20.0
    lf["lin"] = 1.0
    # Thaw parameters
    lf.thaw(["off", "lin"])
    
    # Assume we know about a relation between 'lin' and 'off'
    # In particular, lin = 9.0 * off. We use the function getLinearRelation
    # to obtain a function object defining the relation.
    lf.relate("lin", ["off"], getLinearRelation(9))
    
    # Start fitting
    lf.fit(x, y)
    
    # Investigate the result
    lf.parameterSummary()
    plt.plot(x, y, 'bp')
    plt.plot(x, lf.model, 'r--')
#     plt.show()
  
  def sanity_CombiningModels(self):
    # Import numpy and matplotlib
    from numpy import arange, sqrt, exp, pi, random, ones
    import matplotlib.pylab as plt
    # ... and now the funcFit package
    from PyAstronomy import funcFit as fuf
    
    # Creating Gaussians with some noise
    # Choose some parameters...
    gPar1 = {"A":-5.0, "sig":10.0, "mu":20.0, "off":1.0, "lin":0.0}
    gPar2 = {"A":+10.0, "sig":10.0, "mu":-20.0, "off":0.0, "lin":0.0}
    # Calculate profile
    x = arange(100) - 50.0
    y = gPar1["off"] + gPar1["A"] / sqrt(2*pi*gPar1["sig"]**2) \
        * exp(-(x-gPar1["mu"])**2/(2*gPar1["sig"]**2))
    y -= gPar2["off"] + gPar2["A"] / sqrt(2*pi*gPar2["sig"]**2) \
        * exp(-(x-gPar2["mu"])**2/(2*gPar2["sig"]**2))
    # Add some noise
    y += random.normal(0.0, 0.01, x.size)
    # Let us see what we have done...
    plt.plot(x, y, 'bp')
    
    # Now let us come to the fitting
    # First, we create two Gauss1d fit objects
    gf1 = fuf.GaussFit1d()
    gf2 = fuf.GaussFit1d()
    
    # Assign guess values for the parameters
    gf1["A"] = -0.3
    gf1["sig"] = 3.0
    gf1["off"] = 0.0
    gf1["mu"] = +5.0
    
    gf2["A"] = 3.0
    gf2["sig"] = 15.0
    gf2["off"] = 1.0
    gf2["mu"] = -10.0
    
    # Which parameters shall be variable during the fit?
    # 'Thaw' those (the order is irrelevant)
    gf1.thaw(["A", "sig", "mu"])
    gf2.thaw(["sig", "mu", "off"])
    
    # Our actual model is the sum of both Gaussians
    twoG = gf1 + gf2
    
    # Show a description of the model depending on the
    # names of the individual components
    print()
    print("Description of the model: ", twoG.description())
    print()
    
    # Note that now the parameter names changed!
    # Each parameter is now named using the "property"
    # (e.g., 'A' or 'sig') as the first part, the component
    # "root name" (in this case 'Gaussian') and a component
    # number in parenthesis.
    print("New parameter names and values: ")
    twoG.parameterSummary()
    
    # We forgot to thaw the amplitude of the second Gaussian, but
    # we can still do it, but we have to refer to the correct name:
    # either by using the (new) variable name:
    twoG.thaw("A_Gaussian(2)")
    # or by specifying property name, root name, and component number
    # separately (note that a tuple is used to encapsulate them):
    twoG.thaw(("A", "Gaussian", 2))
    # We decide to rather freeze the offset of the second
    # Gaussian (we could have used a tuple here, too).
    twoG.freeze("off_Gaussian(2)")
    
    # Start fit as usual
    twoG.fit(x,y,yerr=ones(x.size)*0.01)
    
    # Write the result to the screen and plot the best fit model
    print()
    print("--------------------------------")
    print("Parameters for the combined fit:")
    print("--------------------------------")
    twoG.parameterSummary()
    
    # Show the data and the best fit model
    plt.plot(x, twoG.model, 'r--')
#     plt.show()

  def sanity_CustomObjectiveFunctions(self):
    # Import numpy and matplotlib
    from numpy import arange, exp, random, ones, sum, abs
    import matplotlib.pylab as plt
    # Import funcFit
    from PyAstronomy import funcFit as fuf
    
    # Define parameters of faked data
    A = 1.0
    tau = 10.
    off = 0.2
    t0 = 40.
    
    # Calculate fake data set
    x = arange(100)
    y = A*exp(-(x-t0)/tau) * (x>t0) + off
    y += random.normal(0., 0.1, 100)
    yerr = ones(100)*0.01
    
    # Exponential decay model
    edf = fuf.ExpDecayFit1d()
    
    # Define free quantities
    edf.thaw(["A", "tau", "off", "t0"])
    # Let the amplitude be positive
    edf.setRestriction({"A":[0.0,None]})
    # Define initial guess
    edf.assignValue({"A":1.0, "tau": 15., "off":0.2, "t0":50.})
    
    # Do not use chi square, but the linear deviation from model
    # to evaluate quality of fit.
    # Use the "MiniFunc" decorator to define your custom objective
    # function. This decorator takes the fitting object as an
    # argument. The function has to accept two arguments: the
    # fitting object and the list of free parameters.
    @fuf.MiniFunc(edf)
    def mini(edf, P):
      m = sum(abs(edf.model - edf.y)/edf.yerr)
      print("mini - current parameters: ", P, ", value is: ", m)
      return m
    
    # Carry out fit WITH SELF-DEFINED OBJECTIVE FUNCTION
    edf.fit(x, y, yerr=yerr, miniFunc=mini)
    
    # Show parameter values and plot best-fit model.
    edf.parameterSummary()
    plt.errorbar(x,y,yerr)
    plt.plot(x, edf.model, 'r-')
#     plt.show()

  def sanity_Overbinning(self):
    # Import numpy and matplotlib
    from numpy import arange, sqrt, exp, pi, random, ones
    import matplotlib.pylab as plt
    # ... and now the funcFit package
    from PyAstronomy import funcFit as fuf
    
    # Creating a Gaussian with some noise
    # Choose some parameters...
    gPar = {"A":-5.0, "sig":10.0, "mu":10.0, "off":1.0, "lin":0.0}
    # Calculate profile
    x = arange(20)/20.0 * 100.0 - 50.0
    y = gPar["off"] + gPar["A"] / sqrt(2*pi*gPar["sig"]**2) \
        * exp(-(x-gPar["mu"])**2/(2*gPar["sig"]**2))
    # Add some noise
    y += random.normal(0.0, 0.01, x.size)
    # Let us see what we have done...
    plt.plot(x, y, 'bp')
    
    # First, we create a "GaussFit1d_Rebin" class object (note that the
    # class object has still to be instantiated, the name is arbitrary).
    GaussFit1d_Rebin = fuf.turnIntoRebin(fuf.GaussFit1d)
    # Do the instantiation and specify how the overbinning should be
    # carried out.
    gf = GaussFit1d_Rebin()
    gf.setRebinArray_Ndt(x, 10, x[1]-x[0])
    # See what parameters are available
    print("List of available parameters: ", gf.availableParameters())
    # Set guess values for the parameters
    gf["A"] = -10.0
    gf["sig"] = 15.77
    gf["off"] = 0.87
    gf["mu"] = 7.5
    # Let us see whether the assignment worked
    print("Parameters and guess values: ")
    print("  A   : ", gf["A"])
    print("  sig : ", gf["sig"])
    print("  off : ", gf["off"])
    print("  mu  : ", gf["mu"])
    print("")
    
    # Now some of the strengths of funcFit are demonstrated; namely, the
    # ability to consider some parameters as free and others as fixed.
    # By default, all parameters of the GaussFit1d are frozen.
    
    # Show values and names of frozen parameters
    print("Names and values if FROZEN parameters: ", gf.frozenParameters())
    
    # Which parameters shall be variable during the fit?
    # 'Thaw' those (the order is irrelevant)
    gf.thaw(["A", "sig", "off", "mu"])
    
    # Let us assume that we know that the amplitude is negative, i.e.,
    # no lower boundary (None) and 0.0 as upper limit.
    gf.setRestriction({"A":[None,0.0]})
    
    # Now start the fit
    gf.fit(x, y, yerr=ones(x.size)*0.01)
    
    # Write the result to the screen and plot the best fit model
    gf.parameterSummary()
    # Plot the final best-fit model
    plt.plot(x, gf.model, 'rp--')
    # Show the overbinned (=unbinned) model, indicate by color
    # which point are averaged to obtain a point in the binned
    # model.
    for k, v in gf.rebinIdent.items():
      c = "y"
      if k % 2 == 0: c = "k"
      plt.plot(gf.rebinTimes[v], gf.unbinnedModel[v], c+'.')
    
    # Show the data and the best fit model
#     plt.show()

  def sanity_simultaneousFit(self):
    from PyAstronomy import funcFit as fuf
    import numpy
    import matplotlib.pylab as plt
    
    # Set up two different x axes.
    x1 = numpy.arange(100.)/100. - 0.5
    x2 = numpy.arange(150.)/150. - 0.25
    
    # Getting the models ...
    gauss = fuf.GaussFit1d()
    calor = fuf.CauchyLorentz1d()
    # and assign parameters.
    gauss.assignValue({"A":0.02, "sig":0.1, "mu":0.0, "off":1.0, "lin":0.0})
    calor.assignValue({"A":0.07, "g":0.1, "mu":0.2, "off":1.0, "lin":0.0})
    
    # Create noisy data.
    y1 = gauss.evaluate(x1) + numpy.random.normal(0., 0.01, 100)
    y2 = calor.evaluate(x2) + numpy.random.normal(0., 0.01, 150)
    
    # Plot the noisy data.
    plt.subplot(2,1,1)
    plt.errorbar(x1, y1, yerr=numpy.ones(100)*0.01)
    plt.subplot(2,1,2)
    plt.errorbar(x2, y2, yerr=numpy.ones(150)*0.01)
    
    # Now, get ready two fit the data sets simultaneously.
    sf = fuf.SyncFitContainer()
    # Tell the class about the two components and save the
    # component numbers assigned to them:
    gaussCno = sf.addComponent(gauss)
    calorCno = sf.addComponent(calor)
    
    print("Component numbers in the syncFit container:")
    print("  Gauss: ", gaussCno, ",  Cauchy-Lorentz: ", calorCno)
    print()
    
    # See what happened to the parameters in the
    # simultaneous fitting class.
    # The variable names have changed.
    sf.parameterSummary()
    
    # Thaw all parameters (for later fit) ...
    sf.thaw(list(sf.parameters()))
    # but not the linear term.
    sf.freeze(["lin_Gaussian[s1]", "lin_CauLor[s2]"])
    
    # Tell the class about the identity of parameters,
    # either by using the "property name" of the parameter:
    sf.treatAsEqual("off")
    # or by specifying the names explicitly.
    sf.treatAsEqual(["g_CauLor[s2]", "sig_Gaussian[s1]"])
    
    # See what happened to the parameters in the
    # simultaneous fitting class.
    print()
    print("Parameters after 'treatAsEqual' has been applied:")
    sf.parameterSummary()
    
    # Randomize starting values.
    for fp in sf.freeParamNames():
      sf[fp] = sf[fp] + numpy.random.normal(0., 0.05)
    
    # Set up the data appropriately.
    data = {gaussCno:[x1, y1], calorCno:[x2, y2]}
    yerr = {gaussCno: numpy.ones(100)*0.01, \
            calorCno: numpy.ones(150)*0.01}
    
    # Start the fit.
    sf.fit(data, yerr=yerr)
    
    # Show the best-fit values.
    print()
    print("Best-fit parameters:")
    sf.parameterSummary()
    
    # Plot the best-fit model(s).
    plt.subplot(2,1,1)
    plt.plot(x1, sf.models[gaussCno], 'r--')
    plt.subplot(2,1,2)
    plt.plot(x2, sf.models[calorCno], 'r--')
    
#     plt.show()

  def sanity_2dCircularFit(self):
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit as fuf
    
    # Get the circular model and assign
    # parameter values
    c = fuf.Circle2d()
    c["r"] = 1.0
    c["t0"] = 0.0
    c["per"] = 3.0
    
    # Evaluate the model at a number of
    # time stamps
    t = np.linspace(0.0, 10.0, 20)
    pos = c.evaluate(t)
    
    # Add some error to the "measurement"
    pos += np.reshape(np.random.normal(0.0, 0.2, pos.size), pos.shape)
    err = np.reshape(np.ones(pos.size), pos.shape) * 0.2
    
    # Define free parameters and fit the model
    c.thaw(["r", "t0", "per"])
    c.fit(t, pos, yerr=err)
    c.parameterSummary()
    
    # Evaluate the model at a larger number of
    # points for plotting
    tt = np.linspace(0.0, 10.0, 200)
    model = c.evaluate(tt)
    
    # Plot the result
    plt.errorbar(pos[::,0], pos[::,1], yerr=err[::,1], \
                 xerr=err[::,0], fmt='bp')
    plt.plot(model[::,0], model[::,1], 'r--')
#    plt.show()

  def sanity_2dGaussFit(self):
    from PyAstronomy import funcFit as fuf
    import numpy as np
    import matplotlib.pylab as plt
    
    # Constructing the individual coordinate axes
    x = np.linspace(-2.,2.,50)
    y = np.linspace(-2.,2.,50)
    # Applying funcFit's "coordinateGrid" helper function
    # to built appropriate array-index -> coordinate mapping
    # needed for nD fitting.
    g = fuf.coordinateGrid(x, y)
    
    # Create the 2d-Gaussian model and assign
    # some model parameters.
    gf = fuf.GaussFit2d()
    gf["sigx"] = 0.75
    gf["sigy"] = 0.4
    gf["A"] = 1.0
    gf["rho"] = 0.4
    
    # Get the "data" by evaluating the model
    # and adding some noise. Note that the coordinate
    # mapping (array g) is passed to evaluate here.
    im = gf.evaluate(g)
    im += np.reshape(np.random.normal(0.0, 0.1, 2500), (50,50))
    err = np.ones((50,50))*0.1
    
    # Thaw parameters and fit
    gf.thaw(["A", "rho"])
    gf.fit(g, im, yerr=err)
    
    # Show the resulting parameter values ...
    gf.parameterSummary()
    
    # ... and plot the result.
    plt.title("Image data")
    plt.imshow(np.transpose(im), origin="lower")
#    plt.show()
    plt.title("Residuals")
    plt.imshow(np.transpose(im - gf.evaluate(g)), origin="lower")
#    plt.show()

  def sanity_2gGaussFitTupleExample(self):
    from PyAstronomy import funcFit as fuf
    import numpy as np
    import matplotlib.pylab as plt
    
    # Constructing the individual coordinate axes
    x = np.linspace(-2.,2.,50)
    y = np.linspace(-2.,2.,50)
    
    # Create the 2d-Gaussian model and assign
    # some model parameters.
    gf = fuf.GaussFit2dTuple()
    gf["sigx"] = 0.75
    gf["sigy"] = 0.4
    gf["A"] = 1.0
    gf["rho"] = 0.4
    
    # Get the "data" by evaluating the model
    # and adding some noise. Note that the coordinate
    # mapping (array g) is passed to evaluate here.
    im = gf.evaluate((x,y))
    im += np.reshape(np.random.normal(0.0, 0.1, 2500), (50,50))
    err = np.ones((50,50))*0.1
    
    # Thaw parameters and fit
    gf.thaw(["A", "rho"])
    gf.fit((x,y), im, yerr=err)
    
    # Show the resulting parameter values ...
    gf.parameterSummary()
    
    # ... and plot the result.
    plt.title("Image data")
    plt.imshow(np.transpose(im), origin="lower")
#    plt.show()
    plt.title("Residuals")
    plt.imshow(np.transpose(im - gf.evaluate((x,y))), origin="lower")
#    plt.show()

  def sanity_coordinateGridExample(self):
    from PyAstronomy import funcFit as fuf
    import numpy as np
    
    # Constructing the two individual coordinate axes
    x = np.linspace(-2.,2.,50)
    y = np.linspace(-2.,2.,50)
    
    # Applying funcFit's "coordinateGrid" helper function
    # to built appropriate array-index -> coordinate mapping
    # needed for nD fitting.
    g = fuf.coordinateGrid(x, y)
    
    print("(x, y) coordinates at index (11, 28): ", g[11,28])

  def sanity_CashStatisticsExample(self):
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit as fuf
    
    # Get a Gaussian fitting object and
    # set some parameters
    g = fuf.GaussFit1d()
    g["A"] = 5.1
    g["sig"] = 0.5
    g["mu"] = 3.94
    
    # Generate some data with Poisson statistics
    x = np.linspace(0.0, 7., 50)
    y = np.zeros(len(x))
    for i in range(len(x)):
      y[i] = np.random.poisson(g.evaluate(x[i]))
    
    # Choose free parameters and "disturb" the
    # starting parameters for the fit a little.
    g.thaw(["A", "sig", "mu"])
    for par in g.freeParamNames():
      g[par] += np.random.normal(0.0, g[par]*0.1)
    
    # Fit using Cash statistic and print out
    # result.
    g.fit(x, y, miniFunc="cash79")
    g.parameterSummary()
    
    # Plot the result
    plt.plot(x, y, 'bp')
    plt.plot(x, g.evaluate(x), 'r--')
#    plt.show()
  
  def sanity_steppar1(self):
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit as fuf
    
    # Set up a Gaussian model
    # and create some "data"
    x = np.linspace(0,2,100)
    gf = fuf.GaussFit1d()
    gf["A"] = 0.87
    gf["mu"] = 1.0
    gf["sig"] = 0.2
    y = gf.evaluate(x)
    y += np.random.normal(0.0, 0.1, len(x))
    
    # Thaw parameters, which are to be fitted. Note
    # that those parameters will also be fitted during
    # the stepping; no further parameters will be thawed.
    gf.thaw(["A", "mu", "sig"])
    # ... and "disturb" starting values a little.
    gf["A"] = gf["A"] + np.random.normal(0.0, 0.1)
    gf["mu"] = gf["mu"] + np.random.normal(0.0, 0.1)
    gf["sig"] = gf["sig"] + np.random.normal(0.0, 0.03)
    # Find the best fit solution
    gf.fit(x, y, yerr=np.ones(len(x))*0.1)
    
    # Step the amplitude (area of the Gaussian) through
    # the range 0.8 to 0.95 in 20 steps. Note that the
    # last part of `ranges` ('lin') is optional. You may
    # also use `log`; in this case, the stepping would be
    # equidistant in the logarithm.
    # In each step of `A`, "mu" and "sig" will be fitted,
    # because they had been thawed earlier.
    sp = gf.steppar("A", ranges={"A":[0.8, 0.95, 20, 'lin']})
    # Extract the values for the Gaussian normalization
    # (amplitude) ...
    As = list(map(lambda x:x[0], sp))
    # ... and chi square.
    chis = list(map(lambda x:x[1], sp))
    
    
    # Find minimum chi square
    cmin = min(chis)
    
    # Plot A vs. chi square
    plt.title('A vs. $\chi^2$ with 68% and 90% confidence levels')
    plt.xlabel("A")
    plt.ylabel("$\chi^2$")
    plt.plot(As, chis, 'bp-')
    plt.plot(As, [cmin+1.0]*len(As), 'k--')
    plt.plot(As, [cmin+2.706]*len(As), 'k:')
#    plt.show() 

  def sanity_steppar2(self):
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit as fuf
    
    # Set up a Gaussian model
    # and create some "data"
    x = np.linspace(0,2,100)
    gf = fuf.GaussFit1d()
    gf["A"] = 0.87
    gf["mu"] = 1.0
    gf["sig"] = 0.2
    y = gf.evaluate(x)
    y += np.random.normal(0.0, 0.1, len(x))
    
    # Thaw parameters, which are to be fitted ...
    gf.thaw(["A", "mu", "sig"])
    # ... and "disturb" starting values a little.
    gf["A"] = gf["A"] + np.random.normal(0.0, 0.1)
    gf["mu"] = gf["mu"] + np.random.normal(0.0, 0.1)
    gf["sig"] = gf["sig"] + np.random.normal(0.0, 0.03)
    # Find the best fit solution
    gf.fit(x, y, yerr=np.ones(len(x))*0.1)
    
    # Step the amplitude (area of the Gaussian) and the
    # center ("mu") of the Gaussian through the given
    # ranges.
    sp = gf.steppar(["A", "mu"], ranges={"A":[0.8, 0.95, 20], \
                    "mu":[0.96,1.05,15]})
    
    # Get the values for `A`, `mu`, and chi-square
    # from the output of steppar.
    As = list(map(lambda x:x[0], sp))
    mus = list(map(lambda x:x[1], sp))
    chis = list(map(lambda x:x[2], sp))
    
    # Create a chi-square array using the
    # indices contained in the output.
    z = np.zeros((20, 15))
    for s in sp:
      z[s[3]] = s[2]
    
    # Find minimum chi-square and define levels
    # for 68%, 90%, and 99% confidence intervals.
    cm = min(chis)
    levels = [cm+2.3, cm+4.61, cm+9.21]
    
    # Plot the contours to explore the confidence
    # interval and correlation.
    plt.xlabel("mu")
    plt.ylabel("A")
    plt.contour(np.sort(np.unique(mus)), np.sort(np.unique(As)), z, \
                levels=levels)
    # Plot the input value
    plt.plot([1.0], [0.87], 'k+', markersize=20)
#    plt.show()

  def sanity_errorConfInterval(self):
    """
      Checking example of errorConfInterval
    """
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit as fuf
    
    # Set up a Gaussian model
    # and create some "data"
    x = np.linspace(0,2,100)
    gf = fuf.GaussFit1d()
    gf["A"] = 0.87
    gf["mu"] = 1.0
    gf["sig"] = 0.2
    y = gf.evaluate(x)
    y += np.random.normal(0.0, 0.1, len(x))
    
    # Thaw parameters, which are to be fitted. Note
    # that those parameters will also be fitted during
    # the stepping; no further parameters will be thawed.
    gf.thaw(["A", "mu", "sig"])
    # ... and "disturb" starting values a little.
    gf["A"] = gf["A"] + np.random.normal(0.0, 0.1)
    gf["mu"] = gf["mu"] + np.random.normal(0.0, 0.1)
    gf["sig"] = gf["sig"] + np.random.normal(0.0, 0.03)
    # Find the best fit solution
    gf.fit(x, y, yerr=np.ones(len(x))*0.1)
    
    # Step the amplitude (area of the Gaussian) through
    # the range 0.8 to 0.95 in 20 steps. Note that the
    # last part of `ranges` ('lin') is optional. You may
    # also use `log`; in this case, the stepping would be
    # equidistant in the logarithm.
    # In each step of `A`, "mu" and "sig" will be fitted,
    # because they had been thawed earlier.
    sp = gf.steppar("A", ranges={"A":[0.8, 0.95, 20, 'lin']})
    # Extract the values for the Gaussian normalization
    # (amplitude) ...
    As = [x[0] for x in sp]
    # ... and chi square.
    chis = [x[1] for x in sp]
    
    # Calculate the confidence interval automatically
    cfi90 = gf.errorConfInterval("A", dstat=2.706)
    print("90% Confidence interval: ", cfi90["limits"])
    print("  corresponding objective function values: ", cfi90["OFVals"])
    print("  number of iterations needed: ", cfi90["iters"])
    
    cfi68 = gf.errorConfInterval("A", dstat=1.0)
    print("68% Confidence interval: ", cfi68["limits"])
    print("  corresponding objective function values: ", cfi68["OFVals"])
    print("  number of iterations needed: ", cfi68["iters"])
    
    # Plot A vs. chi square
    plt.title('A vs. $\chi^2$ 90% (black) and 68% (blue) confidence intervals')
    plt.xlabel("A")
    plt.ylabel("$\chi^2$")
    plt.plot(As, chis, 'bp-')
    # Indicate confidence levels by vertical lines
    plt.plot(As, [cfi90["OFMin"] +1.0]*len(As), 'g:')
    plt.plot(As, [cfi90["OFMin"]+2.706]*len(As), 'g:')
    # PLot lines to indicate confidence intervals
    plt.plot([cfi90["limits"][0]]*2, [min(chis), max(chis)], 'k--')
    plt.plot([cfi90["limits"][1]]*2, [min(chis), max(chis)], 'k--')
    plt.plot([cfi68["limits"][0]]*2, [min(chis), max(chis)], 'b--')
    plt.plot([cfi68["limits"][1]]*2, [min(chis), max(chis)], 'b--')
    
#     plt.show()


  def sanity_conditionalRestrictions(self):
    """
      Check the conditional restriction example.
    """
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit as fuf
    
    # Get fitting object for a Gaussian ...
    g = fuf.GaussFit1d()
    # .. and define the parameters
    g["A"] = 0.97
    g["mu"] = 0.1
    g["sig"] = 0.06
    
    # Generate some "data" with noise included
    x = np.linspace(-1.0,1.0,200)
    y = g.evaluate(x) + np.random.normal(0.0, 0.1, len(x))
    yerr = np.ones(len(x)) * 0.1
    
    
    def myRestriction(A, sig):
      """
        A conditional restriction.
        
        Returns
        -------
        Penalty : float
            A large value if condition is violated
            and zero otherwise.
      """
      if A > 10.0*sig:
        return np.abs(A-10.0*sig + 1.0)*1e20
      return 0.0
    
      
    # Add the conditional restriction to the model and save
    # the unique ID, which can be used to refer to that
    # restriction.
    uid = g.addConditionalRestriction(["A", "sig"], myRestriction)
    print("Conditional restriction has been assigned the ID: ", uid)
    print()
    
    # Now see whether the restriction is really in place
    g.showConditionalRestrictions()
    
    # Define free parameters ...
    g.thaw(["A", "mu", "sig"])
    # ... and fit the model (restriction included)
    g.fit(x, y, yerr=yerr)
    
    # Save the resulting best-fit model
    restrictedModel = g.model.copy()
    
    # Remove the conditional restriction and re-fit
    g.removeConditionalRestriction(uid)
    g.fit(x, y, yerr=yerr)
    
    # Save new model
    unrestrictedModel = g.model.copy()
    
    # Plot the result
#     plt.errorbar(x, y, yerr=yerr, fmt='b.')
#     plt.plot(x, restrictedModel, 'r--', label="Restricted")
#     plt.plot(x, unrestrictedModel, 'g--', label="Unrestricted")
#     plt.legend()
#     plt.show()



class MCMCExampleSanity(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    try:
      os.remove("mcmcExample.tmp")
    except:
      print("Could not remove file: mcmcExample.tmp")
    try:
      os.remove("mcmcTA.tmp")
    except:
      print("Could not remove file: mcmcTA.tmp")
    try:
      os.remove("mcmcSample.tmp")
    except:
      print("Could not remove file: mcmcSample.tmp")
    try:
      os.remove("chain.emcee")
    except:
      pass

  def sanity_MCMCSampler(self):    
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
#     plt.show()



  def sanity_MCMCPriorExample(self):
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
    
    #plt.show()
    
  def sanity_autoMCMCExample1(self):
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
#     plt.show()


  def sanity_autoMCMCExample2(self):
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
#    plt.show()


  def sanity_TAtut_createTrace(self):
    """
      TA tutorial, all examples
    """
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
#    plt.show()
    
    # Thaw the parameters, which we wish to vary
    # during the sampling
    vp.thaw(["A", "al", "mu", "ad"])
    
    # Use current parameters as starting point for the sampling
    X0 = vp.freeParameters()
    print("Starting point for sampling: ", X0)
    
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

    
    ######## Second example
    
    from PyAstronomy import funcFit as fuf
    
    # Create an instance of TraceAnalysis
    # telling it which file to use
    ta = fuf.TraceAnalysis("mcmcTA.tmp")
    
    # Have a look at the deviance to check if and when
    # the chains reached equilibrium.
    ta.plotTrace("deviance")
#    ta.show()
    
    # Say, we are sure that after 500 iterations, the chain
    # reached equilibrium. We use this as the burn-in phase
    ta.setBurn(500)
    
    # Have a second look at the deviance, this time considering
    # the burn-in. Note that the first 500 iterations are not
    # removed from the chain. They are just not considered any
    # more.
    ta.plotTrace("deviance")
#    ta.show()

    
    ######## Third example
    
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

    ######## Fourth example

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
#     ta.show()
    
    ######## Fifth example

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
#     plt.show()
    
    ######## Sixth example
    
    from PyAstronomy import funcFit as fuf
    
    # Create an instance of TraceAnalysis
    # telling it which file to use
    ta = fuf.TraceAnalysis("mcmcTA.tmp")
    
    # Use the burn-in from the previous example
    ta.setBurn(500)
    
    # Investigate a trace
    ta.plotTrace("mu")
#    ta.show()
    # and its distribution.
    ta.plotHist("mu")
#    ta.show()
    # Combine trace and distribution
    ta.plotTraceHist("mu")
#    ta.show()
    # Plot correlations
    ta.plotCorr(parsList=["mu", "ad", "al"])
#    ta.show()

  def sanity_MCMCautoParameters(self):
    """
      Checking sanity of MCMCautoParameters
    """
    from PyAstronomy import funcFit as fuf
    import numpy as np
    import matplotlib.pylab as plt
    
    x = np.linspace(0,30,1000)
    gauss = fuf.GaussFit1d()
    gauss["A"] = 1
    gauss["mu"] = 23.
    gauss["sig"] = 0.5
    yerr = np.random.normal(0., 0.05, len(x))
    y = gauss.evaluate(x) + yerr
    # This step is not necessary if <picky>=False in MCMCautoParameters.
    gauss.thaw(["A","mu","sig"])
    X0, lims, steps = gauss.MCMCautoParameters({"A":[0,10],"mu":3, "sig":[0.1,1.0]})
    gauss.fitMCMC(x, y, X0, lims, steps, yerr=yerr, iter=1000)
    
#     plt.plot(x, y, 'k+')
#     plt.plot(x, gauss.evaluate(x), 'r--')
#     plt.show()

  def sanity_EMCEEfirstexample(self):

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
    
    # Start the sampling (ps could be used to continueb the sampling)
    ps = gf.fitEMCEE(x, y, yerr=ones(x.size)/snr, sampleArgs=sampleArgs)
    
    # Plot the distributions of the chains
    # NOTE: the order of the parameters in the chain object is the same
    #       as the order of the parameters returned by freeParamNames()
    for i, p in enumerate(gf.freeParamNames()):
      plt.subplot(len(gf.freeParamNames()), 1, i+1)
      plt.hist(gf.emceeSampler.flatchain[::,i], label=p)
      plt.legend()
#     plt.show()


  def sanity_EMCEEpriorexample(self):
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
#     plt.show()
    
    
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
#     plt.show()


  def sanity_InstatiatePrior(self):
    from PyAstronomy import funcFit as fuf
    
    # Instantiate prior
    gp = fuf.FuFPrior("gaussian", sig=0.1, mu=1.0)
    
    # Current values (arbitrary)
    cvals = {"a":1.4, "b":0.86, "c":1.1}
    
    # Get log(prior) for parameter "b"
    print(gp(cvals, "b"))