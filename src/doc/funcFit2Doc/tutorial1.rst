The funcFit2 tutorial
=======================

This tutorial is supposed to enable you to exploit the funcFit2 functionality
by presenting an anthology of examples.

Some knowledge of numpy, SciPy, and matplotlib is helpful to dive in.

.. contents:: Sections of this tutorial

.. _matplotlib: http://matplotlib.sourceforge.net/
.. _pymc: https://github.com/pymc-devs/pymc
.. _SciPy: www.scipy.org/
.. _numpy: numpy.scipy.org/
.. _XSPEC: http://heasarc.nasa.gov/xanadu/xspec/
.. _emcee: http://dan.iel.fm/emcee/current/

Prerequisites
-------------------
To run the example in this tutorial you need to have installed the following packages:
 * numpy_
 * SciPy_
 * matplotlib_

After installing PyAstronomy (PyA), funcFit2
is available for being used. 
As a first step, let us import the
package and check the status of scipy and emcee in out installation.
To that end,

::

  from PyAstronomy import funcFit2 as fuf2
  
  fuf2.status()

Depending on your installation the output should look like:

::

    Status of funcFit2:
    --------------------------
    Is scipy.optimize available?  yes
      Version:  1.1.0
    Is emcee available?  yes
      Version:  2.2.1


Fitting a Gaussian
-------------------------------

Now let us dive into business. The first example below demonstrates how to \
exploit the functionality of funcFit to fit a Gaussian to one-dimensional mock data. \
It shows how free parameters can be specified and restrictions can be applied.

::

    from __future__ import print_function, division
    # Import numpy and matplotlib
    from numpy import arange, sqrt, exp, pi, random, ones_like
    import matplotlib.pylab as plt
    # ... and now the funcFit package
    from PyAstronomy import funcFit2 as fuf2
    
    # Before we can start fitting, we need something to fit.
    # So let us create some data...
    
    # Creating a Gaussian with some noise
    # Choose some parameters...
    gPar = {"A":-5.0, "sig":10.0, "mu":10.0, "off":1.0, "lin":0.0}
    # Calculate profile
    x = arange(100) - 50.0
    y = gPar["off"] + gPar["A"] / sqrt(2*pi*gPar["sig"]**2) \
        * exp(-(x-gPar["mu"])**2/(2*gPar["sig"]**2))
    # Add some noise...
    y += random.normal(0.0, 0.01, x.size)
    # ...and save the error bars
    yerr = ones_like(x)*0.01
    # Let us see what we have done...
    plt.plot(x, y, 'bp')
    
    # Now we can start exploiting the funcFit2 functionality to
    # fit a Gaussian to our data. In the following lines, we
    # create a fitting object representing a Gaussian and set guess parameters.
    
    # Create a model object
    gf = fuf2.GaussFit1d()
    # See what parameters are available
    print("List of available parameters: ", gf.availableParameters())
    # Set guess values for the parameters
    gf["A"] = -10.0
    gf["sig"] = 15.77
    gf["off"] = 0.96
    gf["mu"] = 7.5
    # Let us see whether the assignment worked
    print("Parameters and guess values: ")
    print("  A   : ", gf["A"])
    print("  sig : ", gf["sig"])
    print("  off : ", gf["off"])
    print("  mu  : ", gf["mu"])
    print("")
    
    # Evaluate the model for the staring values
    startmodel = gf.evaluate(x)
    
    # Which parameters shall be varied during the fit?
    # 'Thaw' those (the order is irrelevant)
    gf.thaw(["A", "sig", "off", "mu"])
    
    # Show values and names of frozen parameters (i.e., parameters not varied in a fit)
    print("Names and values of FROZEN parameters: ", gf.frozenParameters())
    
    # Use a convenience function to carry out optimization using scipy's fmin
    fr = fuf2.fitfmin1d(gf, x, y, yerr=yerr)
    
    print()
    print("Fit result from scipy's fmin: ", fr)
    print()
    
    # Get a summary of current parameters from the model
    gf.parameterSummary()
    
    plt.plot(x, startmodel, 'g:', label="Model for starting values")
    plt.plot(x, gf.evaluate(x), 'r--', label="Best-fit model")
    plt.legend()
    plt.show()
    

On the objective function
----------------------------------

So what do we minimize in the previous example or elsewhere?

::

    from __future__ import print_function, division
    from PyAstronomy import funcFit2 as fuf2
    
    # Create a model object
    gf = fuf2.GaussFit1d()
    
    print("Information on the objective function:")
    print("    ", gf.objfInfo())

Whose answer reads:

::

    Information on the objective function:
          Objective function: -ln(Likelihood) 

The default objective function is the negative (natural) logarithm of
the likelihood.





Restrictions
----------------------------

A *restriction* limits the valid range of values of a parameter. Restrictions are common in modeling
and problems of optimization. 
For example, the width (standard deviation) of a Gaussian should be positive or certain spectral
lines must only occur in absorption or emission for physical reasons.

Restrictions can be handled in many ways:

- The restriction can be absorbed in the definition of the model, e.g., by using the absolute value of the
  standard deviation in calculating a Gaussian curve.
- Some optimization algorithm allow to specify boundaries (or more general constraints) for the parameter
  values. One example of such an algorithm is scipy's "fmin_l_bfgs_b".
- Restrictions can be implemented by penalizing the objective function when the boundaries are violated.
  If combined with an optimization algorithm based on gradient descent, it is often helpful to implement
  "soft edges" for penalty, i.e., a strong but finite gradient in the objective, which allows the algorithm
  to "find its way back".
- In a Bayesian analysis, restrictions are naturally absorbed in the priors, when the posterior is the objective.


Algorithm allowing boundary conditions (fmin_l_bfgs_b)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here directly invoke the
`fmin_l_bfgs_b <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html>`_
as implemented in scipy to carry out an optimization with boundary conditions.

.. note:: Please mind the citation request for use of the algorithm explained
          `here <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html>`_.
  
::
    
    from __future__ import print_function, division
    from numpy import arange, sqrt, exp, pi, random, ones_like
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit2 as fuf2
    import scipy.optimize as sco
    
    
    # Creating a Gaussian with some noise
    # Choose some parameters...
    gPar = {"A":-5.0, "sig":10.0, "mu":10.0, "off":1.0, "lin":0.0}
    # Calculate profile
    x = arange(100) - 50.0
    y = gPar["off"] + gPar["A"] / sqrt(2*pi*gPar["sig"]**2) \
        * exp(-(x-gPar["mu"])**2/(2*gPar["sig"]**2))
    # Add some noise...
    y += random.normal(0.0, 0.01, x.size)
    # ...and save the error bars
    yerr = ones_like(x)*0.01
    # Let us see what we have done...
    plt.plot(x, y, 'bp')
    
    # Create a model object
    gf = fuf2.GaussFit1d()
    
    # Set guess values for the parameters
    gf.assignValues({"A":-3, "sig":10.77, "off":0.96, "mu":10.5})
    
    # 'Thaw' those (the order is irrelevant)
    gf.thaw(["mu", "sig", "off", "A"])
    
    # We need the order to get the order of bounds right
    # This is not necessarily the order in which they are thawed!
    print("Free parameter names and their order: ", gf.freeParamNames())
    
    # Use fmin_l_bfgs_b with amplitude restricted to the (-2,0) interval
    fr = sco.fmin_l_bfgs_b(gf.objf, gf.freeParamVals(), args=(x,y,yerr), bounds=((-2.,0), (None,None), (None,None), (None,None)), approx_grad=True)
    # Set the parameter values to best-fit
    gf.setFreeParamVals(fr[0])
    
    gf.parameterSummary()
    plt.plot(x, gf.evaluate(x), 'r--')
    plt.show()


Built-in restrictions with convenience function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Convenience functions can automatically grab the restrictions from the model
and hand them to the minimization algorithm.

::
    
    from __future__ import print_function, division
    from numpy import arange, sqrt, exp, pi, random, ones_like
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit2 as fuf2
    import scipy.optimize as sco
    
    
    # Creating a Gaussian with some noise
    # Choose some parameters...
    gPar = {"A":1.0, "sig":10.0, "mu":10.0, "off":1.0, "lin":0.0}
    # Calculate profile
    x = arange(100) - 50.0
    y = gPar["off"] + gPar["A"] / sqrt(2*pi*gPar["sig"]**2) \
        * exp(-(x-gPar["mu"])**2/(2*gPar["sig"]**2))
    # Add some noise...
    y += random.normal(0.0, 0.002, x.size)
    # ...and save the error bars
    yerr = ones_like(x)*0.002
    # Let us see what we have done...
    plt.plot(x, y, 'bp')
    
    # Create a model object
    gf = fuf2.GaussFit1d()
    
    # Set guess values for the parameters
    gf.assignValues({"A":3, "sig":3.77, "off":0.96, "mu":10.5})
    
    # 'Thaw' those (the order is irrelevant)
    gf.thaw(["mu", "sig", "off", "A"])
    
    # Restrict the valid range for sigma
    gf.setRestriction({"sig":[0,5]})
    
    # The convenience function 'fitfmin_l_bfgs_b1d' automatically channels
    # the restrictions from the model to the algorithm.
    fuf2.fitfmin_l_bfgs_b1d(gf, x, y, yerr=yerr)
    
    gf.parameterSummary()
    plt.plot(x, gf.evaluate(x), 'r--')
    plt.show()


Use penalties for boundary violations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Penalties can be used to implement restrictions. Here we use `objfPenalize`
to add penalties to the objective function.

::
    
    from __future__ import print_function, division
    from numpy import arange, sqrt, exp, pi, random, ones_like
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit2 as fuf2
    import scipy.optimize as sco
    
    
    # Creating a Gaussian with some noise
    # Choose some parameters...
    gPar = {"A":1.0, "sig":10.0, "mu":10.0, "off":1.0, "lin":0.0}
    # Calculate profile
    x = arange(100) - 50.0
    y = gPar["off"] + gPar["A"] / sqrt(2*pi*gPar["sig"]**2) \
        * exp(-(x-gPar["mu"])**2/(2*gPar["sig"]**2))
    # Add some noise...
    y += random.normal(0.0, 0.002, x.size)
    # ...and save the error bars
    yerr = ones_like(x)*0.002
    # Let us see what we have done...
    plt.plot(x, y, 'bp')
    
    # Create a model object
    gf = fuf2.GaussFit1d()
    
    # Set guess values for the parameters
    gf.assignValues({"A":3, "sig":3.77, "off":0.96, "mu":9.5})
    
    # 'Thaw' those (the order is irrelevant)
    gf.thaw(["mu", "sig", "off", "A"])
    
    # Restrict parameter ranges
    gf.setRestriction({"sig":[0,7]})
    
    # Use chi-square is objective
    gf.objfnChiSquare()
    # Apply penalties for violating boundaries
    gf.objfPenalize()
    
    # Use a minimization algorithm not accounting for boundaries
    # with penalties
    fr = fuf2.fitfmin1d(gf, x, y, yerr=yerr)
    print("Fit result: ", fr)
    
    gf.parameterSummary()
    plt.plot(x, gf.evaluate(x), 'r--')
    plt.show()


Use a custom objective function
-------------------------------------

::
    
    from __future__ import print_function, division
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit2 as fuf2
    import scipy.optimize as sco
    
    
    # Creating a Gaussian with some noise
    # Choose some parameters...
    gPar = {"A":1.0, "sig":10.0, "mu":10.0, "off":1.0, "lin":0.0}
    # Calculate profile
    x = np.arange(100) - 50.0
    y = gPar["off"] + gPar["A"] / np.sqrt(2*np.pi*gPar["sig"]**2) \
        * np.exp(-(x-gPar["mu"])**2/(2*gPar["sig"]**2))
    # Add some noise...
    y += np.random.normal(0.0, 0.002, x.size)
    # ...and save the error bars
    yerr = np.ones_like(x)*0.002
    # Let us see what we have done...
    plt.plot(x, y, 'bp')
    
    # Create a model object
    gf = fuf2.GaussFit1d()
    
    # Set guess values for the parameters
    gf.assignValues({"A":3, "sig":3.77, "off":0.96, "mu":9.5})
    
    # 'Thaw' those (the order is irrelevant)
    gf.thaw(["mu", "sig", "off", "A"])
    
    
    def myobjf(m, pars, x, y, **kwargs):
        """
        Calculate the absolute distance between model and data points and area of the Gaussian
        """
        model = m.evaluate(x)
        r = np.sum(np.abs(model - y))
        # Weird dependence on A
        r += 0.1*np.abs(m["A"])
        # Some informative output
        fpv = tuple(zip((m.freeParamNames(), m.freeParamVals())))
        print("Parameter values: ", fpv)
        print("    identical to those received: ", pars)
        print("Objective: ", r)
        return r
    
    # Assign objective function. As a result of this assignment,
    # the parameter values of the model will be set to those handed
    # to the objective function prior to call
    gf.objf = myobjf
    
    fr = sco.fmin(gf.objf, gf.freeParamVals(), args=(x,y))
    print("Fit result: ", fr)
    # Set the parameter values to best-fit
    gf.setFreeParamVals(fr)
    
    gf.parameterSummary()
    plt.plot(x, gf.evaluate(x), 'r--')
    plt.show()


Applying relations
--------------------------

Relations define functional dependences between different parameter values.

::
    
    from __future__ import print_function, division
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit2 as fuf2
    import scipy.optimize as sco
    
    
    # Creating a Gaussian with some noise
    # Choose some parameters...
    gPar = {"A":1.0, "sig":10.0, "mu":10.0, "off":1.0, "lin":0.0}
    # Calculate profile
    x = np.arange(100) - 50.0
    y = gPar["off"] + gPar["A"] / np.sqrt(2*np.pi*gPar["sig"]**2) \
        * np.exp(-(x-gPar["mu"])**2/(2*gPar["sig"]**2))
    # Add some noise...
    y += np.random.normal(0.0, 0.002, x.size)
    # ...and save the error bars
    yerr = np.ones_like(x)*0.002
    # Let us see what we have done...
    plt.plot(x, y, 'bp')
    
    # Create a model object
    gf = fuf2.GaussFit1d()
    
    # Set guess values for the parameters
    gf.assignValues({"A":3, "sig":3.77, "off":0.96, "mu":9.5})
    
    # 'Thaw' those (the order is irrelevant)
    gf.thaw(["mu", "sig", "off", "A"])
    
    def relat(sig, off):
        """ Combine values of sig and off """
        return 0.1*sig - off
    
    # 'A' is a function of 'sig' and 'off'.
    gf.relate("A", ["sig", "off"], func=relat)
    
    fr = sco.fmin(gf.objf, gf.freeParamVals(), args=(x,y))
    print("Fit result: ", fr)
    # Set the parameter values to best-fit
    gf.setFreeParamVals(fr)
    
    gf.parameterSummary()
    plt.plot(x, gf.evaluate(x), 'r--')
    plt.show()


Custom models
---------------------

Straight line
~~~~~~~~~~~~~~~~~

::

    from __future__ import print_function, division
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit2 as fuf2
    import scipy.optimize as sco
    
    
    class LinMod(fuf2.MBO2):
        """ Linear model with additional jitter """
        
        def __init__(self):
            # 'pars' specifies parameter names in the model
            fuf2.MBO2.__init__(self, pars=["const", "slope"], rootName="LinMod")
            # Use likelihood based on Gaussian errors with std yerr
            self.setlogL("1dgauss")
        
        def evaluate(self, x):
            """ Evaluate line """
            return self["const"] + x * self["slope"]
    
    
    # Instantiate model
    lm = LinMod()
    lm["slope"] = 1.1
    lm["const"] = -0.5
    
    # Get some 'data' and add Gaussian noise with STD 10
    x = np.arange(15.)
    y = lm.evaluate(x) + np.random.normal(0,1,len(x))
    yerr = np.ones_like(x)
    
    lm.thaw(["slope", "const"])
    
    fr = sco.fmin(lm.objf, x0=lm.freeParamVals(), args=(x,y,yerr))
    lm.setFreeParamVals(fr)
    
    lm.parameterSummary()
    
    plt.errorbar(x, y, yerr=yerr, fmt='b+')
    plt.plot(x, lm.evaluate(x), 'r--')
    plt.show()



Straight line with jitter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A linear model with custom specification of the likelihood, which allows to take into account
an additional contribution to the uncertainty of the data (jitter). The jitter can be treated
as a free parameter in the fit. 

::

    from __future__ import print_function, division
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit2 as fuf2
    import scipy.optimize as sco
    
    
    class LinMod(fuf2.MBO2):
        """ Linear model with additional jitter """
        
        def __init__(self):
            fuf2.MBO2.__init__(self, pars=["const", "slope", "jitter"], rootName="LinMod")
        
        def evaluate(self, x):
            """ Evaluate line """
            return self["const"] + x * self["slope"]
    
        def logL(self, x, y, yerr, **kwargs):
            """ ln(Likelihood) with jitter as additional term """
            yr = np.sqrt(yerr**2 + self["jitter"]**2)
            m = self.evaluate(x)
            lnl = -len(x)/2.0*np.log(2.*np.pi) - np.sum(np.log(yr)) - 0.5 * np.sum((m-y)**2/(yr**2))
            return lnl
    
    # Instantiate model
    lm = LinMod()
    lm["slope"] = 1.1
    lm["const"] = -0.5
    
    # Get some 'data' and add Gaussian noise with STD 10
    x = np.arange(150.)
    y = lm.evaluate(x) + np.random.normal(0,10,len(x))
    # Nominal error has STD 1
    yerr = np.ones_like(x)
    
    lm.thaw(["slope", "const", "jitter"])
    
    fr = sco.fmin(lm.objf, x0=lm.freeParamVals(), args=(x,y,yerr))
    lm.setFreeParamVals(fr)
    
    lm.parameterSummary()
    
    plt.errorbar(x, y, yerr=yerr, fmt='b+')
    plt.plot(x, lm.evaluate(x), 'r--')
    plt.show()

