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
    

Restrictions
----------------------------

A *restriction* limits the valid range of values of a parameter. Restrictions are common in modeling
and problems of optimization. 
For example, the width (standard deviation) of a Gaussian should be positive or certain spectral
lines must only occur in absorption or emission for physical reasons.

Restrictions can be handled in many ways:

- Some optimization algorithm allow to specify boundaries (or more general constraints) for the parameter
  values. One example of such an algorithm is scipy's "fmin_l_bfgs_b".
- Restrictions can be implemented by penalizing the objective function when the boundaries are violated.
  If combined with an optimization algorithm based on gradient descent, it is often helpful to implement
  "soft edges" for penalty, i.e., a strong but finite gradient in the objective, which allows the algorithm
  to "find its way back".
- In a Bayesian analysis, restrictions are naturally absorbed in the priors, when the posterior is the objective.


Algorithm allowing boundary conditions (fmin_l_bfgs_b)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here directly invoke scipy's fmin_l_bfgs_b to carry out an optimization with
boundary conditions.
  
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






