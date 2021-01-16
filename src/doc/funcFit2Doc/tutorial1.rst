The funcFit2 tutorial
=====================

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

.. note:: Please mind the citation request for use of the scipy algorithms explained
          `here <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html>`_.
  

Prerequisites
-------------
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


Example: Fitting a Gaussian
---------------------------

* :doc:`ex_gaussfit1` :download:`(Download notebook) <ex_gaussfit1.ipynb>`

On the objective function
-------------------------

The objective function is the function, the value of which is to be minimized.
Models can but do not have to have a predefined objective function. For example, the
GaussFit model uses chi square as a default objective functions. 

::

    from __future__ import print_function, division
    from PyAstronomy import funcFit2 as fuf2
    
    # Create a model object
    gf = fuf2.GaussFit()
    
    print("Information on the objective function:")
    print("    ", gf.objfInfo())

Whose answer reads:

::

    Information on the objective function:
        chi-square 
          


Restrictions of parameter ranges
--------------------------------

A *restriction* in funcFit2 limits the valid range of values of a parameter. Restrictions are common as boundary
conditions in modeling and problems of optimization. 
For example, the width (standard deviation) of a Gaussian must be positive or certain spectral
lines must only occur in absorption or emission for physical reasons.

Restrictions can be handled in many ways. The possibilities to implement restrictions include:

- **Parameter transformations**
  The restriction can be absorbed in the definition of the model, e.g., by using the absolute value of the
  standard deviation in calculating a Gaussian curve. Also mapping of the real number onto a finite range
  such as the `sigmoid function <https://en.wikipedia.org/wiki/Sigmoid_function>`_ can be used.
- **Optimization algorithms** Some optimization algorithm allow to specify boundaries (or more general constraints)
  for the parameter values. One example of such an algorithm is scipy's
  `fmin_l_bfgs_b <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html>`_.
- **Penalty functions** Restrictions can be implemented by penalizing the objective function when the
  boundaries are violated.
  If combined with an optimization algorithm based on gradient descent, it is often helpful to implement
  "soft edges" for penalty, i.e., a strong but finite gradient in the objective, which allows the algorithm
  to "find its way back".


Use penalties to enforce restrictions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Penalties are funcFit2's default mechanism to account for restrictions. A restriction is specified
by a two-valued list with a lower and an upper limit for the parameter range. None can be used to
indicate that no boundary applies on one (or both) sides. Restrictions can be added
to parameters via the `setRestriction` method.

If a parameter value violates the specified boundary restrictions by some margin x, a value of
abs(x)*penaltyFactor is added to the value of the objective function. The default value of
penaltyFactor is 1e20.

* :doc:`ex_restrictviapenalties` :download:`(Download notebook) <ex_restrictviapenalties.ipynb>`


Algorithm allowing boundary conditions (fmin_l_bfgs_b)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here directly invoke the
`fmin_l_bfgs_b <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html>`_
as implemented in scipy to carry out an optimization with boundary conditions

* :doc:`ex_algorestrict` :download:`(Download notebook) <ex_algorestrict.ipynb>`

Use a convenience function to automatically channel the restrictions from the model to the
algorithm 

* :doc:`ex_algorestrict_con` :download:`(Download notebook) <ex_algorestrict_con.ipynb>`



Use a custom objective function
-------------------------------

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
------------------

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
-------------

Using custom models is easy.

* :doc:`ex_linmod1` :download:`(Download notebook) <ex_linmod1.ipynb>`
* :doc:`ex_linmod_jit` :download:`(Download notebook) <ex_linmod_jit.ipynb>`


Combine models
--------------

::
        
    from __future__ import print_function, division
    import numpy as np
    import matplotlib.pylab as plt
    # ... and now the funcFit package
    from PyAstronomy import funcFit2 as fuf2
    import scipy.optimize as sco
    
    # Instantiate Gaussian model objects
    gf1 = fuf2.GaussFit1d()
    gf2 = fuf2.GaussFit1d()
    # Sum the models (refers to their 'evaluate' methods)
    # Any of +-*/ and ** can be used
    gf = gf1 + gf2
    # Use the Gaussian likelihood
    gf.setlogL("1dgauss")
    
    gf.parameterSummary()
    
    gf["A_GF1d(1)"] = 1
    gf["A_GF1d(2)"] = 2
    gf["mu_GF1d(1)"] = 0.0
    gf["mu_GF1d(2)"] = 3.0
    gf["sig_GF1d(1)"] = 1.0
    gf["sig_GF1d(2)"] = 1.0
    
    # Evaluate model and add noise
    x = np.linspace(-4., 6., 200)
    y = gf.evaluate(x) + np.random.normal(0, 0.02, len(x))
    
    # Re-fit model
    # Use filename-like pattern matching to thaw parameters
    gf.thaw(["A_*", "sig_*", "mu_*"])
    fuf2.fitfmin_l_bfgs_b1d(gf, x, y, yerr=0.02)
    
    gf.parameterSummary()
    
    plt.errorbar(x, y, yerr=0.02, fmt="b+")
    plt.plot(x, gf.evaluate(x), 'r--')
    plt.show()
    

