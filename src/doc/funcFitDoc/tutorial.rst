The funcFit tutorial
=======================

This tutorial is supposed to enable you to exploit the funcFit functionality \
by presenting examples and delineate the basic principles.

Please note that this tutorial is not intended to be an introduction to Python or \
the numpy, SciPy, or matplotlib package. Basic knowledge of these components \
is desirable.

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
 
For MCMC sampling, at least one of the following packages should be installed
 
 * pymc_ (in version *2.x*, unfortunately, version 3.x is not supported)
 * emcee_

We also assume that you have a basic knowledge of Python and these packages (except for
pymc).

.. warning:: funcFit relies on **Python 2.7.x**. The 2.6.x series (and prior) has a bug affecting the copying of
             dynamically created class methods, which has not been (and will not be) corrected. This interferes with many of
             funcFit's algorithms.

What is it good for?
------------------------
Popular Python packages such as SciPy already offer implementations of the most common \
fitting algorithms. Although a great deal can be achieved with them alone, some \
specific tasks need considerable effort to be achieved. What, for example, if you wish to \
try different sets of parameters, maybe freeze some of them at a specific value for now and let them \
vary next? \
This is possible, yet probably not quite convenient. What if the parameters are \
restricted or even related by some function? 

These are the problems to be approached by the funcFit package. Basically, it is a wrapper \
around minimization routines defined in other packages (e.g., SciPy).
The whole *magic* \
is to provide a framework to handle parameters and provide a convenient interface; \
those of you familiar with XSPEC_ will probably experience a déjà vu. \
Read on to find out how it \
works and how it can make your work easier.


.. note:: In this tutorial, it is assumed that you also have matplotlib installed,
          which provides a neat interface for plotting under Python.
          If you did not already install it, have a look at their we page (matplotlib_).

Diving into the first example
-------------------------------

After you have installed PyAstronomy (PyA), the funcFit package
is ready for being used. As a very first step, let us import the
package and see whether we will be able to actually fit something.
Therefore, execute the following lines as a script or from the
Python interactive command shell: 

::

  from PyAstronomy import funcFit as fuf
  
  fuf.status()

Depending on your installation the output should look like:

::

    Status of funcFit:
    --------------------------
    Is scipy.optimize available?  yes
    Is pymc available?  yes
      pymc is available in version:  2.3.4
    Is emcee available?  yes
    

Now let us dive into the business of fitting. The first example shown below demonstrates how to \
exploit the functionality of funcFit to fit a Gaussian to artificially created data. \
It shows how free parameters can be specified and restrictions can be applied.

::

    from __future__ import print_function, division
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
    plt.show()



Running the above script yields the following output (numbers may slightly differ):

::

    List of available parameters:  ['A', 'mu', 'lin', 'sig', 'off']
    Parameters and guess values: 
      A   :  -10.0
      sig :  15.77
      off :  0.87
      mu  :  7.5
    
    Names and values of FROZEN parameters:  {'A': -10.0, 'mu': 7.5, 'lin': 0.0, 'sig': 15.77, 'off': 0.87}
    Optimization terminated successfully.
             Current function value: 111.455503
             Iterations: 176
             Function evaluations: 310
    ----------------------------------
    Parameters for Component: Gaussian
    ----------------------------------
    Parameter:   A  Gaussian, [  A], value:     -4.92037, free:  True, restricted:  True, related: False
        Restriction: [None,  0]
    Parameter:  mu  Gaussian, [ mu], value:      9.83938, free:  True, restricted: False, related: False
    Parameter: lin  Gaussian, [lin], value:            0, free: False, restricted: False, related: False
    Parameter: sig  Gaussian, [sig], value:      9.97104, free:  True, restricted: False, related: False
    Parameter: off  Gaussian, [off], value:     0.999786, free:  True, restricted: False, related: False


Some points in the example shall be emphasized:
  * The names of the parameters are defined by the fitting object (in this case GaussFit1d),
  * Parameter values can be set and obtained using brackets,
  * Individual parameters can be thawed or frozen depending on the needs of the user,
  * Restrictions on the parameter ranges can be applied either on both or just on side of the range,
  * After the fit, the best-fit values become the current parameters, i.e., they can be obtained using the bracket operator,
  * After the fit, the best-fit model can be accessed through the `model` property. 

The central step of the script is the call to *fit*. The method takes at least two arguments: the
x-axis and corresponding y-axis values; errors on the y-axis values can be given optionally
via the *yerr* keyword as shown in the example. In default configuration, the *fit* method
uses the *fmin* routine provided by SciPy.optimize to minimize either the sum of quadratic residuals
of no error is provided, or :math:`\chi^2` if errors (yerr)
are given.

.. note:: Restrictions are implemented using a **penalty function**. The steepness of the penalty
          may be changed by the *setPenaltyFactor* method or by accessing the `penaltyFactor`
          property directly.


Introducing a custom model
-------------------------------

The funcFit package comes with some fitting models, but in many cases it will be necessary
to use custom models. Introducing a new model is easy in funcFit and will be demonstrated
in the next example. Here we implement a straight line and fit it to some artificial data.

::

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
    plt.show()

This example resembles the first one, but here we defined a custom fitting
model at the top instead of using the *GaussFit1d* class as in the first example.

A new fitting model is a class, which inherits from the *OneDFit* class. Additionally, two
methods (*__init__* and *evaluate*) must be implemented.
In the example, we
provide a minimal constructor (__init__ method), which only consists of a call to the
base class (OneDFit) constructor. The argument is a list of strings with the names of the
variables characterizing the model. The *evaluate* method takes a single argument, which is
an array of values at which to evaluate the model. It returns the function values at the
given position. Note how, e.g., *self["off"]*, is used to get the current value if the offset
variable in *evaluate*.

Applying *relations*
----------------------
In funcFit *relations* refer to a functional dependence between two or more model parameters.
To demonstrate the application of such a relation, we slightly extend the previous example.
In particular, we will assume that the gradient of our line is a multiple of the offset.

::
  
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
    plt.show()
  

The output of the script reads (numbers may differ):

::

    Optimization terminated successfully.
             Current function value: 251539.530679
             Iterations: 27
             Function evaluations: 54
    ---------------------------------
    Parameters for Component: unnamed
    ---------------------------------
    Parameter: lin  , [lin], value:       3.5004, free: False, restricted: False, related:  True
         Relation: lin = f(off)
    Parameter: off  , [off], value:     0.388933, free:  True, restricted: False, related: False


.. note:: The `lin` parameter is no longer free, as it depends on `off`.

The *relate* method takes three arguments. The first is the name of
the dependent variable (in this case "lin"). The second is a list
containing the names of the independent variables (in this case
only "off"). The third argument is a callable object, which provides
the numerical relation between the independent and the dependent
variables (there may be more than one independent variable).


Combining models  
------------------

The funcFit package allows to combine two models. That means that models (then becoming model components) can
be added, subtracted, divided, multiplied, and even used as exponents. This can be very useful in creating
more complex models and requires only little effort. The following example shows how two Gaussians models
can be summed.

::
       
    from __future__ import print_function, division
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
    plt.show()

    
.. note:: `twoG` contains copies (not references) two its "ancestors" (`gf1` and `gf2`). You can, thus, continue using those as usual.

When the models are combined (added in this case), funcFit adds "component identifiers" to the variable names to ensure that they
remain unique. A component identifier is simply an appendix to the variable name consisting of an underscore, the model name,
and a number. The combined
model behaves exactly like the individual ones. It should also be noted that model characteristics such as relations, restrictions, etc.,
are preserved in the combined model.


Applying conditional restrictions
--------------------------------------

Via `conditional restrictions` complex penalty (or reward) functions can be
defined, which keep the fit out or force it into a specific subspace of the
parameter space. Conditional restrictions are self-defined callables such
as function, which take a number of parameters and return a float, which
specifies the penalty. The latter is added to the objective function.

Conditional restrictions are referred to by a unique ID, which is generated
as soon as it is added to the model. Note that this ID does not change, when
models are combined.    

::
    
    from __future__ import print_function, division
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
    plt.errorbar(x, y, yerr=yerr, fmt='b.')
    plt.plot(x, restrictedModel, 'r--', label="Restricted")
    plt.plot(x, unrestrictedModel, 'g--', label="Unrestricted")
    plt.legend()
    plt.show()



Minimize the Cash statistic
----------------------------

In many cases, the use of the :math:`\chi^2` statistic is inappropriate.
If, for instance, the data consist of only a few counts per bin, using
the Cash statistic (Cash 1979, ApJ 228, 939) can be more appropriate.
Built-in statistics can be used by specifying the `miniFunc` parameter
on call to fit, as is demonstrated in the following example.

::

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
  plt.show()
 

Using "steppar" to determine confidence intervals
---------------------------------------------------

The "steppar" command can be used to analyze the
behavior of the objective function (e.g., :math:`\chi^2`)
as the parameter values are varied. In particular, the specified
parameter(s) are set to a number of values and the remaining
free parameters are fitted.

The example below shows how to determine a confidence
interval for the normalization of a Gaussian.

::

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
  plt.show()



The next example demonstrates how to step two parameters
through given ranges and plot the resulting confidence
contours.

::

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
  plt.show()


Use `errorConfInterval` to determine confidence intervals
-----------------------------------------------------------

The `steppar` example shows how confidence intervals may be
estimated by exploring the behavior of the objective function
manually. The `errorConfInterval` strives to find the confidence
interval automatically.

::
     
    from __future__ import print_function, division
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
    
    plt.show()



Using custom objective functions
----------------------------------

By default, funcFit minimizes :math:`\chi^2` when an error is given and the quadratic model deviation otherwise. It may, however,
be necessary to minimize something else such as the likelihood for instance. The following example shows
how to apply a custom objective function, in this case, we simply use the linear deviation between
model and data (weighted by the error) to define the fit quality. 

::

    from __future__ import print_function, division
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
    plt.show()


Some points may be highlighted in this example:
  * You may have noticed that although the parameter `P` is given to the *mini* function, it is not
    used there. You cannot leave it out, however, because the decorator, in fact, creates a more complex
    object, which needs this information.
  * The penalty assignment (for restricted parameters) is done automatically. You do not have
    to include it in your objective function.
  * The custom objective function has to be specified on call to the fit routine (*miniFunc* keyword).


Using an overbinned model 
-----------------------------

In some cases it may be necessary to evaluate a model at more points than actually required by, e.g.,
an observation. The final model is than obtained by averaging a number of points. This may be
necessary to take finite integration times of your instrument into account as can be the case
in planetary transit modeling.

The *turnIntoRebin* method of funcFit provides a convenient way to work with
such "overbinned" models; a demonstration is given in the example below.

::

    from __future__ import print_function, division
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
    plt.show()


This example is very similar to the very first one. Some differences shall, however, be
emphasized:
  * Obtaining the model:
    In this example, we used a model that we called *GaussFit1d_Rebin*. 
    We created the model by calling the *turnIntoRebin* method giving
    GaussFit1d (by name NOT instance, i.e., we use the class object) as the parameter. The return value of this function is another class object,
    in particular, GaussFit1d extended by the overbinning functionality. In the next line,
    we instantiate this extended model and use it, just as we would use the original model.
  * In the end, the overbinned model and the final averaged model are juxtaposed
    to highlight the effect.
 


Fit two models simultaneously
--------------------------------

The following example demonstrates how the *SyncFitContainer* class can be used
to fit two different models with a partly overlapping parameter set, but differing
x-axes simultaneously.
 
::
    
    from __future__ import print_function, division
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
    
    plt.show()

