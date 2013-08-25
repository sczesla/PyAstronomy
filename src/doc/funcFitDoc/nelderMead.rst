Nelder-Mead downhill simplex
==============================

.. currentmodule:: PyAstronomy.funcFit
.. autoclass:: NelderMead
   :members:
   :private-members:

Example: Application of the NelderMead class
----------------------------------------------

::

  from PyAstronomy import funcFit as fuf
  import numpy as np
  
  # Construct the fitter
  nm = fuf.NelderMead()
  
  # Get a GaussFit object ...
  gf = fuf.GaussFit1d()
  gf.thaw(["A", "mu"])
  # ... and define some initial values.
  gf["A"] = 5.0
  gf["mu"] = 2.0
  gf["sig"] = 1.0
  
  # Construct some "data" ...
  x = np.arange(100)/20.
  y = gf.evaluate(x) + np.random.normal(0.0, 0.08, len(x))
  yerr = np.ones(len(x)) * 0.08
  # ... and define a "data set" object
  ds = fuf.FufDS(x, y, yerr=yerr)
  
  # On purpose, we vail our knowledge about the correct
  # parameters
  gf["A"] = 1.0
  gf["mu"] = 1.0
  
  # We carry out the fit. In particular, the squared
  # distance between the model and the data is minimized
  # (sqrdiff).
  bfpars = nm.fit(gf, ds, objf="sqrdiff", maxIter=100)
  # Have a look at the result
  print "Best-fit parameters: ", bfpars
  print
  # Have a look at the entire fit object
  gf.parameterSummary()
  
  
  print
  # Again, we change the starting parameters. Note that ...
  gf["A"] = 1.0
  # ... mu has been set to zero.
  gf["mu"] = 0.0
  
  # To allow the construction of an appropriate start simplex,
  # we need to provide an "initDelta" for mu, which should
  # reflect the typical scale of the problem.
  bfpars = nm.fit(gf, ds, objf="chisqr", initDelta={"mu":0.1})
  print "Second fit, best-fit parameters: ", bfpars