
Analytical transit model (Mandel & Agol 2002)
======================================================

The *forTrans* module provides easy access to the FORTRAN
routines by Mandel & Agol 2002 [#ManAg02]_, which can be used to
calculate transit light curves. While the FORTRAN code ensures
fast computation of the light curve, PyAstronomy supplies
funcFit's fitting framework.

.. note:: The code contained within the files *occultnl.f* and *occultquad.f* is used
          with kind permission of the authors. It is described in [#ManAg02]_; further
          explanations are given in Eastman et al. 2013 (PASP 125, 83).

.. [#ManAg02] Mandel & Agol 2002, "Analytic Light Curves for Planetary Transit Searches", ApJ, 580, 171.

Quadratic limb-darkening - The MandelAgolLC class
-------------------------------------------------

.. currentmodule:: PyAstronomy.modelSuite.XTran.forTrans.mandelAgol
.. autoclass:: MandelAgolLC



Non-linear limb-darkening - The MandelAgolNLLC class
----------------------------------------------------

.. currentmodule:: PyAstronomy.modelSuite.XTran.forTrans.mandelAgolNL
.. autoclass:: MandelAgolNLLC



Example code - Calculate model light curve
------------------------------------------

::

  # Import some unrelated modules
  import numpy
  import matplotlib.pylab as mpl
  
  # ... and now the forTrans module
  from PyAstronomy.modelSuite import forTrans as ft

  # Load a simple timer for demonstration purposes
  import time as timerMod

  # Create MandelAgolLC object
  ma = ft.MandelAgolLC()

  # Set parameters
  ma["per"] = 0.2
  ma["i"] = 90.0
  ma["a"] = 6.5
  ma["T0"] = 0.5
  ma["p"] = 0.16
  ma["linLimb"] = 0.47
  ma["quadLimb"] = 0.24
  ma["b"] = 0.
  
  # Choose some (large) time axis
  time = numpy.linspace(0,5,100000)
  
  # Start the clock
  tstart = timerMod.time()
  
  # ... and calculate model
  y = ma.evaluate(time)
  
  # Stop the clock
  tstop = timerMod.time()
  print "Elapsed time [sec]:", tstop -tstart
  
  # Let's see what happened ...
  mpl.plot(time, y,'.')
  mpl.show()
  
  


Example code - Comparing limb-darkening laws
--------------------------------------------

::

  # Import some modules
  import numpy
  import matplotlib.pylab as mpl
  
  # ... and now the forTrans module
  from PyAstronomy.modelSuite import forTrans as ft
  
  # First, let's compute a transit model using
  # quadratic limb-darkening prescription.
  ma = ft.MandelAgolLC()
  
  # Set parameters. The LD coefficients are taken
  # from Claret 2011 for a solar-metallicity star
  # with Teff=6000 K and logg=4.5.
  ma["per"] = 0.2
  ma["i"] = 90.0
  ma["a"] = 6.5
  ma["T0"] = 0.5
  ma["p"] = 0.16
  ma["linLimb"] = 0.0479
  ma["quadLimb"] = 0.2716
  ma["b"] = 0.

  # Choose some time axis
  time = numpy.linspace(0,0.2,1000)
  
  # ... and calculate model
  yQLD = ma.evaluate(time)
  
  # Now, let's compute a transit model with
  # non-linear limb-darkening prescription
  # for the same stellar parameters.
  maNL = ft.MandelAgolNLLC()
  maNL["per"] = 0.2
  maNL["i"] = 90.0
  maNL["a"] = 6.5
  maNL["T0"] = 0.5
  maNL["p"] = 0.16
  maNL["a1"] = 0.5335
  maNL["a2"] = 0.0793
  maNL["a3"] = -0.3466
  maNL["a4"] = 0.1609
  maNL["b"] = 0.
  
  yNLLD = maNL.evaluate(time)
  
  # Let's compare both models...
  mpl.plot(time, yQLD,'-',label="Quadratic LD")
  mpl.plot(time, yNLLD,'d',label="Non-linear LD")
  mpl.legend()
  mpl.show()
  
  

