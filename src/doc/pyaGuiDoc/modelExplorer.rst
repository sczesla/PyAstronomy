Interactive model explorer
=================================

.. p23ready

The "FuncFit model explorer" is a GUI tool, which can be used
to manipulate and plot `funcFit` models interactively.
Individual model parameters can be selected and their values can
be manipulated via the *mouse wheel* or changed directly.

See the following example to explore its possibilities.

Example---Manipulating a Gaussian
------------------------------------ 

This example demonstrates the basic usage of the model explorer.

::

  import numpy as np
  from PyAstronomy import pyaGui
  from PyAstronomy import funcFit as fuf
  
  # Create a Gaussian fitting instance
  # and set some parameters
  gg = fuf.GaussFit1d()
  gg["A"] = 1.0
  gg["sig"] = 0.5
  
  # Let A and mu be free during a fit
  gg.thaw(["A", "mu"])
  
  # Create some artificial data
  x = np.linspace(-2., 2., 100)
  yerr = np.ones(len(x))*0.01
  y = gg.evaluate(x) + np.random.normal(0.,0.01, len(x))
  
  # In order to use the interactive explorer, you
  # need a class, which plots the model. The default
  # class for this purpose is "FFModelPlotFit", which
  # needs the x and y values. Optionally, you can specify
  # errors via `yerr`. Depending on the setting for
  # "withResiduals", the residuals will be shown or not.
  mp = pyaGui.FFModelPlotFit(x, y, yerr=yerr, withResiduals=True)
  
  # Use the function ffmodelExplorer (note the lowercase letters)
  # to create an instance of the FFModelExplorer class, which
  # needs to be given the model (gg in this case) and
  # the plotter (and fitter), which we created above.
  # 
  g = pyaGui.ffmodelExplorer(gg, mp)
  g.show()



Implementation
----------------

.. currentmodule:: PyAstronomy.pyaGui
.. autoclass:: FFModelExplorerList
   :members:

.. autoclass:: FFModelExplorerDropDownMenu
   :members:

.. autoclass:: FFModelPlotFit
   :members: