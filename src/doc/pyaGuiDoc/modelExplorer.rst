Interactive model explorer
=================================

The "FuncFit model explorer" is a GUI tool, which can be used
to manipulate and plot `funcFit` models interactively.
Individual model parameters can be selected, and their values can
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
  mp = pyaGui.FFModelPlotFit(x, y, yerr=yerr, withResiduals=False)
  
  # Create an instance of the FFModelExplorer class, which
  # needs to be given the model (gg in this case) and
  # the plotter (and fitter), which we created above.
  g = pyaGui.FFModelExplorer(gg, mp)
  g.show()


.. note:: Plotting and fitting is not done by the class `FFModelPlotFit`, which
          provides the basic capabilities to handle one-dimensional models.
          In more general cases, custom classes can be used. 


Example---Providing a custom class for plotting and fitting
------------------------------------------------------------

A *custom* class for plotting and fitting must provide two methods:
  1) `plot`
  2) `fit`.

The *plot* method takes two arguments: The first is an instance of
a matplotlib Figure class on which drawing takes place and the second
is the model class under consideration.

The *fit* class takes a single argument, which is the model class instance.

In the example below, we demonstrate how such a class can be implemented.
In particular, the class *CustomPlotFit* uses an arbitrary number of panels
to show the plot.

::

  import numpy as np
  from PyAstronomy import pyaGui
  from PyAstronomy import funcFit as fuf


  class CustomPlotFit:
    """
      A custom class for plotting and fitting.
      
      This class uses several panels
      can be used for plotting.
      
      Parameters
      ----------
      x, y : arrays
          The data.
      n : int
          The number of panels.
    """
    
    def __init__(self, x, y, n):
      self.x = x
      self.y = y
      self.n = n
      self.a = None
      self.color = {0:"b", 1:"r", 2:"g"}
      
    def plot(self, f, odf):
      """
        Plot the model into n panels.
        
        Parameters
        ----------
        f : Matplotlib Figure
            Figure on which to plot.
        odf : Model class instance
            The model to be plotted.
      """
      
      if self.a is None:
        # Create axis instances for plotting
        self.a = []
        for i in range(self.n):
          self.a.append(f.add_subplot(self.n, 1, i+1))
      
      model = odf.evaluate(self.x)
      for i in range(self.n):
        self.a[i].cla()
        self.a[i].plot(self.x, self.y, self.color[i % 3]+'p')
        self.a[i].plot(self.x, model, 'k--')
    
    def fit(self, odf):
      """
        Fit the model.
        
        Parameters
        ----------
        odf : Model class instance
            The model to be fitted.
      """
      odf.fit(self.x, self.y)
  
  
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
  
  # Use the custom class for plotting and fitting, which
  # will show the plot in 3 panels.
  cm = CustomPlotFit(x, y, 3)
  
  # Create an instance of the FFModelExplorer class, which
  # needs to be given the model (gg in this case) and
  # the plotter (and fitter), which we created above.
  g = pyaGui.FFModelExplorer(gg, cm)
  g.show()


Implementation
----------------

.. currentmodule:: PyAstronomy.pyaGui
.. autoclass:: FFModelExplorer
   :members:

.. autoclass:: FFModelPlotFit
   :members: