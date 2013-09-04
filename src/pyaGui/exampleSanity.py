import unittest

class ExampleSanity(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  
  def sanity_modelExplorer(self):
    """
      Check sanity of the model explorer example
    """
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
    # The 'version' flag specifies whether a "list" view
    # or a "dropdown" menu is used to specify the currently
    # adapted parameter.
    g = pyaGui.ffmodelExplorer(gg, mp, version="list")
#     g.show()


  def sanity_modelExplorerCustomPlot(self):
    """
      Check custom plot example of model explorer.
    """
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
#     g.show()