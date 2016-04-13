from __future__ import print_function, division
import unittest

class ExampleSanity(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_pyaPickerExample(self):
    from PyAstronomy import pyaGui
    import numpy as np
    
    # Data for the plot
    x = np.linspace(0., 10., 100.)
    y = np.exp(-x/10.)
    
    # Create Picker instance
    pp = pyaGui.Picker()
    
#     # Create the plot in which points
#     # are to be selected
#     pp.a.plot(x, y, 'b.-')
#     
#     points = pp.pick()
#     
#     print("Number of selected points: ", len(points))
#     print("Selected points (x, y):")
#     for p in points:
#       print("%g, %g" % p) 
  
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
    g = pyaGui.ffmodelExplorer(gg, mp)
#    g.show()


#  def sanity_modelExplorerCustomPlot(self):
#    """
#      Check custom plot example of model explorer.
#    """
#    import numpy as np
#    from PyAstronomy import pyaGui
#    from PyAstronomy import funcFit as fuf
#    
#    
#    class CustomPlotFit:
#      """
#        A custom class for plotting and fitting.
#        
#        This class uses several panels
#        can be used for plotting.
#        
#        Parameters
#        ----------
#        x, y : arrays
#            The data.
#        n : int
#            The number of panels.
#      """
#      
#      def __init__(self, x, y, n):
#        self.x = x
#        self.y = y
#        self.n = n
#        self.a = None
#        self.color = {0:"b", 1:"r", 2:"g"}
#        
#      def plot(self, f, odf):
#        """
#          Plot the model into n panels.
#          
#          Parameters
#          ----------
#          f : Matplotlib Figure
#              Figure on which to plot.
#          odf : Model class instance
#              The model to be plotted.
#        """
#        
#        if self.a is None:
#          # Create axis instances for plotting
#          self.a = []
#          for i in range(self.n):
#            self.a.append(f.add_subplot(self.n, 1, i+1))
#        
#        model = odf.evaluate(self.x)
#        for i in range(self.n):
#          self.a[i].cla()
#          self.a[i].plot(self.x, self.y, self.color[i % 3]+'p')
#          self.a[i].plot(self.x, model, 'k--')
#      
#      def fit(self, odf):
#        """
#          Fit the model.
#          
#          Parameters
#          ----------
#          odf : Model class instance
#              The model to be fitted.
#        """
#        odf.fit(self.x, self.y)
#    
#    
#    # Create a Gaussian fitting instance
#    # and set some parameters
#    gg = fuf.GaussFit1d()
#    gg["A"] = 1.0
#    gg["sig"] = 0.5
#    
#    # Let A and mu be free during a fit
#    gg.thaw(["A", "mu"])
#    
#    # Create some artificial data
#    x = np.linspace(-2., 2., 100)
#    yerr = np.ones(len(x))*0.01
#    y = gg.evaluate(x) + np.random.normal(0.,0.01, len(x))
#    
#    # Use the custom class for plotting and fitting, which
#    # will show the plot in 3 panels.
#    cm = CustomPlotFit(x, y, 3)
#    
#    # Create an instance of the FFModelExplorer class, which
#    # needs to be given the model (gg in this case) and
#    # the plotter (and fitter), which we created above.
#    g = pyaGui.FFModelExplorer(gg, cm)
##     g.show()

  def sanity_InteractiveContinuum(self):
    """
      Checking sanity of interactive continuum.
    """
    from PyAstronomy import pyaGui as pg
    import numpy as np
    import matplotlib.pylab as plt
    
    # Create some artificial data/spectrum
    x = np.arange(1000)
    y = np.random.random(1000) + 1000.0 + np.sin(x/80.)
    y += np.random.random(len(x))
    
    # Construct an instance of the interactive normalization GUI
    cf = pg.ContinuumInteractive(x, y)
    
    # The 'plot' command can be used to draw something on the canvas
    # on which normalization points are selected. Uncomment the following
    # line to see an example.
    #
    # cf.plot([4,400,600], [1001,1000,999], 'k--')
    
    # The 'plotNorm' command can be used to draw something on the
    # canvas showing the normalized spectrum (if display is enabled).
    # Uncomment the following line to see an example.
    #
    # cf.plotNorm([0.,1000.], [1.001,0.999], 'g:')
    
    # Opens the GUI and starts the interactive session.
#    c = cf.findContinuum()
    
    # The outcome is a dictionary with the following keys:
    #
    #            - points: A list of two-float tuples holding the
    #                      x,y location of the selected points.
    #            - continuum : Array holding the continuum estimate
    #                          at the given x-values.
    #            - splineKind : A string specifying the selected
    #                           spline option.
    #            - normalizedData : An array holding the normalized
    #                               data.
#    plt.title("Normalized data")
#    plt.plot(x, c["normalizedData"], 'b.--')
#    plt.show()

  def sanity_IAGVFitExample(self):
    """
      Checking example of IAGVFit
    """
    from PyAstronomy import pyaGui
    from PyAstronomy import funcFit as fuf
    import numpy as np
    
    # Data for the plot
    x = np.linspace(5000., 5010, 200)
    y = np.ones(len(x))
    yerr = np.ones(len(x)) * 0.01
    y += np.random.normal(0., 0.01, len(x))
    
    gf = fuf.GaussFit1d()
    gf["A"] = -0.3
    gf["mu"] = 5004.
    gf["sig"] = 0.2
    y += gf.evaluate(x)
    
    # Create interactive fitter
    igv = pyaGui.IAGVFit(x, y, yerr=yerr, mode="gauss")
    
#     r = igv.interactiveFit()
    
#     print("Parameters of the fit: ", r[0])
#     print("Parameters of active component: ", r[1])
#     print("No. of components: ", r[2])
