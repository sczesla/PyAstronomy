Interactive data normalization and continuum finder
====================================================

.. p23ready

The `ContinuumInteractive` class provides a GUI, which helps to
interactively normalize data, e.g., for continuum normalization.

Example of usage
-------------------

::

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
  c = cf.findContinuum()
  
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
  plt.title("Normalized data")
  plt.plot(x, c["normalizedData"], 'b.--')
  plt.show()

Class documentation
-----------------------

.. currentmodule:: PyAstronomy.pyaGui
.. autoclass:: ContinuumInteractive
   :members: