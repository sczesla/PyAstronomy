PyA Picker --- Select points from a plot
===============================================

.. p23ready

The "Picker" is a simple tool to select points in a
plot interactively.

.. currentmodule:: PyAstronomy.pyaGui
.. autoclass:: Picker
   :members:

Example use of Picker
----------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyaGui
    import numpy as np
    
    # Data for the plot
    x = np.linspace(0., 10., 100.)
    y = np.exp(-x/10.)
    
    # Create Picker instance
    pp = pyaGui.Picker()
    
    # Create the plot in which points
    # are to be selected
    pp.a.plot(x, y, 'b.-')
    
    points = pp.pick()
    
    print("Number of selected points: ", len(points))
    print("Selected points (x, y):")
    for p in points:
      print("%g, %g" % p) 
