INTEP interpolation algorithm
===============================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: intep

Example
----------

::

  from PyAstronomy import pyasl
  import numpy as np
  import matplotlib.pylab as plt
  
  # Create some tabulated data
  x = np.arange(10.)
  y = np.sin(x/5.*2.*np.pi)
  
  # Choose the values at which to interpolate
  xx = np.arange(120.)/10.-2.0
  # Interpolate and suppress the exception that
  # would indicate that some of our requested
  # values (xx) are beyond the range covered by x.
  yy = pyasl.intep(x, y, xx, boundsError=False)
  
  # Plot the result
  plt.plot(x, y, 'bp')
  plt.plot(xx, yy, 'r--')
  plt.show()