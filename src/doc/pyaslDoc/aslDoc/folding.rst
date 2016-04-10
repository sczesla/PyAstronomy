Folding time series
======================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: foldAt

Folding example
------------------

::
  
  from PyAstronomy.pyasl import foldAt
  import matplotlib.pylab as plt
  import numpy as np
  
  # Generate some data ...
  time = np.random.random(1000) * 100.
  flux = 0.05 * np.sin(time*(2.*np.pi/21.5) + 15)
  # ... and add some noise
  flux += np.random.normal(0, 0.02, len(flux))
  
  # Obtain the phases with respect to some
  # reference point (in this case T0=217.4)
  phases = foldAt(time, 21.5, T0=217.4)
  
  # Sort with respect to phase
  # First, get the order of indices ...
  sortIndi = np.argsort(phases)
  # ... and, second, rearrange the arrays.
  phases = phases[sortIndi]
  flux = flux[sortIndi]
  
  # Plot the result
  plt.plot(phases, flux, 'bp')
  plt.show()