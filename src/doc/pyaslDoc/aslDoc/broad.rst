Instrumental (Gaussian kernel) broadening
===========================================

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: broadGaussFast
.. autofunction:: instrBroadGaussFast

Example of usage
-----------------

::

  from PyAstronomy import pyasl
  import matplotlib.pylab as plt
  import numpy as np
  
  # Set up an input spectrum
  x = np.linspace(5000.0,5100.0,1000)
  y = np.ones(x.size)
  
  # Introduce some delta-peaked lines
  y[165] = 0.7
  y[187] = 0.3
  y[505] = 0.1
  y[610] = 0.1
  y[615] = 0.7
  
  # Apply Gaussian instrumental broadening, setting the resolution to 10000.
  r, fwhm = pyasl.instrBroadGaussFast(x, y, 10000,
            edgeHandling="firstlast", fullout=True)
  
  print "FWHM used for the Gaussian kernel: ", fwhm, " A"
  
  # Plot the output
  plt.plot(x,r, 'r--p')
  plt.plot(x,y, 'b-')
  plt.show()