X-ray luminosity / rotation period relations
=============================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: Pizzolato2003
   :members:

Example: Plot Lx and Lx/Lbol as a function of rotation period
--------------------------------------------------------------

::

  from PyAstronomy import pyasl
  import matplotlib.pylab as plt
  import numpy as np
  
  p = pyasl.Pizzolato2003()
  
  # Define array of rotation periods [days]
  prot = np.arange(0.2, 30, 0.1)
  
  lx = np.zeros(prot.size)
  lxlbol = np.zeros(prot.size)
  
  # B-V color of star
  bv = 0.7
  
  # Obtain ...
  for i in range(prot.size):
    # ... log10 of X-ray luminosity
    lx[i] = p.log10lxbv(bv, prot[i])[0]
    # ... and log10(Lx/Lbol)
    lxlbol[i] = p.log10lxlbolbv(bv, prot[i])[0]
  
  # Plot result
  plt.subplot(2,1,1)
  plt.plot(prot, lx, 'bp-')
  plt.subplot(2,1,2)
  plt.plot(prot, lxlbol, 'bp-')
  plt.show()
