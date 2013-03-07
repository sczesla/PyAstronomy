Access the exoplanet.eu data base
==================================

.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: ExoplanetEU
   :members:
   :inherited-members:

Example of usage
-----------------

::

  from PyAstronomy import pyasl
  import matplotlib.pylab as plt
  
  eu = pyasl.ExoplanetEU()
  
  # See what information is available
  cols = eu.availableColumns()
  print cols
  
  print
  # Get all data and plot planet Mass vs.
  # semi-major axis in log-log plot
  dat = eu.getAllData()
  plt.xlabel("Planet Mass [RJ]")
  plt.ylabel("Semi-major axis [AU]")
  plt.loglog(dat.plMass, dat.sma, 'b.')
  plt.show()