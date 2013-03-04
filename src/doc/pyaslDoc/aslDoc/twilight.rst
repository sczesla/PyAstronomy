Name twilight given solar altitude
====================================

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: twilightName

Example
--------

::

  from PyAstronomy import pyasl
  import numpy as np
  
  for alt in np.linspace(-20., 5., 15):
    print "Altitude = {0:6.2f}, Twilight is called: ".format(alt), \
          pyasl.twilightName(alt)


