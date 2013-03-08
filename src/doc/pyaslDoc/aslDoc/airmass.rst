Calculate the airmass
======================

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: airmassSpherical

Example: Airmass with homogeneous, spherical atmosphere
---------------------------------------------------------

::

  from PyAstronomy import pyasl
  
  obsAlt = 2400.0
  
  for za in range(0,90,10):
    print "Zenith angle: %2d, airmass = %7.2f" % \
      (za, pyasl.airmassSpherical(za, obsAlt))