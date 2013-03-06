Estimate the transit duration
==============================

The following function estimates the duration
of a transit, i.e., the time from first to fourth
contact assuming spherical objects and a circular
orbit.

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: transitDuration

Example: Estimate duration of Earth's transit
-----------------------------------------------

::

  from PyAstronomy import pyasl
  from PyAstronomy import constants as pc
  
  # Earth radius expressed in Jovian radii
  reJ = pc.REarth/pc.RJ
  print "Earth radius in Jovian units: ", reJ
  
  # Estimate the duration of Earth's transit
  td = pyasl.transitDuration(1.0, reJ, 1.0, 90.0, 365.0)
  print "The transit of Earth lasts about: %5.3f days" % td