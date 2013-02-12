Observatory locations
=======================

.. currentModule:: PyAstronomy.pyasl

.. autofunction::  listObservatories
.. autofunction::  observatory

Example of usage
-----------------

::

  from PyAstronomy import pyasl
  
  # List all available observatory data
  pyasl.listObservatories()
  
  print
  print "Data for Kitt Peak National Observatory"
  print pyasl.observatory("kpno")
  print "(longitude and latitude in degrees, altitude in meters, and"
  print "time zone in hours West of Greenwich"