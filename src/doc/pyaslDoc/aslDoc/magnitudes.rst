Magnitude conversions
=======================

.. currentModule:: PyAstronomy.pyasl

Translate absolute magnitude to power scale
--------------------------------------------

.. autofunction:: absMagToPower

Example
~~~~~~~~

::

  from PyAstronomy import pyasl
  
  absMagSun = 4.75
  print "Absolute bolometric magnitude of the Sun: ", absMagSun
  print "  Absolute luminosity [erg/s]: ", pyasl.absMagToPower(absMagSun)