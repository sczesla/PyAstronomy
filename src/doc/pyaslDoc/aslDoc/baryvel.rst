Barycentric velocity correction
=================================

.. currentmodule:: PyAstronomy.pyasl

The :py:func:`baryvel` and :py:func:`baryCorr` function allow to
calculate Earth's helio- and barycentric motion and project it
onto a given direction toward a star.

PyA's `baryvel` function is a port of its IDL Astrolib's counterpart. 

.. autofunction:: baryvel
.. autofunction:: baryCorr

Example
--------

::

  from PyAstronomy import pyasl
  
  jd = 2.476468576e6
  
  heli, bary = pyasl.baryvel(jd, deq=2000.0)
  
  print "Earth's velocity at JD: ", jd
  print "Heliocentric velocity [km/s]: ", heli
  print "Barycentric velocity [km/s] : ", bary
  
  # Coordinates of Sirius
  ra  = 101.28715535
  dec = -16.71611587
  
  vh, vb = pyasl.baryCorr(jd, ra, dec, deq=2000.0)
  print "Barycentric velocity of Earth toward Sirius: ", vb
