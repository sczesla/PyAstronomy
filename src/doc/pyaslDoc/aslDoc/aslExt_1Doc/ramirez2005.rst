Converting between effective temperature and stellar color
===============================================================

In their 2005 publication, Ramírez and Meléndez (ApJ 626, 465-485) present
metallicity-dependent relations between stellar effective temperature
and color. Based on these relations, the class Ramirez2005 allows
to convert between effective temperature and color. All 17 color
indices given by the authors can be used.

Example: 
---------------------------

::

  from PyAstronomy import pyasl

  # Create class instance
  r = pyasl.Ramirez2005()

  # Which color bands are available
  print "Available color bands: ", r.availableBands()
  
  # Convert B-V to effective temperature and back
  bv = 0.75
  feh = 0.0
  teff = r.colorToTeff("B-V", bv, feh)
  bv1 = r.teffToColor("B-V", teff, feh)
  # Watch out for differences between input bv and the output bv1
  print "B-V = ", bv, ", Teff = ", teff, ", bv1 = ", bv1, ", bv-bv1 = ", bv-bv1
  
API documentation
--------------------
  
.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: Ramirez2005
   :members:
   :private-members: