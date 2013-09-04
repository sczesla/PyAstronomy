Access the exoplanets.org data base
====================================

.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: ExoplanetsOrg
   :members:
   :inherited-members:

Example of usage
-----------------

::

  from PyAstronomy import pyasl
  
  # Instantiate the access class
  epl = pyasl.ExoplanetsOrg()
  
  # Show the available columns
  epl.availableColumns()
  
  # Get information in Kepler-5 b
  d = epl.selectByPlanetName("kepler-5 b")
  
  # Print whatever information has been received
  print
  print "Information on Kepler-5 b"
  print
  for k, v in d.items():
    print "%12s  %12s" % (k,str(v))