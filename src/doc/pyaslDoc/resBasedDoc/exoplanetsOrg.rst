Access the exoplanets.org data base
====================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: ExoplanetsOrg
   :members:
   :inherited-members:

Example of usage
-----------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Instantiate the access class
    epl = pyasl.ExoplanetsOrg()
    
    # Show the available columns
    epl.availableColumns()
    
    # Get information in Kepler-5 b
    d = epl.selectByPlanetName("kepler-5 b")
    
    # Print whatever information has been received
    print()
    print("Information on Kepler-5 b")
    print()
    for k, v in list(d.items()):
      print("%12s  %12s" % (k,str(v)))
      