Access the "NASA Exoplanet Archive"
=====================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: NasaExoplanetArchive
   :members:
   :inherited-members:

Example of usage
------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    
    nexa = pyasl.NasaExoplanetArchive()
    
    # See what information is available
    cols = nexa.availableColumns()
    print()
    
    # Get all information for planet 'wasp-12 b'
    # By default, the search is case-insensitive
    print("Entry of Wasp-12 b")
    print(nexa.selectByPlanetName("Wasp-12 b"))
    
    print()
    # Get all data and plot ra vs. dec
    dat = nexa.getAllData()
    plt.plot(dat.ra, dat.dec, 'b.')
    plt.show()
