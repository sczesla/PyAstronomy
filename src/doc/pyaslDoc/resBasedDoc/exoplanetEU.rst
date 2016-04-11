Access the exoplanet.eu data base
==================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

PyA provides the :py:func:`ExoplanetEU2` to access the data provided by
exoplanet.eu. The class download the data as a Virtual Observatory (VO)
table to provide access to it.

.. note:: PyA also provides the :py:func:`ExoplanetEU` (note the missing '2'), which
          is the predecessor of the above implementation. Although functional,
          this implementation should
          be considered deprecated.

ExoplanetEU2
------------------

Example: Usage of ExoplanetEU2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    
    # Instantiate exoplanetEU2 object
    v = pyasl.ExoplanetEU2()
    
    # Show the available data
    v.showAvailableData()
    print()
    
    # Get a list of all available column names
    acs = v.getColnames()
    print("Available column names: " + ", ".join(acs))
    print()
    
    # Select data by planet name (returns a dictionary)
    print(v.selectByPlanetName("CoRoT-2 b"))
    print()
    
    # Get all data as an astropy table
    at = v.getAllDataAPT()
    
    # Export all data as a pandas DataFrame
    pd = v.getAllDataPandas()
    
    # Plot mass vs. SMA
    plt.title("Mass vs. SMA")
    plt.xlabel("[" + v.getUnitOf("mass") + "]")
    plt.ylabel("[" + v.getUnitOf("semi_major_axis") + "]")
    plt.loglog(at["mass"], at["semi_major_axis"], 'b.')
    plt.show()
    


API documentation (ExoplanetEU2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: ExoplanetEU2
   :members:
   :inherited-members:



   
ExoplanetEU (deprecated)
-----------------------------

Example: Using ExoplanetEU
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    
    eu = pyasl.ExoplanetEU()
    
    # See what information is available
    cols = eu.availableColumns()
    print(cols)
    
    print()
    # Get all data and plot planet Mass vs.
    # semi-major axis in log-log plot
    dat = eu.getAllData()
    plt.xlabel("Planet Mass [RJ]")
    plt.ylabel("Semi-major axis [AU]")
    plt.loglog(dat.plMass, dat.sma, 'b.')
    plt.show()


API documentation (ExoplanetEU)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: ExoplanetEU
   :members:
   :inherited-members:
