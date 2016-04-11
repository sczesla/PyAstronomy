Access to Kurucz atmospheric models
======================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

The classes and functions available here provide access
to the model grids made available by Robert L. Kurucz
on: http://kurucz.harvard.edu/grids.html

Example: Get access to the models
---------------------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    km = pyasl.KuruczModels()
    # See what model grids are available
    print(km.availableGrids())
    
    # See whether model grid for log(metallicity) = 0.0
    # is available
    print(km.gridAvailable(0.0))
    
    # Obtain the model grid for solar metallicity
    mg = km.requestModelGrid(0.0)
    
    # See what Teffs and logg are available
    print("Teffs: ", mg.availableTeffs())
    print("Loggs: ", mg.availableLoggs())
    
    print()
    print()
    
    # Use simple access method to obtain a model.
    # The input is: Teff, logg, and log10(metallicity)
    model = pyasl.getKuruczModel(4250, 4.5, 0.1)


Classes and functions
-----------------------

.. autofunction:: getKuruczModel
.. autoclass:: KuruczModels
   :members:
.. autoclass:: KuruczMT
   :members: