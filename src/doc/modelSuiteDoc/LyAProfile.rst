Hydrogen Lyman-alpha line-profile
=====================================

.. p23ready

.. currentmodule:: PyAstronomy.modelSuite
.. autoclass:: LyaTransmission
   :members:

Example
--------

::
  
    from PyAstronomy import modelSuite as ms
    import numpy as np
    import matplotlib.pylab as plt
    
    la = ms.LyaTransmission()
    # Set some parameters
    la["N"] = 5e17
    la["b"] = 12.2
    la["Dfrac"] = 1.9e-5
    
    # Set up wavelength axis ...
    wvl = np.linspace(1214.,1217.,1000)
    # ... and evaluate model
    m = la.evaluate(wvl)
    
    # Plot the result
    plt.plot(wvl, m, 'b.-')
    plt.show()