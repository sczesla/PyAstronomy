Voigt profile with astronomical parameterization
==================================================

.. p23ready

.. currentmodule:: PyAstronomy.modelSuite

.. autoclass:: VoigtAstroP
   :members:

Example
~~~~~~~~~

::
  
    from PyAstronomy import modelSuite as ms
    import numpy as np
    import matplotlib.pylab as plt
    
    # Obtain an object of type VoigtAstroP ...
    v = ms.VoigtAstroP()
    # ... and set some parameters
    v["b"] = 87.7
    v["f"] = 0.5
    v["w0"] = 1214.0
    # Damping constant [cm]
    v["gamma"] = 2e-9
    
    # Generate wavelength axis ...
    wvl = np.linspace(1212.,1216.,200)
    # ... and evaluate model
    m = v.evaluate(wvl)
    
    # Plot result
    plt.plot(wvl, m, 'b.-')
    plt.show()
