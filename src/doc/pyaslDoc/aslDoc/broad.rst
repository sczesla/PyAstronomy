Instrumental (Gaussian kernel) broadening
===========================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: broadGaussFast
.. autofunction:: instrBroadGaussFast

Example of usage
-----------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    import numpy as np
    
    # Set up an input spectrum
    x = np.linspace(5000.0,5100.0,1000)
    y = np.ones(x.size)
    
    # Introduce some delta-peaked lines
    y[165] = 0.7
    y[187] = 0.3
    y[505] = 0.1
    y[610] = 0.1
    y[615] = 0.7
    
    # Apply Gaussian instrumental broadening, setting the resolution to 10000.
    r, fwhm = pyasl.instrBroadGaussFast(x, y, 10000,
              edgeHandling="firstlast", fullout=True)
    
    # Apply Gaussian instrumental broadening, setting the resolution to 10000.
    # Limit the extent of the Gaussian broadening kernel to five standard
    # deviations.
    r2, fwhm = pyasl.instrBroadGaussFast(x, y, 10000,
              edgeHandling="firstlast", fullout=True, maxsig=5.0)
    
    print("FWHM used for the Gaussian kernel: ", fwhm, " A")
    
    # Plot the output
    plt.plot(x,r, 'r--p', label="Broadened curve (full)")
    plt.plot(x, r2, 'k:', label="Broadened curve (5 stds)")
    plt.plot(x,y, 'b-', label="Input")
    plt.legend(loc=4)
    plt.show()