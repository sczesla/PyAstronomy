Deredden flux vector (unred)
=============================

.. p23ready

.. currentModule:: PyAstronomy.pyasl

.. autofunction:: unred

Example of usage
------------------

::
  
    from PyAstronomy import pyasl
    import numpy as np
    import matplotlib.pylab as plt
    
    # Approximate a solar spectrum using a Planck
    # function with a temperature of 5778 K between
    # 3000 A and 8000 A.
    wvl = np.arange(3000., 8000., 1.0)
    flux = pyasl.planck(T=5778., lam=wvl*1e-10)
    
    # Deredden the spectrum assuming ebv=0.1
    fluxUnred = pyasl.unred(wvl, flux, ebv=0.1, R_V=3.1)
    
    # Plot the result
    plt.title("Reddened flux (red) and dereddened flux (blue)")
    plt.plot(wvl, flux, 'r--')
    plt.plot(wvl, fluxUnred, 'b--')
    plt.show()