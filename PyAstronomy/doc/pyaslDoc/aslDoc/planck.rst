Planck's radiation law
=======================

.. currentModule:: PyAstronomy.pyasl

.. autofunction:: planck


Example
-----------

::
    
    from __future__ import print_function
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy.pyasl import planck
    
    # Define wavelength in meters
    lam = np.arange(1000.0*1e-10, 20000.*1e-10, 20e-10)
    
    # Get the Planck spectrum in [W/(m**2 m)] for a temperature of 7000 K
    s7 = planck(7000., lam=lam)
    # Get the Planck spectrum in [W/(m**2 m)] for a temperature of 5000 K
    s5 = planck(5000., lam=lam)
    
    # Convert into erg/(cm**2 * A * s)
    s5erg = s5 * 1e-7
    s7erg = s7 * 1e-7
    
    # Integrate the spectrum and compare with Stefan-Boltzmann law
    i5 = np.sum(s5) * (lam[1] - lam[0])
    i7 = np.sum(s7) * (lam[1] - lam[0])
    
    print("5000 K integral: %.3e W/m**2 (Stefan-Boltzmann predicts %.3e W/m**2)" % (i5, (5.67e-8*5000.**4)))
    print("7000 K integral: %.3e W/m**2 (Stefan-Boltzmann predicts %.3e W/m**2)" % (i7, (5.67e-8*7000.**4)))
    
    plt.xlabel("Wavelength [$\AA$]")
    plt.ylabel("Flux [erg/cm$^2$/A/s]")
    plt.plot(lam*1e10, s5erg, 'r-')
    plt.plot(lam*1e10, s7erg, 'b-')
    plt.show()
 