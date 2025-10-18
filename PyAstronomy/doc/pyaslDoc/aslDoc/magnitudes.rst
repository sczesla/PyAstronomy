Magnitude conversions
=======================

.. p23ready

.. currentModule:: PyAstronomy.pyasl

Translate absolute magnitude to power scale
--------------------------------------------

.. autofunction:: absMagToPower

Example
~~~~~~~~

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    absMagSun = 4.75
    print("Absolute bolometric magnitude of the Sun: ", absMagSun)
    print("  Absolute luminosity [erg/s]: ", pyasl.absMagToPower(absMagSun))


Translate distance module into distance
----------------------------------------

.. autofunction:: absModuleToDist

Example
~~~~~~~~

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Apparent magnitude
    appMag = 11.37
    # Absolute (bolometric) magnitude of Sun
    absMagSun = 4.75
    
    print("Distance of a sun-like star with apparent bolometric ", end=' ')
    print("brightness of 11.37 mag: %5.2f pc" % (pyasl.absModuleToDist(appMag, absMagSun)))
    
    
Convert magnitude into flux density
-------------------------------------

.. autofunction:: magToFluxDensity_bessel98


Example
~~~~~~~~

::
    
    from __future__ import print_function
    from PyAstronomy import pyasl
    import numpy as np
    
    mag_R = 15.5
    
    fd_nu, le = pyasl.magToFluxDensity_bessel98("R", mag_R, "nu")
    fd_lam, _ = pyasl.magToFluxDensity_bessel98("R", mag_R, "lam")
    
    print("R-band magnitude: ", mag_R)
    print("R-band flux density [erg/cm**2/s/Hz]: ", fd_nu)
    print("R-band flux density [erg/cm**2/s/A]: ", fd_lam)
    
    print("Effective wavelength of filter [A]: ", le)
    print("Convert f_nu into f_lam [erg/cm**2/s/A] by multiplication with (c/lam**2): ",
          fd_nu * (299792458e2/(le/1e8)**2) / 1e8)

