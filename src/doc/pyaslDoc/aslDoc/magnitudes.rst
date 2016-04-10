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
