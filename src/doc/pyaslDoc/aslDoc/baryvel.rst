Barycentric velocity correction
=================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

The :py:func:`baryvel` and :py:func:`baryCorr` functions allow to
calculate Earth's helio- and barycentric motion and project it
onto a given direction toward a star. The :py:func:`helcorr`
includes a barycentric correction including the effect caused
by the rotating Earth.

PyA's `baryvel` function is a port of its IDL Astrolib's counterpart.
The `helcorr` function has been ported from the REDUCE package
(see Piskunov & Valenti 2002, A&A 385, 1095).

.. autofunction:: baryvel
.. autofunction:: baryCorr
.. autofunction:: helcorr

Example: Applying baryvel and baryCorr
-----------------------------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    jd = 2.476468576e6
    
    heli, bary = pyasl.baryvel(jd, deq=2000.0)
    
    print("Earth's velocity at JD: ", jd)
    print("Heliocentric velocity [km/s]: ", heli)
    print("Barycentric velocity [km/s] : ", bary)
    
    # Coordinates of Sirius
    ra  = 101.28715535
    dec = -16.71611587
    
    vh, vb = pyasl.baryCorr(jd, ra, dec, deq=2000.0)
    print("Barycentric velocity of Earth toward Sirius: ", vb)


Example: Obtaining a barycentric correction
---------------------------------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Coordinates of European Southern Observatory
    # (Coordinates of UT1)
    longitude = 289.5967661
    latitude = -24.62586583
    altitude = 2635.43
    
    # Coordinates of HD 12345 (J2000)
    ra2000 = 030.20313477
    dec2000 = -12.87498346
    
    # (Mid-)Time of observation
    jd = 2450528.2335
    
    # Calculate barycentric correction (debug=True show
    # various intermediate results)
    corr, hjd = pyasl.helcorr(longitude, latitude, altitude, \
                ra2000, dec2000, jd, debug=True)
    
    print("Barycentric correction [km/s]: ", corr)
    print("Heliocentric Julian day: ", hjd)



