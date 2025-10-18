Transit duration and contact point times
===========================================================

.. p23ready
.. currentmodule:: PyAstronomy.pyasl


.. contents::


Transit and in- and egress duration for circular orbit
--------------------------------------------------------------

Functions to estimate the duration of the transit and the in- and egress for circular orbits based on
input in solar system units (AU, solar and Jovian radius)

* :py:func:`Transit duration (solar system units) <transitDuration>`
* :py:func:`Ingress duration (solar system units) <ingressDuration>`

or parameters expressed in terms of the stellar radius

* :py:func:`Transit duration (stellar units) <transitDuration_Rs>`
* :py:func:`Ingress duration (stellar units) <ingressDuration_Rs>`



Example: Estimate duration of Earth's transit and ingress (solar system units)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    from __future__ import print_function, division
    from PyAstronomy import pyasl
    from PyAstronomy import constants as pc
    
    # Earth radius expressed in Jovian radii
    reJ = pc.REarth/pc.RJ
    print("Earth radius in Jovian units: ", reJ)
    
    # Estimate the duration of Earth's transit
    td = pyasl.transitDuration(1.0, reJ, 1.0, 90.0, 365.0)
    print("The transit of Earth lasts about: %5.3f days" % td)

    # Estimate the duration of in- and egress
    ti = pyasl.ingressDuration(1.0, reJ, 1.0, 90.0, 365.0)
    print("The in- and egress last: %6.4f days or %4.2f hours" % (ti, ti*24))


Example: Estimate transit and ingress duration of HD 189733 b (stellar units)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Semi-major axis in units of stellar radius
    sma = 8.8
    # Radius ratio (Rp/Rs)
    rprs = 0.16
    
    # Estimate the duration of Earth's transit
    td = pyasl.transitDuration_Rs(sma, rprs, 85.7, 2.21858)
    print("The transit of HD 189733 b lasts about: %5.2f hours" % (td*24.))

    # Estimate the duration of in- and egress
    ti = pyasl.ingressDuration_Rs(sma, rprs, 85.7, 2.21858)
    print("The in- and egress of HD 189733 b lasts: %5.2f hours" % (ti*24.))


Function documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: transitDuration
.. autofunction:: transitDuration_Rs
.. autofunction:: ingressDuration
.. autofunction:: ingressDuration_Rs


Estimate times of contact (T1-T4) for eccentric orbit
---------------------------------------------------------

Numerically estimate contact points for primary and secondary eclipse.

Example
~~~~~~~~

::
    
    from __future__ import print_function
    from PyAstronomy import pyasl
    
    # SMA in stellar radii
    sma = 5.67
    # Rp/Rs
    rprs = 0.15
    # Orbital inclination
    inc = 89.2
    # Orbital period (time units are arbitrary but must be consistent)
    p = 2.0
    # Eccentricity
    e = 0.63
    # Argument of periastron (planetary orbit)
    w = 155.
    # Time of periastron passage
    tau = 2412345.346 
    
    # Contact times for primary transit
    pts = pyasl.transit_T1_T4_ell(sma, rprs, inc, p, tau, e, w, transit="p")
    # Contact times for secondary transit
    sts = pyasl.transit_T1_T4_ell(sma, rprs, inc, p, tau, e, w, transit="s")
    
    print("Transit times at arbitrary epoch (N*period may added)")
    print("Primary transit T1-4: ", pts)
    print("Secondary trabnsit T1-4: ", sts)
    print()
    print("Duration of primary and secondary transit: %5.3f, %5.3f " % (pts[3]-pts[0], sts[3]-sts[0]))


Function documentation
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: transit_T1_T4_ell

