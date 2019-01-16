Estimate the transit duration
==============================

.. p23ready
.. currentmodule:: PyAstronomy.pyasl

Functions to estimate the duration of the transit and the in- and egress based on input in solar system units (AU, solar and Jovian radius)

* :py:func:`Transit duration (solar system units) <transitDuration>`
* :py:func:`Ingress duration (solar system units) <ingressDuration>`

or parameters expressed in terms of the stellar radius

* :py:func:`Transit duration (stellar units) <transitDuration_Rs>`
* :py:func:`Ingress duration (stellar units) <ingressDuration_Rs>`


Example: Estimate duration of Earth's transit and ingress (solar system units)
---------------------------------------------------------------------------------

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
-------------------------------------------------------------------------------

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
-----------------------------

.. autofunction:: transitDuration
.. autofunction:: transitDuration_Rs
.. autofunction:: ingressDuration
.. autofunction:: ingressDuration_Rs

