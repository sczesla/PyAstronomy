Estimate the transit duration
==============================

.. p23ready

Functions to estimate the transit duration.

.. currentmodule:: PyAstronomy.pyasl

* :py:func:`Transit duration (physical units) <transitDuration>`
* :py:func:`Transit duration (stellar units) <transitDuration_Rs>`


Example transitDuration: Estimate duration of Earth's transit
-----------------------------------------------------------------

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


Example transitDuration_Rs: Estimate transit duration of HD 189733 b
----------------------------------------------------------------------

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


Function documentation
-----------------------------

.. autofunction:: transitDuration
.. autofunction:: transitDuration_Rs