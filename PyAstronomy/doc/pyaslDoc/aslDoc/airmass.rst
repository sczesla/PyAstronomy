Calculate the airmass
======================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

Plane-parallel atmosphere
--------------------------- 

.. autofunction:: airmassPP

Example: Airmass with plane-parallel atmosphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    print("Airmass for plane-parallel atmosphere")
    for za in range(0,70,10):
      print("Zenith angle: %2d deg, airmass = %7.2f" % \
        (za, pyasl.airmassPP(za)))


Airmass for spherical, homogeneous atmosphere
------------------------------------------------

.. autofunction:: airmassSpherical

Example: Airmass with homogeneous, spherical atmosphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    obsAlt = 2400.0
    
    for za in range(0,90,10):
      print("Zenith angle: %2d, airmass = %7.2f" % \
        (za, pyasl.airmassSpherical(za, obsAlt)))

