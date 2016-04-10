Aitoff projection
==================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: aitoff

.. autofunction:: inverseAitoff

    
Example: Aitoff projection and its inverse
--------------------------------------------------

Carry out Aitoff projection and its inverse.

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Define longitude and latitude in degrees
    l = 130.
    b = -35.
    print("Input - Longitude: %4d deg, Latitude: %4d deg" % (l, b))
    
    print("Aitoff project them and ...")
    x, y = pyasl.aitoff(l, b)
    print(x, y)
    
    print("... get them back.")
    l2, b2 = pyasl.inverseAitoff(x, y)
    print(l2, b2)