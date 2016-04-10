Get angular distance from RA and DEC
======================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: getAngDist

Example
---------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    print("Angular distance between the poles (deg):")
    print(pyasl.getAngDist(98.0, -90.0, 100., +90.0))
    
    print("Angular distance between Vega and Altair (deg)")
    print(pyasl.getAngDist(279.23473479, +38.78368896,297.69582730, +08.86832120))  
