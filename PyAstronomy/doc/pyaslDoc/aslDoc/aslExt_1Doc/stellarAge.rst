Stellar ages
=================

.. p23ready

Below, algorithms for estimating stellar ages based on
rotation and chromopsheric activity are given.

Gyrochronological age
-----------------------

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: gyroAgeBarnes

Example
~~~~~~~~

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Parameters of the Sun (Barnes 2007, p 1174)
    bv = 0.642
    p = 26.09
    
    # Obtain solar age ...
    age = pyasl.gyroAgeBarnes(p, bv)
    # ... and print it
    print("Solar age: {0:4.2f} +/- {1:4.2f} Ga".format(*age))


Chromospheric age
--------------------

.. autofunction:: chromoAgeRHK

Example
~~~~~~~~~

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Approximate chromospheric age of the Sun
    print("Solar age: {0:4.2f} Ga".format(pyasl.chromoAgeRHK(-4.95)))