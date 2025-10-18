Apply proper motion
=======================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: applyProperMotion

Example: Proper motion of GJ 1215
------------------------------------------------

::

    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Year 2018 coordinates and proper motion (mas/yr)
    # of GJ 1215
    ra, dec = 259.4319131981014, +11.6678903453170
    pmra, pmdec = -352.897, -445.558
    
    # Calculate 2050 position
    ra5, dec5 = pyasl.applyProperMotion(ra, dec, pmra, pmdec, dt=2050-2018, fixes=1)
    
    print("Position 2018: %10.6f  % 10.6f" % (ra, dec))
    print("Position 2050: %10.6f  % 10.6f" % (ra5, dec5))
