Galactic space velocity (UVW)
=================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

.. autofunction:: gal_uvw


Example
--------

::

    from __future__ import print_function
    from PyAstronomy import pyasl
    
    ra, dec = 017.4294147291650, +61.5473037149441
    pmra, pmdec = 627.89, 77.84
    d = 144
    rv = -321.4
    
    print("U, V, W [kms/s]: % 5.3f, % 5.3f, % 5.3f " % \
        pyasl.gal_uvw(ra, dec, pmra, pmdec, d, rv))
    print("U, V, W (LSR corrected) [kms/s]:  % 5.3f, % 5.3f, % 5.3f" % \
        pyasl.gal_uvw(ra, dec, pmra, pmdec, d, rv, lsr=(8.5, 13.38, 6.49) ) )
