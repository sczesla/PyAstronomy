Conversion between sexagesimal and decimal coordinate representation
=====================================================================

.. p23ready

The routines below convert between sexagesimal and decimal coordinate
representations.  

Example: Convert between decimal and sexagesimal representation
------------------------------------------------------------------

::

    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Coordinates of HD 1 from SIMBAD
    hd1 = "00 05 08.83239 +67 50 24.0135"
    
    print("Coordinates of HD 1 (SIMBAD): ", hd1)
    
    # Obtain decimal representation
    ra, dec = pyasl.coordsSexaToDeg(hd1)
    print("Coordinates of HD 1 [deg]: %010.6f  %+09.6f" % (ra, dec))
    
    # Convert back into sexagesimal representation
    sexa = pyasl.coordsDegToSexa(ra, dec)
    print("Coordinates of HD 1 [sexa]: ", sexa)


Routines
----------

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: coordsSexaToDeg
.. autofunction:: coordsDegToSexa
.. autofunction:: hmsToDeg
.. autofunction:: degToHMS
.. autofunction:: dmsToDeg
.. autofunction:: degToDMS
