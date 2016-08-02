Position of the Sun
=======================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

.. autofunction:: sunpos

Example
-----------

::
    
    from __future__ import print_function, division
    import numpy as np
    from PyAstronomy import pyasl
    import datetime
    
    # Convert calendar date into JD
    # use the datetime package
    jd = datetime.datetime(2013, 4, 16)
    jd = pyasl.jdcnv(jd)
    print("JD = " + str(jd))
    pos = pyasl.sunpos(jd, full_output=True)
    print("Coordinates of the Sun (ra, dec): %g, %g" % (pos[1], pos[2]))
    print("Solar elongation = %g and obliquity = %g" % (pos[3], pos[4]))
    
    # Get the Sun's RA and DEC values for a period of time.
    startjd = datetime.datetime(2013, 4, 16)
    endjd = datetime.datetime(2013, 6, 16)
    # Convert into Julian dates
    startjd = pyasl.jdcnv(startjd)
    endjd = pyasl.jdcnv(endjd)
    print()
    pos = pyasl.sunpos(startjd, end_jd=endjd, jd_steps=10, plot=False, full_output=True)
    
    for i in range(len(pos[0])):
      print("At JD = %g: ra = %g, dec = %g" % (pos[0][i], pos[1][i], pos[2][i]))
    
