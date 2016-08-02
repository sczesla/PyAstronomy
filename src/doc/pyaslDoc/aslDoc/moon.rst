Lunar phase and position
==========================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

.. autofunction:: moonpos

Example: Finding the position of the Moon
------------------------------------------

::
    
    from __future__ import print_function, division
    import datetime
    from PyAstronomy import pyasl
    import numpy as np
    
    # Convert calendar date to JD
    # using the datetime package
    jd = datetime.datetime(2013, 4, 16)
    jd = pyasl.jdcnv(jd)
    jd = np.arange(jd, jd + 20, 1)
    # Calculate Moon positions
    res = pyasl.moonpos(jd)
    
    print("%15s  %8s  %8s  %11s  %8s  %8s" % \
      ("JD", "RA", "DEC", "DIST", "GEOLON", "GEOLAT"))
    print("%15s  %8s  %8s  %11s  %8s  %8s" % \
      ("[d]", "[deg]", "[deg]", "[km]", "[deg]", "[deg]"))
    for i in range(jd.size):
      print("%15.4f  %8.4f  %8.4f  %11.4f  %8.4f  %8.4f" % \
        (jd[i], res[0][i], res[1][i], res[2][i], res[3][i], res[4][i]))


.. autofunction:: moonphase
      

Example: Find the Moon's illuminated fraction
--------------------------------------------------

::
    
    from __future__ import print_function, division
    import datetime
    from PyAstronomy import pyasl
    import numpy as np
    
    # Convert calendar date to JD
    # using the datetime package
    jd = datetime.datetime(2013, 4, 16)
    jd = pyasl.jdcnv(jd)
    jd = np.arange(jd, jd+20,1)
    mp = pyasl.moonphase(jd)
    
    print("%15s  %3s" % ("JD", "Phase"))
    for i in range(jd.size):
      print("%15.4f  %3d%%" % (jd[i], mp[i]*100.))
