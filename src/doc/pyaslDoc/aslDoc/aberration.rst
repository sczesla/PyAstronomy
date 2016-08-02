Aberration
============

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: co_aberration

Example: 
--------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import datetime
    import numpy as np
    
    # Convert calendar date to JD
    # use the datetime package
    jd = datetime.datetime(2013, 4, 16)
    jd = pyasl.jdcnv(jd)
    # Specify RA and DEC
    ra = 10.
    dec = 30.
    print("Get change in RA and DEC due to annual aberration") 
    print("  for JD = " + str(jd) + ":", \
          np.ravel(pyasl.co_aberration(jd, ra, dec)))
    
    print()
    print("Get change for several RAs and DECs for the same JD")
    ra = np.arange(10.,50.,10.)
    dec = np.arange(30.,70.,10.)
    res = pyasl.co_aberration(np.repeat(jd, ra.size), ra, dec)
    print(res[0], res[1])
    
    print()
    print("Get change for several RAs and DECs for different JDs")
    jds = np.arange(jd,jd+ra.size,1)
    res = pyasl.co_aberration(jds, ra, dec)
    print("JD             delta(RA)   delta(DEC)")
    for i in range(ra.size):
      print("%12.5f   %8.5f   %8.5f" %(jds[i], res[0][i], res[1][i]))

