Nutation
===========

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

.. autofunction:: nutate

Example: `nutate`
--------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import datetime
    import numpy as np
    
    # Convert calendar date into JD
    # use the datetime package
    jd = datetime.datetime(2013, 4, 16)
    jd = pyasl.jdcnv(jd)
    print("Nutation for the date.")
    res = pyasl.nutate(jd)
    print("JD = " + str(jd) + ", Longitude = " + str(res[0]) + \
          ", Obliquity = " + str(res[1]))
    
    # Get nutation for an array of JDs.
    startjd = datetime.datetime(2013, 4, 16)
    endjd = datetime.datetime(2013, 6, 16)
    startjd = pyasl.jdcnv(startjd)
    endjd = pyasl.jdcnv(endjd)
    jds = np.arange(startjd, endjd, .5)
    print()
    print("Plot the results")
    res = pyasl.nutate(jds, plot=True)
    
    print("Longitude: ", res[0])
    print("Obliquity: ", res[1])


---------------------------------------------------------------------------

.. autofunction:: co_nutate

Example: `co_nutate`
-------------------------
  
::

    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import datetime
    import numpy as np
    
    # Convert calendar date into JD
    # use the datetime package
    jd = datetime.datetime(2013, 4, 16)
    jd = pyasl.jdcnv(jd)
    # Specify RA and DEC (degrees)
    ra = 10.
    dec = 30.
    print("Get change in RA and DEC due to Earth's nutation for JD = " \
          + str(jd))
    print(pyasl.co_nutate(jd, ra, dec))
    
    print()
    print("Get change for several RAs and DECs for the same JD")
    ra = np.arange(0.,160.,20.)
    dec = np.arange(-80.,80.,20.)
    res = pyasl.co_nutate(np.repeat(jd, ra.size), ra, dec)
    print(res[0], res[1])
    
    print()
    print("Get change for several RAs and DECs for different JDs")
    jds = np.arange(jd,jd+ra.size,1)
    res = pyasl.co_nutate(jds, ra, dec)
    print("JD             delta(RA)   delta(DEC)")
    for i in range(ra.size):
      print("%12.5f   %8.5f   %8.5f" % (jds[i], res[0][i], res[1][i]))
  
