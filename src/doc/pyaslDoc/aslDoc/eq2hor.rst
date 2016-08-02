Horizontal coordinates
=========================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

Convert hour angle and declination into horizontal coordinates
----------------------------------------------------------------

.. autofunction:: hadec2altaz

Example: 
~~~~~~~~~

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import datetime
    import numpy as np
    
    # Hour angle 0. means transiting the local meridian.
    ha = 0.
    # Declination of object
    dec = 30.
    # Latitude of the observer (here Hamburger Sternwarte)
    lat = +53.48
    print("Get altitude and azimuth of object in DEGREES")
    print(pyasl.hadec2altaz(ha, dec, lat))
    
    # List of coordinates
    ha = np.arange(0.,20.,5.)
    dec = np.arange(30.,50.,5.)
    lat = np.zeros(dec.size)+53.48
    print()
    print("Get altitude and azimuth for a list of objects from same observer latitude")
    altaz = pyasl.hadec2altaz(ha, dec, lat)
    print("alt: ", altaz[0])
    print("az: ", altaz[1])
    




Convert celestial coordinates (RA/DEC) into local horizon coordinates (ALT/AZ)
---------------------------------------------------------------------------------

.. autofunction:: eq2hor

Example: 
~~~~~~~~~

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    import datetime
    import numpy as np
    
    # Convert calendar date to JD
    # use the datetime package
    jd = datetime.datetime(2013, 4, 16)
    jd = pyasl.jdcnv(jd)
    # Specific RA and DEC
    ra = 10.
    dec = 30.
    print()
    print("Get horizontal coordinates (alt, az, ha) from JD, RA,")
    print("  and DEC for the Hamburger Sternwarte")
    print(pyasl.eq2hor(jd, ra, dec, observatory="HS"))
    
    print()
    print("From a list of Julian dates ...")
    jds = np.arange(jd,jd+1,.2)
    ras = np.zeros(jds.size) + ra
    decs = np.zeros(jds.size) + dec
    alt, az, ha = pyasl.eq2hor(jds, ras, decs, lon=-70.4042, lat=-24.6272, alt=2635.)
    
    for i in range(alt.size):
      print("JD = %g : alt = % g,  az = % g,  ha = % g" % (jds[i], alt[i], az[i], ha[i]))
    
    
    print()
    print("For one object and different times at the VLT...")
    jds = np.arange(jd-.25,jd+.25,.01)
    ras = np.zeros(jds.size) + 130.
    decs = np.zeros(jds.size) - 30.
    res = pyasl.eq2hor(jds, ras, decs, lon=-70.4042, lat=-24.6272, alt=2635.)
    
    plt.plot(jds, res[0])
    plt.xlabel("Julian date")
    plt.ylabel("Altitude [deg]")
    plt.show()
