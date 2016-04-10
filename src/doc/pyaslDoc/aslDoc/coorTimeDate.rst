Julian/calendar date and heliocentric correction
===================================================================

.. p23ready

Below you find the documentation of the following routines:

 * :py:func:`daycnv <PyAstronomy.pyasl.asl.astroTimeLegacy.daycnv>`
 * :py:func:`bprecess <PyAstronomy.pyasl.asl.astroTimeLegacy.bprecess>`
 * :py:func:`premat <PyAstronomy.pyasl.asl.astroTimeLegacy.premat>`
 * :py:func:`precess <PyAstronomy.pyasl.asl.astroTimeLegacy.precess>`
 * :py:func:`precess_xyz <PyAstronomy.pyasl.asl.astroTimeLegacy.precess_xyz>`
 * :py:func:`xyz <PyAstronomy.pyasl.asl.astroTimeLegacy.xyz>`
 * :py:func:`helio_jd <PyAstronomy.pyasl.asl.astroTimeLegacy.helio_jd>`
 * :py:func:`jdcnv <PyAstronomy.pyasl.asl.astroTimeLegacy.jdcnv>`
 * :py:func:`get_juldate <PyAstronomy.pyasl.asl.astroTimeLegacy.get_juldate>`
 * :py:func:`juldate <PyAstronomy.pyasl.asl.astroTimeLegacy.juldate>`
 * :py:func:`localTime <PyAstronomy.pyasl.asl.localTime>`

Example: Convert from JD into calendar date (daycnv)
--------------------------------------------------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Convert JD to calendar date
    jd = 2440000.0 + 18614./(24.*3600.)
    print("year = %4d, month = %2d, day = %2d, hour = %5.3f" \
          % tuple(pyasl.daycnv(jd)))
    print()
    print("year = %4d, month = %2d, day = %2d, hour = %2d, minute = %2d, seconds = %2d, microseconds = %6d" \
          % tuple(pyasl.daycnv(jd, mode='dtlist')))
    print()
    dt = pyasl.daycnv(jd, mode='dt')
    print("Datetime object: ", dt)



Example: daycnv, helio_jd, jdcnv, get_juldate, and juldate
------------------------------------------------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import datetime
    
    # Convert JD to calendar date
    jd = 2440000.0
    print("year = %4d, month = %2d, day = %2d, hour = %5.3f" \
          % tuple(pyasl.daycnv(jd)))
    print()
    
    # Convert calendar date to JD
    dt = datetime.datetime(1968, 5, 23, 12)
    print("Input date: ", dt)
    print("Corresponding Julian date: ", pyasl.jdcnv(dt))
    print("Corresponding reduced Julian date: ", pyasl.juldate(dt))
    print()
    
    # Get current Julian date and convert to heliocentric date
    cjd = pyasl.get_juldate()
    print("Current Julian date: ", cjd)
    print("Current (reduced) heliocentric Julian date (ra=100 deg, dec=37 deg): ", \
          pyasl.helio_jd(cjd-2.4e6, 100.0, 37.0))


Routines
---------------

.. automodule:: PyAstronomy.pyasl.asl.astroTimeLegacy
   :members:

.. autofunction:: PyAstronomy.pyasl.asl.localTime