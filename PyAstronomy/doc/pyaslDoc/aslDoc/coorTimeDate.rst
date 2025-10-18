Julian/calendar date and heliocentric correction
===================================================================

.. p23ready

Below you find the documentation of the following (mostly legacy) routines,
which provide several conversions involving Gregorian and
Julian dates as well as a heliocentric correction.

 * :py:func:`daycnv (Convert Julian into Gregorian dates) <PyAstronomy.pyasl.asl.astroTimeLegacy.daycnv>`
 * :py:func:`jdcnv (Convert Gregorian into Julian dates) <PyAstronomy.pyasl.asl.astroTimeLegacy.jdcnv>`
 * :py:func:`xyz (Geocentric 3d solar velocity and position) <PyAstronomy.pyasl.asl.astroTimeLegacy.xyz>`
 * :py:func:`helio_jd (heliocentric correction) <PyAstronomy.pyasl.asl.astroTimeLegacy.helio_jd>`
 * :py:func:`get_juldate (Current Julian date) <PyAstronomy.pyasl.asl.astroTimeLegacy.get_juldate>`
 * :py:func:`juldate (Compute RJD) <PyAstronomy.pyasl.asl.astroTimeLegacy.juldate>`
 * :py:func:`localTime (Approximation of local time) <PyAstronomy.pyasl.asl.localTime>`
 * :py:func:`weekday <PyAstronomy.pyasl.asl.astroTimeLegacy.weekday>`
 * :py:func:`bprecess (coordinate precession) <PyAstronomy.pyasl.asl.astroTimeLegacy.bprecess>`
 * :py:func:`premat (Precession matrix) <PyAstronomy.pyasl.asl.astroTimeLegacy.premat>`
 * :py:func:`precess (coordinate precession) <PyAstronomy.pyasl.asl.astroTimeLegacy.precess>`
 * :py:func:`precess_xyz (precess rectengular coordinates) <PyAstronomy.pyasl.asl.astroTimeLegacy.precess_xyz>`


.. _timescales:

Notes on timescales in astronomy
----------------------------------

Accurate representation of time requires to specify the time scale.
An overview about astronomical time scales is given by the IAU in the context of
their SOFA collection (`see <https://www.iausofa.org/current_C.html#Documents>`_); see particularly
their cookbook on *Time Scale and Calendar Tools*. Also Eastman et al. 2010 (PASP 122, 935) present
a comprehensive summary. 
Modern IAU definitions of time systems are based on
General Relativity (GR) considerations. 
Below is a very informal presentation of some frequently encountered terms.

The main aspects of GR relevant here are:
(a) Sitting in a gravitational well dilates time, (b) relative motion dilates time. This means that the notion if time
depends on one's location in the potential well of the Earth, Sun, etc.. 

	* TAI : Atomic time (enumerates constant seconds ticking on surface of Earth). This is (almost) proper time on the surface of the Earth.
	* TT : ``Terrestial time'' = TAI + 32.184 s (convention).
	* UTC : Essentially TAI **but** with leap seconds introduced to adjust to variable rotation of Earth (solar time). JD is **not** a unique representation of UTC because of leap seconds (days may be longer or shorter). 
	* TCG : ``Geocentric coordinate time'' is time, which would be recorded in the *Geocentric celestial reference system*, a metric (spacetime) comoving with the Earth but nonrotating and without the potential well of the Earth. Time ticks a bit faster here than TT and the relation is linear. 
	* TCB : ``Barycentric coordinate time'' is time, which would be recorded in the *barycentric celestial reference system*, a metric (spacetime) comoving with the barycenter of the solar system but without its potential wells (no Sun, no planets). Time ticks again a bit faster here and relation to TCG is complex.
	* TDB : ``Barycentric dynamical time'' is essentially TCB but slowed down by a linear transformation to make it (on average) go in step with TT.

To change the timescale, `astropy.time <https://docs.astropy.org/en/stable/time/index.html>`_ provides all the necessary tools
in Python, ultimately based on the IAU's SOFA routines.

*Which timescale to use?* The answer naturally depends on your needs. For event outside of the solar system TDB is a common choice.
At any rate, try to get rid of UTC's leap second discontinuity nightmare asap unless it is specifically required.

.. note:: Despite the name ``barycentric'' in TCB and TDB, the scale has nothing to do with light travel time corrections. In fact,
		  BJD, HJD, or geocentric Julian Dates can be expressed in any timescale.
		  


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
    print("year = %4d, month = %2d, day = %2d, hour = %5.3f"
          % tuple(pyasl.daycnv(jd)))
    print()
    
    # Convert calendar date to JD
    dt = datetime.datetime(1968, 5, 23, 12)
    print("Input date: ", dt)
    print("Corresponding Julian date: ", pyasl.jdcnv(dt))
    print("Corresponding reduced Julian date: ", pyasl.juldate(dt))
    print("Corresponds to weekday: ", pyasl.weekday(dt))
    print()
    
    # Get current Julian date and convert to heliocentric date
    cjd = pyasl.get_juldate()
    print("Current Julian date: ", cjd)
    print("Current (reduced) heliocentric Julian date (ra=100 deg, dec=37 deg): ",
          pyasl.helio_jd(cjd-2.4e6, 100.0, 37.0))


Routines
---------------

.. automodule:: PyAstronomy.pyasl.asl.astroTimeLegacy
   :members:

.. autofunction:: PyAstronomy.pyasl.asl.localTime