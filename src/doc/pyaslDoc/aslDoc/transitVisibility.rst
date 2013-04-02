Transit visibility
====================

.. currentmodule:: PyAstronomy.pyasl

The function :py:func:`transitTimes` provides a table containing, among others,
mid-transit times, altitude above horizon, and Moon distance and phase; i.e., the most relevant
information to evaluate whether a transit is observable at a particular point.

Using :py:func:`transitVisibilityPlot`, the data can be plotted to obtain a
visual impression of the visibility.

Example of usage
------------------

The following example demonstrates the use of the functions. Both
code and oputput is shown.

Code
~~~~~~

::

  from PyAstronomy import pyasl
  import datetime as dt
  
  # Get the data for WASP-7 from NEXA data base
  nexa = pyasl.NasaExoplanetArchive()
  dat = nexa.selectByPlanetName("Wasp-7 b")
  
  # Convert July 14th, 2018, 10pm to a Julian date
  d = dt.datetime(2018, 7, 14, 22)
  jd = pyasl.jdcnv(d)
  
  # Calculate transit data for transits between July 14th and
  # July 24th, 2018.
  dat = pyasl.transitTimes(jd, jd+10.0, dat, nexaInput=True, \
                           observatory="eso", obsOffset=1./24.)
  # Plot the result
  pyasl.transitVisibilityPlot(dat, markTransit=True)


Text output
~~~~~~~~~~~~

::

  Specified time span
  Start date (DDDD-MM-YY and fractional hours): 2018-07-14 21.962
  End date (DDDD-MM-YY and fractional hours): 2018-07-24 21.983
  
  Transit duration:  226.512  minutes
  Off-transit time before and after transit:  60.0  minutes
  
  No. Tmid [HJD]      Obs. start [UT] [ALT, DIR(AZI)]     Transit mid [UT] [ALT, DIR(AZI)]     Obs. end [UT] [ALT, DIR(AZI)]   twilight (SUN ALT)                   moon distance     moon phase
    1 2458310.43133   10. 7. 19:27    [ 70°,E(126°)]      10. 7. 22:21     [ 68°,W(236°)]      11. 7.  1:14  [ 35°,W(240°)]     astron. twilight (-82°,-49°,-12°)   (127°,129°,131°)    7%
    2 2458315.38599   15. 7. 18:22    [ 62°,E(119°)]      15. 7. 21:15     [ 76°,S(221°)]      16. 7.  0:09  [ 43°,W(242°)]                night (-73°,-63°,-25°)   (147°,146°,144°)   12%
    3 2458320.34065   20. 7. 17:17    [ 53°,E(117°)]      20. 7. 20:10     [ 80°,S(182°)]      20. 7. 23:03  [ 52°,W(242°)]                night (-59°,-76°,-39°)   ( 86°, 85°, 83°)   61%


Graphical output
~~~~~~~~~~~~~~~~~

The figure shows the product of the `transitVisibilityPlot` function.
In this case, three transits are shown in the form of altitude (or
airmass) as a function of UT and local time. The numbers indicate
to which transit the individual lines refer. The color of the thin
lines indicates the light conditions (e.g., twilight or astronomical night).
The green overlay indicates the time of transit; here, we have one hour of
additional time before and after the transit, which was specified on
the call to `transitTimes`. On the right, you see the dates at which the
*observation* (not necessarily the transit) starts.

.. image:: transVis-WASP-7b.png
   :height: 600px


API
--------

.. autofunction:: transitTimes

.. autofunction:: transitVisibilityPlot
