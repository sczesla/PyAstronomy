Periodograms
===============

Using the classes presented below, several periodograms can be
calculated:

  * :class:`Fourier transform <PyAstronomy.pyTiming.pyPeriod.Fourier>`
  * :class:`Lomb Scargle <PyAstronomy.pyTiming.pyPeriod.LombScargle>`
  * :class:`Generalized Lomb Scargle <PyAstronomy.pyTiming.pyPeriod.Gls>`

All the presented classes inherit from the PeriodBase
class.

.. currentmodule:: PyAstronomy.pyTiming.pyPeriod

The Fourier Transform
------------------------------------------

.. autoclass:: Fourier
   :members:
   :inherited-members:


The Lomb-Scargle-Periodogram (*fast*)
------------------------------------------

.. autoclass:: LombScargle
   :members:
   :inherited-members:


The generalized Lomb-Periodogram
------------------------------------------

.. autoclass:: Gls
   :members:
   :inherited-members:


References
-------------

.. [vdK] van der Klis, *Fourier Techniques In X-Ray Timing*
.. [Leahy83] Leahy et al. 1983, "On searches for pulsed emission with application to four globular cluster X-ray sources - NGC 1851, 6441, 6624, and 6712", 1983ApJ...266..160L 
.. [ZK09] Zechmeister & Kuerster 2009, "The generalised Lomb-Scargle periodogram. A new formalism for the floating-mean and Keplerian periodograms", 2009A&A...496..577Z
.. [HB86] Horne & Baliunas 1986, "A prescription for period analysis of unevenly sampled time series", 1986ApJ...302..757H
.. [Scargle82] Scargle 1982, "Studies in astronomical time series analysis. II - Statistical aspects of spectral analysis of unevenly spaced data", 1982ApJ...263..835S
.. [NR] *Numerical Recipes in C*, Second edition