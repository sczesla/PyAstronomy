Harmonic timing analysis using periodograms
=================================================

The search for harmonic signals in the presence of noise is a
fundamental problem of timing analysis. This package provides
a collection of periodogram implementations to approach this
problem.

The cornerstone of this package is the implementation of the "Generalized Lomb-Scargle
:ref:`(GLS) <GLS>` periodogram". The class is also available as a
stand-alone module via the page of `Mathias Zechmeister <http://www.astro.physik.uni-goettingen.de/~zechmeister/>`_.

.. toctree::
   :maxdepth: 2
   
   gls.rst


Classical Lomb-Scargle and FFT
---------------------------------------

PyAstronomy provides a class to calculate the Fast Fourier transform (Fourier).
and a class to calculate the classical Lomb-Scargle periodogram (LombScargle).
The latter can be treated as a special case of the GLS and the :ref:`GLS <GLS>`
class can be used to obtain it.

.. toctree::
   :maxdepth: 2
   
   baseClasses.rst
   periodograms.rst
   examples.rst

   
References
-------------

.. [vdK] van der Klis, *Fourier Techniques In X-Ray Timing*
.. [Leahy83] Leahy et al. 1983, "On searches for pulsed emission with application to four globular cluster X-ray sources - NGC 1851, 6441, 6624, and 6712", 1983ApJ...266..160L 
.. [ZK09] Zechmeister & Kuerster 2009, "The generalised Lomb-Scargle periodogram. A new formalism for the floating-mean and Keplerian periodograms", 2009A&A...496..577Z
.. [HB86] Horne & Baliunas 1986, "A prescription for period analysis of unevenly sampled time series", 1986ApJ...302..757H
.. [Scargle82] Scargle 1982, "Studies in astronomical time series analysis. II - Statistical aspects of spectral analysis of unevenly spaced data", 1982ApJ...263..835S
.. [NR] *Numerical Recipes in C*, Second edition
