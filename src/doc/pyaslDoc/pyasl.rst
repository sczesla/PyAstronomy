Python AstroLib
==========================

.. _AstroLib : http://idlastro.gsfc.nasa.gov/

The aim of this project is to bundle a number of helpful rather small-scale routines.
This comprises porting the well-known "IDL Astronomy User's Library" (AstroLib_)
to Python. As the corresponding IDL code has been in use for decades now, it is "tested by application".
Therefore, we try to stay as close as possible to the IDL implementation. Furthermore, helpful routines
will be added.

.. note:: This is not an attempt to port the entire IDL library at once. We do what we need! Help is
          appreciated.

License note
----------------
The "IDL Astronomy User's Library" is released in "public domain". That means that nobody holds the
copyright and, thus, nobody has the right to sublicense it (not even the authors).
As we remain as close as possible in code
and ship the documentation mainly as is, it would be quite cumbersome to actually
assign copyright to parts written by us. Therefore, code which was directly
ported from the IDL Astronomy User's Library is released in public domain and
is marked by a note in the respective documentation saying:

::

  This function was ported from the IDL Astronomy User's Library.

All other code not carrying this note (may include fresh implementations of functionality
already provided by the IDL Astronomy User's Library, newly introduced functions, etc.)
is released under the MIT license.

Available functionality
---------------------------  

Useful algorithms
~~~~~~~~~~~~~~~~~~

.. currentModule:: PyAstronomy.pyasl

.. toctree::
   :maxdepth: 1
   
   aslDoc/pyasl_wvlconv.rst
   aslDoc/pyaslAitoff.rst
   aslDoc/coorTimeDate.rst
   aslDoc/decimalYear.rst
   aslDoc/baryvel.rst
   aslDoc/folding.rst
   aslDoc/keplerOrbit.rst
   aslDoc/statTest.rst
   aslDoc/binning.rst
   aslDoc/planck.rst
   aslDoc/planetPhase.rst
   aslDoc/dopplerShift.rst
   aslDoc/rotBroad.rst
   aslDoc/broad.rst
   aslDoc/fluxConversion.rst
   aslDoc/svd.rst
   aslDoc/unredDoc.rst
   aslDoc/smooth.rst
   aslDoc/estimateSNR.rst
   aslDoc/outlier.rst
   aslDoc/intep.rst
   aslDoc/magnitudes.rst
   aslDoc/sunpos.rst
   aslDoc/nutate.rst
   aslDoc/aberration.rst
   aslDoc/altitude.rst
   aslDoc/eq2hor.rst
   aslDoc/observatory.rst
   aslDoc/moon.rst
   aslDoc/lineWidth.rst
   aslDoc/cardinalPoint.rst
   aslDoc/twilight.rst
   aslDoc/angularDistance.rst
   aslDoc/airmass.rst
   aslDoc/crosscorr.rst

Transit routines
~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   aslDoc/transitDuration.rst
   aslDoc/inTran.rst
   aslDoc/transitVisibility.rst


Resource based helpers
~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
   
   resBasedDoc/baraffe98tracks.rst
   resBasedDoc/nasaExoplanetArchive.rst
   resBasedDoc/exoplanetEU.rst
