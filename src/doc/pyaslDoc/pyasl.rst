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


Available functionality
---------------------------  

Below, you find a list of available functionality roughly sorted
by category.

Uncategorized algorithms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentModule:: PyAstronomy.pyasl

.. toctree::
   :maxdepth: 1
   
   aslDoc/keplerOrbit.rst
   aslDoc/magnitudes.rst
   aslDoc/sindex.rst
   aslDoc/atomicNo.rst
   aslDoc/aslExt_1Doc/abundances.rst
   aslDoc/aslExt_1Doc/pizzolato2003Doc.rst
   aslDoc/aslExt_1Doc/ramirez2005.rst
   aslDoc/aslExt_1Doc/stellarAge.rst

General data analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
   
   aslDoc/folding.rst
   aslDoc/binning.rst
   aslDoc/smooth.rst
   aslDoc/estimateSNR.rst
   aslDoc/outlier.rst
   aslDoc/intep.rst
   aslDoc/quadextreme.rst
   aslDoc/svd.rst

Spectral routines
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
   
   aslDoc/pyasl_wvlconv.rst
   aslDoc/baryvel.rst
   aslDoc/planck.rst
   aslDoc/crosscorr.rst
   aslDoc/dopplerShift.rst
   aslDoc/broadening.rst
   aslDoc/unredDoc.rst
   aslDoc/lineWidth.rst
   aslDoc/fluxConversion.rst
   aslDoc/readFitsSpec.rst
   aslDoc/aslExt_1Doc/transCurves.rst


Coordinates, time, and observation tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
   
   aslDoc/coorTimeDate.rst
   aslDoc/decimalYear.rst
   aslDoc/pyaslAitoff.rst
   aslDoc/sunpos.rst
   aslDoc/nutate.rst
   aslDoc/aberration.rst
   aslDoc/altitude.rst
   aslDoc/eq2hor.rst
   aslDoc/moon.rst
   aslDoc/cardinalPoint.rst
   aslDoc/twilight.rst
   aslDoc/angularDistance.rst
   aslDoc/posAngle.rst
   aslDoc/airmass.rst
   aslDoc/observatory.rst
   aslDoc/coordinates.rst

Transit and planet routines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   aslDoc/transitDuration.rst
   aslDoc/inTran.rst
   aslDoc/transitVisibility.rst
   aslDoc/planetPhase.rst

Statistics
~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   aslDoc/statTest.rst
   aslDoc/aslExt_1Doc/expCorrRN.rst

Resource based helpers
~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
   
   resBasedDoc/baraffe98tracks.rst
   resBasedDoc/nasaExoplanetArchive.rst
   resBasedDoc/exoplanetEU.rst
   resBasedDoc/exoplanetsOrg.rst
   resBasedDoc/sweet_cat.rst
   resBasedDoc/kuruczModels.rst
   resBasedDoc/fip.rst
   

   
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