The model classes
=============================

The *PalLC* and *PalLCCirc* classes encapsulates all relevant calculations;
the latter considers only circular orbits, while the former, *PAlLC*,
uses the `KeplerEllipse` class from PyA's astrolib to consider a
Kepler orbit.
 
They inherits from funFit's OneDFit class.

.. currentmodule:: PyAstronomy.modelSuite.XTran.palTrans

Circular orbit
-----------------

.. autoclass:: PalLC
   :members:

Full Keplerian orbit
----------------------

.. autoclass:: PalLCKep
   :members:
