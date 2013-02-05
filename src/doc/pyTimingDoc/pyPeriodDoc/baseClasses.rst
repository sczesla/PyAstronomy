Base classes
=================

Some commonly needed functionality and data structures are implemented
by the below listed classes. These are needed throughout the package.

The *TimeSeries* base class
------------------------------

This class provides an interface for a time series, viz., a time, flux, and (optionally)
an error array.

.. currentmodule:: PyAstronomy.pyTiming.pyPeriod
.. autoclass:: TimeSeries
   :members:

The PeriodBase Base class
------------------------------------------

This class provides functionality and an interface needed in every
periodogram class. The latter are to inherit from this class so that
the functionality and interface can be shared.

.. autoclass:: PeriodBase
   :members:
