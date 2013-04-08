Decimal representation of year
===============================

.. currentmodule:: PyAstronomy.pyasl

.. autofunction:: decimalYear

Example:
----------

::

  import datetime as dt
  from PyAstronomy import pyasl
  
  # Convert July 2nd, 1998, 12:30:59 into decimal
  # representation
  d = dt.datetime(1998,7,2,12,30,59)
  
  print "Decimal representation: ", pyasl.decimalYear(d)