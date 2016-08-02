Decimal representation of year
===============================

.. p23ready

Convert between decimal year representation and gregorian date.

.. note:: Due to the nature of the gregorian calendar, the decimal year
          does not represent a uniform measure of time, because some years
          are longer than others. 

.. currentmodule:: PyAstronomy.pyasl

.. autofunction:: decimalYear
.. autofunction:: decimalYearGregorianDate


Example:
----------

::
    
    from __future__ import print_function, division
    import datetime as dt
    from PyAstronomy import pyasl
      
    # Convert July 2nd, 1998, 12:30:59 into decimal
    # representation
    d = dt.datetime(1998,7,2,12,30,59)
    
    # Obtain float representation of decimal year
    decy = pyasl.decimalYear(d)
    print("Decimal representation: ", decy)
    
    # Convert back into gregorian date first
    print("The decimal year %10.5f correspond to " % decy + \
          pyasl.decimalYearGregorianDate(decy, "yyyy-mm-dd hh:mm:ss"))
    print(" ... or equivalently (y, m, d, h, m, s, ms): ", \
          pyasl.decimalYearGregorianDate(decy, "tuple"))
