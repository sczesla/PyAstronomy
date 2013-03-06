Check whether time falls in transit window
===========================================

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: isInTransit


Example: Check individual point in time
----------------------------------------

::

  from PyAstronomy import pyasl
  
  # Time of interest
  time = 2476357.756234
  # Define some (arbitrary) transit parameters
  T0 = 2475123.01245
  period = 3.4789112
  duration = 2.2/24.0
  
  # Check whether the time is in-transit
  print "Time is within transit? ",
  if not pyasl.isInTransit(time, T0, period, duration/2.0):
    print "No"
  else:
    print "Yes"


Example: Checking a series of times
------------------------------------
    
::

  from PyAstronomy import pyasl
  import numpy as np
  
  # Times of interest
  times = 2476357.756234 + np.linspace(0.0, 5.0, 300)
  # Define some (arbitrary) transit parameters
  T0 = 2475123.01245
  period = 3.4789112
  duration = 2.2/24.0
  
  # Check whether the time is in-transit
  print "Indices if time points within transit: ",
  print pyasl.isInTransit(times, T0, period, duration/2.0)
  
  print
  print "For each time point, a flag indicating whether it"
  print "is in- or off-transit:"
  print pyasl.isInTransit(times, T0, period, duration/2.0, boolOutput=True)