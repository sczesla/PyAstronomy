Needful things
================

.. currentmodule:: PyAstronomy.pyaC

A collection of helpers, which have no place elsewhere.

Nested loop iterator
---------------------

.. autoclass:: NestedLoop

Example
~~~~~~~~

::
  
  from PyAstronomy import pyaC
  
  nl = pyaC.NestedLoop([4,2,5], lowerLimits=[0,0,-5])
  
  for indices in nl:
    print indices


Invert index selection
-----------------------

.. autofunction:: invertIndexSelection

Example
~~~~~~~~

::

  import numpy as np
  from PyAstronomy import pyaC as pc
  
  # Create "data" values and select some
  x = np.exp(np.arange(20.0)/20.0)
  indi = np.where(np.logical_and(x > 1.4, x < 1.7))
  
  print "Selected indices and values:"
  print "  indices: ", indi
  print "  values : ", x[indi]
  
  indiInv = pc.invertIndexSelection(x, indi)
  
  print
  print "Inverted selection:"
  print "  indices: ", indiInv
  print "  values : ", x[indiInv]
  
  # Check that the same result is obtained by simply
  # passing the length of the array `x`
  
  indiInvLen = pc.invertIndexSelection(len(x), indi)
  
  print
  print "Are indiInv and indiInvLen are the same? "
  print "  ", np.all(indiInvLen == indiInv)

  