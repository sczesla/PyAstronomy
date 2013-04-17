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


Fuzzy word matching
--------------------

.. autofunction:: fuzzyMatch

Example
~~~~~~~~~

::

  from PyAstronomy import pyaC
  
  wordList = ["one", "two", "three", "four", "o-ne"]
  
  r = pyaC.fuzzyMatch("One", wordList)
  print "Exact match: {em:}, close match(es): {cm:}".format(**r)
  
  r = pyaC.fuzzyMatch("One", wordList, cutoff=0.4)
  print "Exact match: {em:}, close match(es): {cm:}".format(**r)
  
  r = pyaC.fuzzyMatch("One", wordList, caseSensitive=False)
  print "Exact match: {em:}, close match(es): {cm:}".format(**r)

  
Convert degree into rad and vice versa
---------------------------------------

.. autofunction:: degtorad
.. autofunction:: radtodeg


Computer numerical derivatives
--------------------------------

.. autofunction:: diffCFD

Example
~~~~~~~~

::

  from PyAstronomy import pyaC
  import matplotlib.pylab as plt
  import numpy as np
  
  x = np.linspace(-10,10,1000)
  
  # Computer polynomial and its derivatives
  # (quasi analytically)
  y = [np.poly1d((0.03, -0.31, 0.4, 0.35, 1.4))]
  for i in range(4):
    y.append(y[-1].deriv())
  
  # Compute derivates numerically and compare to
  # analytic solution
  erro = 2
  for i in range(1,5):
    indi, der = pyaC.diffCFD(x, np.polyval(y[0], x), i, erro)
    plt.plot(x[indi], np.polyval(y[i], x[indi]), 'b.')
    plt.plot(x[indi], der, 'r--')
    plt.show()
