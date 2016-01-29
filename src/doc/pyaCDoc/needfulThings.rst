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
  
  from __future__ import print_function
  from PyAstronomy import pyaC
  
  nl = pyaC.NestedLoop([4,2,5], lowerLimits=[0,0,-5])
  
  for indices in nl:
    print(indices)


Invert index selection
-----------------------

.. autofunction:: invertIndexSelection

Example
~~~~~~~~

::

    from __future__ import print_function
    import numpy as np
    from PyAstronomy import pyaC as pc
    
    # Create "data" values and select some
    x = np.exp(np.arange(20.0)/20.0)
    indi = np.where(np.logical_and(x > 1.4, x < 1.7))
    
    print("Selected indices and values:")
    print("  indices: ", indi)
    print("  values : ", x[indi])
    
    indiInv = pc.invertIndexSelection(x, indi)
    
    print()
    print("Inverted selection:")
    print("  indices: ", indiInv)
    print("  values : ", x[indiInv])
    
    # Check that the same result is obtained by simply
    # passing the length of the array `x`
    
    indiInvLen = pc.invertIndexSelection(len(x), indi)
    
    print()
    print("Are indiInv and indiInvLen are the same? ")
    print("  ", np.all(indiInvLen == indiInv))


Fuzzy word matching
--------------------

.. autofunction:: fuzzyMatch

Example
~~~~~~~~~

::

    from __future__ import print_function
    from PyAstronomy import pyaC
    
    wordList = ["one", "two", "three", "four", "o-ne"]
    
    r = pyaC.fuzzyMatch("One", wordList)
    print("Exact match: {em:}, close match(es): {cm:}".format(**r))
    
    r = pyaC.fuzzyMatch("One", wordList, cutoff=0.4)
    print("Exact match: {em:}, close match(es): {cm:}".format(**r))
    
    r = pyaC.fuzzyMatch("One", wordList, caseSensitive=False)
    print("Exact match: {em:}, close match(es): {cm:}".format(**r))

  
Convert degree into rad and vice versa
---------------------------------------

.. autofunction:: degtorad
.. autofunction:: radtodeg

    
Matrix (2d) output
-------------------

.. autofunction:: matrix2doutput

Example
~~~~~~~~

::

    from __future__ import print_function
    from PyAstronomy import pyaC as PC
    import numpy as np
    
    m = np.zeros((4,3))
    colNames = ["first", "second", "third"]
    rowNames = ["1", "2", "3", "4"]
    
    for j in range(4):
      for i in range(3):
        m[j,i] = (i+1) * 10**(j+1) * np.sin(i*j)
    
    PC.matrix2doutput(m, colNames=colNames, rowNames=rowNames)
    print()
    PC.matrix2doutput(m, rowNames=rowNames)
    print()
    PC.matrix2doutput(m, colsep=" ")
    print()
    PC.matrix2doutput(m, oformat="% 12.5f", colNames=colNames, rowNames=rowNames)
    print()
    PC.matrix2doutput(m, oformat=["% 12.5f", "% 6.3f", "% e"], colNames=colNames)



Simple Input-Output file
-------------------------

.. autoclass:: SimIOF
   :members:

::
  
    from __future__ import print_function
    from PyAstronomy import pyaC as PC
    import numpy as np
    
    f = PC.SimIOF("origin", "test.tmp", 'w')
    
    a = 5.6
    b = 8.7
    c = 5
    
    f.addProp("a", a, fmt="% 4.4e")
    f.addProp(["b", "c"], [b,c])
    
    for x in range(10):
      f.write(str(x) + "\n")
    
    f.close()
    
    g = PC.SimIOF("origin", "test.tmp", 'r')
    
    # See the properties assigned
    print(g.args)
    # Use loadtxt to load the data
    print(np.loadtxt(g))