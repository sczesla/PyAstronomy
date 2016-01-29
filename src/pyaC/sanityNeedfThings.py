from __future__ import print_function, division
import unittest
import os

class SanityOfNeedfulThings(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    if os.path.isfile("test.tmp"):
      os.remove("test.tmp")
  
  def sanity_invertIndexSelectionExample(self):
    """
      Example of invertIndexSelection.
    """
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
  
  def sanity_invertIndexSelection(self):
    """
      Sanity of invertIndexSelection.
    """
    import numpy as np
    from PyAstronomy import pyaC as pc
    x = np.arange(20.0)
    indi = np.where(np.logical_or(x < 5.0, x > 17.0))[0]
    i1 = pc.invertIndexSelection(x, indi)
    self.assertEqual(len(x), len(indi)+len(i1), "indi and inverse do not make up the original array.")
    for i in indi:
      self.assertFalse(i in i1, "found an index of indi in inverse.")
    i2 = pc.invertIndexSelection(len(x), indi)
    self.assertTrue(np.all(i1 == i2), "i1 and i2 are not the same.")
    i3 = pc.invertIndexSelection(x, i1)
    self.assertTrue(np.all(i3 == indi), "re-inverse of i1 does not match indi.")

  def sanity_fuzzyMatchExample(self):
    """
      Check example for fuzzyMatch
    """
    from PyAstronomy import pyaC
    
    wordList = ["one", "two", "three", "four", "o-ne"]
    
    r = pyaC.fuzzyMatch("One", wordList)
    print("Exact match: {em:}, close match(es): {cm:}".format(**r))
    
    r = pyaC.fuzzyMatch("One", wordList, cutoff=0.4)
    print("Exact match: {em:}, close match(es): {cm:}".format(**r))
    
    r = pyaC.fuzzyMatch("One", wordList, caseSensitive=False)
    print("Exact match: {em:}, close match(es): {cm:}".format(**r))
    
  def sanity_matrix2doutputExample(self):
    """
      Sanity of matrix2doutput.
    """
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

  def sanity_SimIOF_Example(self):
    """
      Checking example for simple input/output file.
    """
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