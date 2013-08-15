import unittest
import os

class SanityOfPyasl(unittest.TestCase):
  
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
    print "Exact match: {em:}, close match(es): {cm:}".format(**r)
    
    r = pyaC.fuzzyMatch("One", wordList, cutoff=0.4)
    print "Exact match: {em:}, close match(es): {cm:}".format(**r)
    
    r = pyaC.fuzzyMatch("One", wordList, caseSensitive=False)
    print "Exact match: {em:}, close match(es): {cm:}".format(**r)

  def sanity_numericalDerivative(self):
    """
      Checking accuracy of numerical derivatives.
    """
    # Check polynomial
    from numericalDerivatives import diffCFD
    import numpy as np
    
    # Analytical derivatives
    x = np.linspace(-10,10,1000)
    y = [np.poly1d((0.03, -0.31, 0.4, 0.35, 1.4))]
    for i in range(4):
      y.append(y[-1].deriv())
    
    for erro in [2,4,6]:
      for i in range(1,5):
        indi, der = diffCFD(x, np.polyval(y[0], x), i, erro)
        self.assertLess(np.max(np.abs(der/np.polyval(y[i], x[indi]) - 1.0)), 2e-2*pow(10.0,erro-2))
    
    # Check trigonometric
    y = [np.sin(x/10.0*2*np.pi+3)]
    y.append(2*np.pi/10.0 * np.cos(x/10.0*2*np.pi+3))
    y.append(-(2*np.pi/10.0)**2 * np.sin(x/10.0*2*np.pi+3))
    y.append(-(2*np.pi/10.0)**3 * np.cos(x/10.0*2*np.pi+3))
    y.append((2*np.pi/10.0)**4 * np.sin(x/10.0*2*np.pi+3))
    for erro in [2,4,6]:
      for i in range(1,5):
        indi, der = diffCFD(x, y[0], i, erro)
        self.assertLess(np.max(np.abs(der/y[i][indi] - 1.0)), 1e-3*pow(10.0,erro-2))
    
    # Check exponential
    y = np.exp(x)
    for erro in [2,4,6]:
      for i in range(1,5):
        print i, erro
        indi, der = diffCFD(x, y, i, erro)
        self.assertLess(np.max(np.abs(der/y[indi] - 1.0)), 1e-3*pow(10.0,erro-2))

  def sanity_diffCFDExample(self):
    """
      diffCFD example
    """
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
#      plt.show()

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
    print
    PC.matrix2doutput(m, rowNames=rowNames)
    print
    PC.matrix2doutput(m, colsep=" ")
    print
    PC.matrix2doutput(m, oformat="% 12.5f", colNames=colNames, rowNames=rowNames)
    print
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
    print g.args
    # Use loadtxt to load the data
    print np.loadtxt(g)