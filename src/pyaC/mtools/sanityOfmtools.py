import unittest
import os

class SanityOfmtools(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    if os.path.isfile("test.tmp"):
      os.remove("test.tmp")
  
  def sanity_numericalDerivative(self):
    """
      mtools: Checking accuracy of numerical derivatives.
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
      mtools: diffCFD example
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
