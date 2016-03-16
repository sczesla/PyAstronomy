from __future__ import print_function
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
    from .numericalDerivatives import diffCFD
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
        print(i, erro)
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

  def sanity_ibtrapzExample(self):
    """
      mtools: Checking example of ibtrapz
    """
    from PyAstronomy.pyaC import mtools
    import numpy as np
    
    x = np.arange(-2.,2.01,0.1)
    y = x**3 + 1.7
    
    x0 = -1.375
    x1 = +1.943
    
    # Analytical value of integral
    analyt = 0.25*(x1**4 - x0**4) + 1.7*(x1-x0)
    
    print("Analytical value: ", analyt)
    print("ibtrapz: ", mtools.ibtrapz(x, y, x0, x1))

  def sanity_ibtrapz(self):
    """
      mtools: Checking ibtrapz
    """
    from PyAstronomy.pyaC import mtools
    import numpy as np
    
    x = np.arange(-2.,2.01,0.1)
    y = 2. * x
    
    x0 = -1.375
    x1 = +1.943
    
    # Analytical value of integral
    analyt = x1**2 - x0**2

    self.assertAlmostEqual(analyt, mtools.ibtrapz(x, y, x0, x1), delta=1e-10, msg="ibtrapz incorrect for linear function.")
    self.assertAlmostEqual((-1.9)**2-(-2.0)**2, mtools.ibtrapz(x, y, -2.0, -2.0+0.1), delta=1e-10, msg="ibtrapz incorrect for linear function (-2,-1.9).")
    self.assertAlmostEqual(0.0, mtools.ibtrapz(x, y, -2.0, +2.0), delta=1e-10, msg="ibtrapz incorrect for linear function (-2,+2).")
    
  def sanity_zerocross1dExample(self):
    """
      Checking sanity of zerocross1d example
    """
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import pyaC
    
    # Generate some 'data'
    x = np.arange(100.)**2
    y = np.sin(x)
    
    # Set the last data point to zero.
    # It will not be counted as a zero crossing!
    y[-1] = 0
    
    # Set point to zero. This will be counted as a
    # zero crossing
    y[10] = 0.0
    
    # Get coordinates and indices of zero crossings
    xc, xi = pyaC.zerocross1d(x, y, getIndices=True)
    
    # Plot the data
    plt.plot(x, y, 'b.-')
    # Add black points where the zero line is crossed
    plt.plot(xc, np.zeros(len(xc)), 'kp')
    # Add green points at data points preceding an actual
    # zero crossing.
    plt.plot(x[xi], y[xi], 'gp')
#     plt.show()

  def sanity_zerocross1d(self):
    """
      Checking sanity of zerocross1d
    """
    import numpy as np
    from PyAstronomy import pyaC
    
    x = np.arange(3.)
    y = np.array([0, 3, 0])
    
    xz = pyaC.zerocross1d(x, y)
    self.assertEqual(len(xz), 0, msg="Found zero crossing in 0,3,0 array (problem with first/last point).")
    
    y = np.array([-1., 1., -2.])
    xz, xi = pyaC.zerocross1d(x, y, getIndices=True)
    self.assertEqual(len(xz), 2, msg="Found the following zero crossings in -1,1,-2 array: " + str(xz))
    self.assertEqual(len(xz), len(xi), "Number of values and indicies is not identical.")
    
    self.assertAlmostEqual(np.max(np.abs(np.array([0.5, 1.+1./3.])-xz)), 0.0, delta=1e-8, msg="Found unexpected zero crossings: " + str(xz))
    self.assertEqual(np.max(np.abs(xi - np.array([0,1]))), 0, msg="Found unexpected indices: " + str(xi))
    