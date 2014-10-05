import unittest
from PyAstronomy import pyasl

class SanityOfPyaslExt1(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_pizzolatoExample(self):
    """
      Example of Pizzolato 2003 relations
    """
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    import numpy as np
    
    p = pyasl.Pizzolato2003()
    
    # Define array of rotation periods [days]
    prot = np.arange(0.2, 30, 0.1)
    
    lx = np.zeros(prot.size)
    lxlbol = np.zeros(prot.size)
    
    # B-V color of star
    bv = 0.7
    
    # Obtain ...
    for i in range(prot.size):
      # ... log10 of X-ray luminosity
      lx[i] = p.log10lxbv(bv, prot[i])[0]
      # ... and log10(Lx/Lbol)
      lxlbol[i] = p.log10lxlbolbv(bv, prot[i])[0]
    
    # Plot result
    plt.subplot(2,1,1)
    plt.plot(prot, lx, 'bp-')
    plt.subplot(2,1,2)
    plt.plot(prot, lxlbol, 'bp-')
#     plt.show()

  def sanity_pizzolato(self):
    """
      Pizzolato 2003 sanity
    """
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    import numpy as np
    
    p = pyasl.Pizzolato2003()
    
    x = p.log10lxbv(0.87, 1.0)
    self.assertAlmostEqual(x[0], 29.9, msg="Lx saturation level for b-v=0.87 does not match.", delta=1e-7)
    self.assertAlmostEqual(x[1], 0.3, msg="Lx error for saturation level for b-v=0.87 does not match.", delta=1e-7)
    
    x = p.log10lxlbolmass(1.15, 3.0)
    self.assertAlmostEqual(x[0], -4.1694, msg="Pizzolato relation (m=1.15, pr=3.0) failed (value="+str(x[0])+").", delta=1e-4)