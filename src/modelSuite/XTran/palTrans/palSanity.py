import unittest
import numpy as np


class PalTransExampleSanity(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def _rr(self, low, high):
    return (np.random.random()*(high - low) + low)
  
  def sanity_CircVsKep(self):
    """
      Check that Keplerian and circular solutions do not differ in case of circular orbit.
    """
    from PyAstronomy.modelSuite.XTran import palTrans as pt
    c = pt.PalLC()
    k = pt.PalLCKep()
    # Ensure that `k` has a properly defined circular orbit
    k["w"] = 0.0; k["Omega"] = 0.0; k["e"] = 0.0
    # A time axis
    x = np.linspace(-5., 5., 3000)
    for counter in xrange(5):
      # Get some random parameters
      # p, a, i, linLimb, quadLimb, T0, per, b
      linLimb = self._rr(0.,1.)
      T0 = self._rr(-8., 3.)
      pars = {"p":self._rr(0.02, 0.2),
              "a":self._rr(2., 10),
              "i":self._rr(86., 94.),
              "linLimb":linLimb,
              "quadLimb":self._rr(0,1.-linLimb),
              "per":self._rr(1.,4.),
              "b":self._rr(0., 0.1)
              }
      c.assignValue(pars)
      k.assignValue(pars)
      c["T0"] = T0
      k["T0pa"] = c["T0"] - k["per"]/4.
      yc = c.evaluate(x)
      yk = k.evaluate(x)
#      import matplotlib.pylab as plt
#      plt.cla()
#      
##      plt.plot(x, c._zlist, 'b--')
##      plt.plot(x, k._zlist, 'r--')
##      plt.show()
#      
#      plt.plot(x, yc)
#      plt.plot(x, yk)
#      plt.plot(x, yc/yk)
#      plt.show()
      np.testing.assert_almost_equal(yc, yk, \
                                     err_msg="Circular and Keplerian Pal models differ.\n"+"\n".join(k.parameterSummary(toScreen=False))+"\n\n"+ \
                                     "\n".join(c.parameterSummary(toScreen=False)) , decimal=10, verbose=True)
