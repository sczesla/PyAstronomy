from __future__ import division
import unittest
import numpy as np
import six.moves as smo


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

    # A time axis
    x = np.linspace(-5., 5., 300)
    for counter in smo.range(5):
      # Get some random parameters
      # p, a, i, linLimb, quadLimb, T0, per, b
      linLimb = self._rr(0.,1.)
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
      c["T0"] = 0.0
      k["tau"] = 0.0 
      yc = c.evaluate(x)
      yk = k.evaluate(x)

      np.testing.assert_almost_equal(yc, yk, \
                                     err_msg="Circular and Keplerian Pal models differ.\n"+"\n".join(k.parameterSummary(toScreen=False))+"\n\n"+ \
                                     "\n".join(c.parameterSummary(toScreen=False)) , decimal=10, verbose=True)
