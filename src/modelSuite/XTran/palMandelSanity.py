from __future__ import division
import unittest
import numpy as np

class palMandelSanity(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def _rr(self, low, high):
    return (np.random.random()*(high - low) + low)
  
  def sanity_comapare(self):
    """
      Compare Mandel's occultquad and the Pal results. 
    """
    from PyAstronomy.modelSuite import palTrans
    from PyAstronomy.modelSuite import forTrans
    pal = palTrans.PalLC()
    onl = forTrans.MandelAgolLC()
    x = np.linspace(0,3,500)
    for counter in range(5):
      linLimb = self._rr(0, 1)
      par = {"p":self._rr(0.,0.2),
             "a":self._rr(2.,5.),
             "i":self._rr(86, 94),
             "linLimb":linLimb,
             "quadLimb":self._rr(0, 1.-linLimb),
             "T0":self._rr(0,1),
             "per":self._rr(2,3),
             "b":self._rr(0,0.1)
             }
      pal.assignValue(par)
      onl.assignValue(par)
      yp = pal.evaluate(x)
      ym = onl.evaluate(x)
      
      np.testing.assert_almost_equal(yp, ym, 5, "The Pal and Mandel models differ."+"\n\n"+"\n".join(pal.parameterSummary(False)) + \
                                     "\n\n"+"\n".join(onl.parameterSummary(False)) , True)
      
  def sanity_basicCirc(self):
    """
      Checking basic MA and Pal transit configuration. Must be identical.
    """
    from PyAstronomy.modelSuite import XTran as xt
    mac = xt.forTrans.MandelAgolLC()
    mak = xt.forTrans.MandelAgolLC(orbit="keplerian")
    
    mac["p"] = 0.16
    mac["per"] = 2.1
    mac["linLimb"] = 0.37
    mac["a"] = 5.4
    
    for p in ["p", "per", "linLimb", "a"]:
      mak[p] = mac[p]
    
    t = np.linspace(0,4,1001)
    cc = mac.evaluate(t)
    ck = mak.evaluate(t)
    self.assertAlmostEqual(np.max(np.abs(cc-ck)), 0.0, delta=1e-12, \
                           msg="MandelAgol: Default for circ and kep orbit not identical")
    
    pa = xt.palTrans.PalLC()
    for p in ["p", "per", "linLimb", "a"]:
      pa[p] = mac[p]
    cp = pa.evaluate(t)
    self.assertAlmostEqual(np.max(np.abs(cc-cp)), 0.0, delta=1e-12, \
                           msg="MandelAgol vs. Pal: Default for circ orbit not identical")