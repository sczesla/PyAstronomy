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