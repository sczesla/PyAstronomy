from __future__ import print_function, division
import unittest
import numpy
from PyAstronomy import funcFit as fuf
import six

class FuncFitSanity(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def parameterAssignment1_sanity(self):
    gf = fuf.GaussFit1d()
    origVars = ["A", "mu", "sig", "off", "lin"]
    vals     = [1.0, 2.0,  3.0,   4.0,    5.0]
    for k, v in zip(origVars, vals):
      gf[k] = v
    for k, v in zip(origVars, vals):
      self.assertEquals(v, gf[k])
    gf.assignValue(dict(zip(origVars, numpy.zeros(len(origVars)))))
    self.assertEquals(numpy.sum(list(gf.parameters().values())), 0.0)
  
  def combine1_sanity(self):
    gf = fuf.GaussFit1d()
    gff = gf + gf + gf
    for p in six.iterkeys(gff.parameters()):
      gff[p] = numpy.random.random()
      gff.thaw(p); gff.thaw([p,p])
      gff.freeze(p); gff.freeze([p,p])
      gff.setRestriction({p:[None, None]})
    for prop in ["A", "mu", "sig", "off", "lin"]:
      for c in [1,2,3]:
        gff[prop, "Gaussian", c] = numpy.random.random()
        s = (prop, "Gaussian", c)
        gff.thaw(s); gff.thaw([s,s])
        gff.freeze(s); gff.freeze([s,s])
        gff.setRestriction({s:[None, 10.]})

  def description_sanity(self):
    """
      Check sanity of 'description'
    """
    gf = fuf.GaussFit1d()
    gff = gf + gf
    self.assertEqual(gff.description(), "Gaussian(No. 1) + Gaussian(No. 2)", "Wrong description: " + gff.description())
        
class MultiVoigtSanity(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def oneVsMulti_sanity(self):
    """
      Checking MultiVoigt1d vs. Voigt1d
    """
    from PyAstronomy import funcFit as fuf
    import numpy as np
    
    v1 = fuf.Voigt1d()
    vm = fuf.MultiVoigt1d(3)
    
    x = np.linspace(-10.,10.,200)
    
    ps = {"lin":0.0067, "off":-12.753, "A":9.0, "ad":2.17, "al":1.78, "mu":-0.771}
    
    v1.assignValues(ps)
    
    vm["lin"] = v1["lin"]
    vm["off"] = v1["off"]
    
    for i in range(1, 4):
      num = str(i)
      for p in ["ad", "al", "mu"]: 
        vm[p+num] = v1[p]
      vm["A"+num] = v1["A"]/3.0
    
    self.assertAlmostEqual(np.max(np.abs(vm.evaluate(x) - v1.evaluate(x))), 0.0, delta=1e-12, \
                           msg="MultiVoigt and Voigt1d deviate with 'amplitude separation'.")
    
    
    v1 = fuf.Voigt1d()
    vm = fuf.MultiVoigt1d(3)
    
    x = np.linspace(-10.,10.,200)
    
    ps = {"lin":0.0067, "off":-12.753, "A":9.0, "ad":2.17, "al":1.78, "mu":-0.771}
    
    v1.assignValues(ps)
    
    vm["lin"] = v1["lin"]
    vm["off"] = v1["off"]
    
    for i in range(1, 4):
      num = str(i)
      for p in ["ad", "al", "mu"]: 
        vm[p+num] = v1[p]
      if i == 1:
        vm["A"+num] = v1["A"]
        
    self.assertAlmostEqual(np.max(np.abs(vm.evaluate(x) - v1.evaluate(x))), 0.0, delta=1e-12, \
                           msg="MultiVoigt and Voigt1d deviate with one nonvanishing profile in multiVoigt.")
      
    
    