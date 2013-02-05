import unittest
import numpy
from PyAstronomy import funcFit as fuf

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
    self.assertEquals(numpy.sum(gf.parameters().values()), 0.0)
  
  def combine1_sanity(self):
    gf = fuf.GaussFit1d()
    gff = gf + gf + gf
    for p in gff.parameters().keys():
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