from __future__ import print_function, division
import numpy
from PyAstronomy.funcFit.utils import ic
from PyAstronomy.pyaC import pyaErrors as PE

try:
  import pymc
except ImportError:
  pass
  
def bayesFactors(model1, model2):
  """
    Computes the Bayes factor for two competing models.
    
    The computation is based on Newton & Raftery 1994 (Eq. 13).
    
    Parameters
    ----------
    model1 : PyAstronomy.funcFit.TraceAnalysis instance
        TraceAnalysis for the first model.
    model2 : PyAstronomy.funcFit.TraceAnalysis instance
        TraceAnalysis for the second model.
    
    Returns
    -------
    BF : float
        The Bayes factor BF (neither log10(BF) nor ln(BF)).
        Note that a small value means that the first model
        is favored, i.e., BF=p(M2)/p(M1).
    
  """
  if not ic.check["pymc"]:
    raise(PE.PyARequiredImport("pymc is required to evaluate the bayesFactor", \
                               solution="Install pymc"))
  
  logp1 = model1["deviance"]/(-2.)
  logp2 = model2["deviance"]/(-2.)

  # Given a number of numbers x1, x2, ... whose logarithms are given by l_i=ln(x_i) etc.,
  # logsum calculates: ln(sum_i exp(l_i)) = ln(sum_i x_i)
  bf = numpy.exp(  pymc.flib.logsum(-logp1) - numpy.log(len(logp1)) 
                - (pymc.flib.logsum(-logp2) - numpy.log(len(logp2))) )
  return bf