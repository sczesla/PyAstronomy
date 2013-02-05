import numpy

def bayesFactors(model1,model2):
  """
    Computes the Bayes factor for two competing models
    according to the prescription by Newton & Raftery 1994.
    
    Parameters
    ----------
    model1 : PyAstronomy.funcFit.TraceAnalysis instance
    model2 : PyAstronomy.funcFit.TraceAnalysis instance
    
    Returns
    -------
    The Bayes factor BF (neither log10(BF) nor ln(BF)
    
  """
  import pymc
  
  logp1 = model1["deviance"]/(-2.)
  logp2 = model2["deviance"]/(-2.)

  bf = numpy.exp(  pymc.flib.logsum(-logp1)-numpy.log(len(logp1)) 
                - (pymc.flib.logsum(-logp2) - numpy.log(len(logp2))) )
  return bf