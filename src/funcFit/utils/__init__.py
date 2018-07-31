from PyAstronomy.funcFit import ic

if ic.check["numpy"] and ic.check["scipy"] and ic.check["matplotlib"] and ic.check["matplotlib.pylab"]:
  from .anaMCMCTraces import TraceAnalysis, hpd, quantiles
  from .bayesFactors import bayesFactors, pD_pr2
  __all__ = ["TraceAnalysis", "bayesFactors", "hpd", "quantiles", "pD_pr2"]