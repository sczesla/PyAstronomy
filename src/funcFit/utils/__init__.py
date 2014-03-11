from PyAstronomy.pyaC import ImportCheck

ic = ImportCheck(["numpy", "scipy", "pymc", "matplotlib", "pyfits", "matplotlib.pylab"])

if ic.check["pymc"] and ic.check["numpy"] and ic.check["scipy"] and ic.check["matplotlib"] and ic.check["matplotlib.pylab"]:
  from anaMCMCTraces import TraceAnalysis
  from bayesFactors import bayesFactors
  __all__ = ["TraceAnalysis", "bayesFactors"]