from PyAstronomy.funcFit import ic

if ic.check["numpy"] and ic.check["scipy"]:
    from .anaMCMCTraces import TraceAnalysis, hpd, quantiles
    from .bayesFactors import bayesFactors
    __all__ = ["TraceAnalysis", "bayesFactors", "hpd", "quantiles"]