import numpy as np
from math import factorial as fac
from PyAstronomy import pyaC
import six.moves as smo


def _cirob(i, ns, nb, f):
    """
    Compute Ci factor.
    
    Ci is calculated via its inverse, making the calculation
    numerically more robust.
    
    Parameters
    ----------
    i : int
        Parameter for which to calculate Ci.
    ns : int
        The number of recorded events in the source region.
    nb : int
        The number of recorded events in the background region.
    f : float
        The ratio of the areas of source and background region.
    
    Returns
    -------
    Ci : float
        C_i factor
    """
    # Compute one over ci
    oci = 0.0
    
    g = (1.+1./f)
    
    for j in smo.range(1, ns+1):
        oci += g**(j-i) * pyaC.farat(ns+nb-j-1, ns+nb-i-1) * pyaC.farat(ns-i, ns-j)
    
    return 1./oci


def ps_pdf_lams(lam_s, ns, nb, f):
    """
    The posterior probability density for Poisson source with background.
    
    The posterior probability density implemented here was given by
    Loredo in his treatise "From Laplace to Supernova SN 1987A: Bayesian
    Inference in Astrophysics" published in "Maximum Entropy and Bayesian
    Methods" (p. 81-142).
    
    The situation is such that nb background events are recorded in a
    a background region and ns events in a source region. Both regions are
    exposed for a certain time T, and the ratio between their areas is f, where
    f > 1 indicates that the background region is larger than the source region.
    In his derivation, Loredo 1990 considered different exposure times for the
    source and background regions, which is equivalent to having different region
    sizes with the same exposure time. The inverse rates are used as non-informative
    priors in the derivation.
    
    The resulting posterior density reads: 
    
    .. math::
    
        p(s|nb) = \\sum_{i=1}^{n_s} C_i \\frac{\\lambda_s^{i-1} e^{-\\lambda_s}}{(i-1)!}
    
    where
    
    .. math::
    
        C_i = \\frac{(1+1/f)^i \, \\frac{(ns+nb-i-1)!}{(ns-i)!}}{\\sum_{j=1}^{ns}(1+1/f)^j \, \\frac{(ns+nb-j-1)!}{(n-j)!}}
    
    Note that source region exposure time 't' used by Loredo is taken to be one numerically.
    Thus, rates refer to the source region exposure time.
    
    Parameters
    ----------
    lam_s : float or array
        The source count rate(s) at which to evaluate the posterior. The rate
        refers to the exposure time of the source and background regions, which
        is explicitly included here.
    ns : int
        The number of recorded events in the source region.
    nb : int
        The number of recorded events in the background region.
    f : float
        The ratio of the areas of source and background region.
    
    Returns
    -------
    pdf : float or array
        The posterior probability density.
    """
    r = 0.0
    for i in smo.range(1, ns+1):
        r += _cirob(i, ns, nb, f) * (lam_s)**(i-1) * np.exp(-lam_s) / fac(i-1)
    return r
    
