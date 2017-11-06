import numpy as np
from math import factorial as fac
from PyAstronomy import pyaC
import six.moves as smo
import scipy.special as ss
import scipy.integrate as scinteg
import scipy.interpolate as sci


class PoisPost:
    
    def __init__(self):
        pass
    
    def zetaj(self, j, nQ, NB, alpha, gamma, p):
        """
        Zeta factor (Eq. 18).
        
        P(n_S = j) = zeta_j
        """
        numerator = p**(nQ - j) * pyaC.farat(nQ-j+NB+alpha-1, nQ-j+1-1) * pyaC.farat(j-gamma+1-1, j+1-1)
        denum = np.sum([p**(nQ - r) * pyaC.farat(nQ-r+NB+alpha-1, nQ-r+1-1) * pyaC.farat(r-gamma+1-1, r+1-1) for r in smo.range(nQ+1)])
        return numerator/denum
    
    def p_ns_x(self, x, nQ, NB, alpha, gamma, f):
        if (x < 0) or (x > nQ):
            return (0.0, None)
        # Limit for integration
        lims = [0.0, 4.0*nQ]
        
        def normsum(ls):
            ns = np.zeros(nQ+1)
            for i in xrange(nQ+1):
                ns[i] = ls**i / ss.gamma(i+1) * np.exp(-ls)
            return np.sum(ns)
        
        def func(ls):
            result = self.post_lam_s(ls, nQ, NB, alpha, gamma, f) * ls**x / ss.gamma(x+1) * np.exp(-ls) / normsum(ls)
            return result
        
        p = scinteg.quad(func, lims[0], lims[1])
        return p
    
    def post_lam_s(self, lams, nQ, NB, alpha, gamma, f):
        """
        Parameters
        ----------
        lams : float or array
            Wavelengths
        nQ : int
            Counts in source region
        NB : int
            Counts in BG region
        alpha : float
            BG prior
        gamma : float
            Source prior
        f: float
            SRC/BG region
        
        Returns
        -------
        Post : array
            Posterior density for lams
        """
        p = f / (1.0 + f)
        r0 = np.array([self.zetaj(j, nQ, NB, alpha, gamma, p) * lams**(j-gamma)*np.exp(-lams) / ss.gamma(j-gamma+1) for j in smo.range(nQ+1)])
        r1 = np.sum(r0, axis=0)
        return r1
    
    def p_ns(self, nQ, NB, alpha, gamma, f):
        
        # Pre-calculate posterior
        lams = np.arange(0.0001, 4.1*nQ, 0.02)
        pl = self.post_lam_s(lams, nQ, NB, alpha, gamma, f)  
        fpl = sci.interp1d(lams, pl)
        
        # Limit for integration
        lims = [0.001, 4.0*nQ]
        
        def normsum(ls):
            # exp(-ls) cancels
            ns = np.zeros(nQ+1)
            for i in xrange(nQ+1):
                ns[i] = ls**i / ss.gamma(i+1)
            return np.sum(ns)
        
        def func(ls, x):
            # exp(-ls) cancels
            result = fpl(ls) * ls**x / ss.gamma(x+1) / normsum(ls)
            return result
        
        pns = np.zeros(nQ+1)
        for x in xrange(nQ+1):
            pns[x] = scinteg.quad(func, lims[0], lims[1], args=(x,))[0]
        
        return lams, pl,  pns



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
    f < 1 indicates that the background region is larger than the source region (i.e.,
    the scaling of the background).
    In his derivation, Loredo 1990 considered different exposure times for the
    source and background regions, which is equivalent to having different region
    sizes with the same exposure time. The inverse rates are used as non-informative
    priors in the derivation.
    
    The resulting posterior density reads: 
    
    .. math::
    
        p(\\lambda_s|ns, nb) = \\sum_{i=1}^{n_s} C_i \\frac{\\lambda_s^{i-1} e^{-\\lambda_s}}{(i-1)!}
    
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
        is implicitly included here.
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
    
