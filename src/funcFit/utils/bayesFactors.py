from __future__ import print_function, division
import numpy as np
from PyAstronomy.funcFit.utils import ic
from PyAstronomy.pyaC import pyaErrors as PE
import scipy.special as ss
import scipy.stats as stats
import six.moves as smo


def pD_pr5(logpost, logftheta):
    """
    Estimate the marginal likelihood of the data
    
    This function estimates the marginal likelihood of the data using a modified version of
    the harmonic mean of the likelihood estimator (KR95, Eq. 12); note that the actual estimate
    is the inverse of the expression given in KR95, Eq. 12.
    
    The choice of the weighting density f is non-trivial. KR95 argue that high efficiency is
    likely attained if f is roughly proportional to the posterior density. Depending on the
    choice of f, results may be very good or poor.
    
    KR95: Kass, Robert E., and Adrian E. Raftery. "Bayes Factors." Journal of the American Statistical Association 90, no. 430
    (1995): 773-95. doi:10.2307/2291091.
    
    Parameters
    ----------
    logPost : array
        Natural logarithm of the posterior probability used for estimation, e.g., computed
        using a Markov Chain.
    ftheta : array
        Natural logarithm of the weighting density (f) corresponding to the respective
        value of the posterior probability.
        
    Returns
    -------
    ln(pD) : float
        Natural logarithm of the marginal likelihood of the data.
    """
    return -(-np.log(len(logpost)) + ss.logsumexp(logftheta-logpost))
    
    
def pD_pr5_mvg(ta, cs=1.0, cm=None, mean=None):
    """
    Use pD_pr5 to estimate marginal likelihood of prior using multivariate Gaussian is weighting density
    
    Use a multivariate Gaussian with covariance matrix and mean estimated from the posterior samples
    to estimate the marginal likelihood of the data.
    
    KR95: Kass, Robert E., and Adrian E. Raftery. "Bayes Factors." Journal of the American Statistical Association 90, no. 430
    (1995): 773-95. doi:10.2307/2291091.
    
    Parameters
    ----------
    ta : TraceAnalysis object
        Holds the chains, parameter names, etc.
    cm : (n x n) array, optional
        The covariance matrix (n denotes the number of chains). If not given,
        the `correlationMatrix` method of TraceAnalysis is used to estimate it.
        The order of parameters is defined by the result of `availableParameters`
        if TraceAnalysis.
    cs : float, optional
        To define the weighting density, the covariance matrix is scaled by this
        factor. The default is 1.
    mean : 1-d array of length n, optional
        The means to be used for the multivariate Gaussian (see cm for comments on
        the order of parameters).
    
    Returns
    -------
    ln(pD) : float
        Natural logarithm of the marginal likelihood of the data.
    mean : array
        The means used in defining the multivariate Gaussian
    cm : array
        The covariance matrix (scaled by `cs`) used in the definition of
        the multivariate Gaussian.
    """
    aps = ta.availableParameters()
    N = len(ta[aps[0]])
    if cm is None:
        cm = ta.correlationMatrix(covariance=True, toScreen=False)[1]
    if mean is None:
        mean = np.array([np.mean(ta[p]) for p in aps])
        
    mvn = stats.multivariate_normal(mean, cm*cs)
    
    lnf = np.zeros(N)
    for i in smo.range(N):
        lnf[i] = mvn.logpdf(np.array([ta[p][i] for p in aps]))
    
    return pD_pr5(ta["lnpost"], lnf), mean, cm*cs
    


def pD_pr2(logl):
    """
    Estimate the marginal likelihood of the data using the Harmonic Mean Estimator (HME)
    
    This function estimates the marginal likelihood of the data using the harmonic
    mean of the likelihood (KR95, Eq. 11).
    
    The resulting estimate is known to be unstable and is often dominated by extreme points.
    Although the estimate can be shown to converge to the "correct" value, its distribution
    does not generally correspond to a Gaussian shrinking in width as the number of points
    increases (Gaussian central limit theorem). In fact, its convergence can be shown to be
    exceedingly slow when the prior is less informative than the likelihood, which is the
    case in almost all practical applications (see Wolpert and Schmidler 2012, Statistica
    Sinica 22, 1233-1251). Therefore, the HME should be used with care.
    
    KR95: Kass, Robert E., and Adrian E. Raftery. "Bayes Factors." Journal of the American
    Statistical Association 90, no. 430 (1995): 773-95. doi:10.2307/2291091.
    
    Parameters
    ----------
    logl : array
        The (natural) logarithm of the likelihood (e.g., obtained from a Markov Chain)
    
    Returns
    -------
    ln(pD) : floats
        The (natural) logarithm of the estimation of the marginal
        likelihood of the data
    """
    return np.log(len(logl)) - ss.logsumexp(-logl)


def bayesFactors(model1, model2):
    """
    Computes the Bayes factor for two competing models.

    The computation is based on the Harmonic Mean Estimator (Newton & Raftery 1994,
    Eq. 13). Note that the convergence of this estimator can be extremely poor (see
    documentation of pD_pr2).

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
    logp1 = model1["deviance"] / (-2.)
    logp2 = model2["deviance"] / (-2.)

    # Given a number of numbers x1, x2, ... whose logarithms are given by l_i=ln(x_i) etc.,
    # logsum calculates: ln(sum_i exp(l_i)) = ln(sum_i x_i)
    bf = np.exp(ss.logsumexp(-logp1) - np.log(len(logp1))
                   - (ss.logsumexp(-logp2) - np.log(len(logp2))))
    return bf
