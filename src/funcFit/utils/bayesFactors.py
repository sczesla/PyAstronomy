from __future__ import print_function, division
import numpy as np
from PyAstronomy.funcFit.utils import ic
from PyAstronomy.pyaC import pyaErrors as PE
import scipy.special as ss

def pD_pr2(logl):
    """
    Estimate the marginal probability of the data
    
    This function estimates the marginal probability of the data using the harmonic
    mean of the likelihood (KR95, Eq. 11).
    
    The resulting estimate is known to be unstable and is often dominated by extreme points.
    Although the estimate can be shown to converge to the "correct" value, its distribution
    does not generally correspond to a Gaussian shrinking in width as the number of points
    increases (Gaussian central limit theorem). See KR95 for further discussion.
    
    KR95: Kass, Robert E., and Adrian E. Raftery. "Bayes Factors." Journal of the American Statistical Association 90, no. 430
    (1995): 773-95. doi:10.2307/2291091.
    
    Parameters
    ----------
    logl : array
        The (natural) logarithm of the likelihood (e.g., obtained from a Markov Chain)
    
    Returns
    -------
    ln(pD) : floats
        The (natural) logarithm of the estimation of the marginal
        probability of the data
    """
    return np.log(len(logl)) - ss.logsumexp(-logl)


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
    logp1 = model1["deviance"] / (-2.)
    logp2 = model2["deviance"] / (-2.)

    # Given a number of numbers x1, x2, ... whose logarithms are given by l_i=ln(x_i) etc.,
    # logsum calculates: ln(sum_i exp(l_i)) = ln(sum_i x_i)
    bf = np.exp(ss.logsumexp(-logp1) - np.log(len(logp1))
                   - (ss.logsumexp(-logp2) - np.log(len(logp2))))
    return bf
