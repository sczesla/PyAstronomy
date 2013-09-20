from scipy.stats import f
import numpy
from scipy.special import erf
from scipy.optimize import fsolve
from PyAstronomy.pyaC import pyaErrors as PE

def chi2test():
  pass

def bayesFactors(model1,model2):
  raise(PE.PyAUnclassifiedError("The function 'bayesFactors' has moved to funcFit.utils."))


def ftest(chi1, chi2, dof1, dof2, compute_ratio_of_variance=False):
  """
    Performs an F-test. 
    
    Parameters
    ----------
    chi1 : float
        The chi-squared statistic or variance for model 1.
    chi2 : float
        The chi-squared statistic or variance for model 2.
        Should be better (lower) than `chi1`.
    dof1 : integer
        The number of degrees of freedom of model 1
    dof2 : integer
        The number of degrees of freedom of model 2.
        Should be lower than `dof1`.
    compute_ratio_of_variance : boolean, optional
        Determines how the F-statistics is computed.
        If False (default), the "ANOVA F-test in regression analysis
        of nested nonlinear models" is applied.
        If True, "Fisher's classical test for the equality
        of two variances" is used.
   
    Returns
    -------
    Test results : dictionary
        - "F statistic" - Value of F statistic.
        - "p-value" - Probability to obtain an F-value as large as
          or larger than observed assuming H_0 to be true
        - "Gaussian sigma level" - Same probability converted to
          Gaussian standard deviations ("sigmas")
   
   
    Notes
    -----

    The "method of least squares" is widely used in parameter estimation.
    Much of its appeal lies in its simplicity: Minimizing the sum of squared
    differences between the model and the data.

    The backbone of least squares is the classical multiple regression
    analysis using linear models to relate several independent variables to
    a dependent variable. However, many physical models are nonlinear;
    this is the realm of nonlinear least squares regression methods.
    
    The aim of model fitting is to explain as much of the variation in the
    dependent variable as possible from information contained in the
    independent variables. The contributions of the independent variables
    to the model are measured by partitions of the total sum of squares of
    the difference of data and model ("analysis of variance", ANOVA).

    **ANOVA F-test in regression analysis for nested nonlinear models**
    
    The sum of squares for any hypothesis can be determined from the
    difference between the residual sums of squares (:math:`RSS`) of two models: The
    so-called "full" model, in the context of which the null hypothesis is to be
    tested, and the "reduced" model, which is derived from
    the full model by imposing additional constraints specified by the null hypothesis;
    setting one parameter to zero is a common example.
    The reduced model is a special case of the full model and, hence, its
    residual sum of squares must not be lower than the residual
    sum of squares for the full model.
    
    Accounting for :math:`k` independent constraints on the full model and :math:`N` data
    points, we compute
    
    .. math:: Q = RSS(reduced) - RSS(full)
    
    which has :math:`k` degrees of freedom. Here, :math:`k` corresponds to the
    difference in the number of parameters---or equivalently, the difference
    in the degrees of freedoms---between the full and reduced model.
    The :math:`F`-test of the null hypothesis can
    then be computed as
    
    .. math:: \\hat{F} = (Q/k) / s^2
    
    where :math:`s^2` is an unbiased estimate of :math:`\\sigma^2`, e.g., derived from the
    full model. In the case of error-weighted data, the :math:`F` statistic reads
    
    .. math:: \\hat{F} = \\frac{(\chi^2_{reduced}-\chi^2_{full})/(\\nu_{reduced}-\\nu_{full})}{\chi^2_{full}/\\nu_{full}}
    
    where :math:`\\nu` denotes the number of degrees of freedom, which may be calculated as
    :math:`N-p-1` given a model with :math:`p` parameters and an additional constant.
    The expectation value of :math:`F` is :math:`1`, thus, if :math:`\\hat{F}=1` 
    there is no significant difference
    between the RSS of both models. If, however, :math:`F` deviates from :math:`1`, we can compute
    the probability for :math:`F` to equal or exceed the obtained value :math:`\\hat{F}` by
    
    .. math:: Prob(F\\geq \\hat{F}) = CDF(F( \\nu_{reduced}-\\nu_{full}, \\nu_{full} ))
    
    For details on the principles of the formalism,
    see Sect. 4.5.4 in Rawlings' "Applied Regression Analysis".
    
    In case of nonlinear models, the distribution of the least-square
    parameter estimates only approximately follows the normal distribution.
    Then, the :math:`F` statistic (in this case called "Wald statistic") also
    holds approximately.
    The validity of the method depends on how well the model is represented
    by a linear approximation in the parameters ("parameter effects curvature").
    For details, we refer to Sect. 15.3 in Rawlings' "Applied Regression Analysis"
    and references therein.
    
    **Fisher's classical test for the equality of two variances**
    
    If the two vying models cannot be treated as *nested*, in the sense that
    the full model encompasses the reduced model, we can directly compare
    the model variances with an F-test.
    
    Let :math:`s_1` and :math:`s_2` denote unbiased estimates of the variances
    of two independent, normally distributed populations of size :math:`N_1` and :math:`N_2`,
    respectively. Under the null hypothesis
    both variances are equal:
    
    .. math:: H_0: s_1^2=s_2^2, \hspace{20pt} H_1: s_1^2>s_2^2
    
    Then the quantity
    
    .. math:: \\hat{F} = \\frac{s_1^2}{s_2^2}
    
    is :math:`F`-distributed with :math:`N_1-1` and :math:`N_2-1` degrees of freedom.
    In the case of weighted errors we have
    
    .. math:: \\hat{F} = \\frac{\\chi_1^2/\\nu_1}{\\chi_2^2/\\nu_2}
    
    and :math:`\\hat{F}` is distributed according to an
    :math:`F`-distribution with :math:`\\nu_1` and :math:`\\nu_2`
    degrees of freedom.
    
    **References**
    
    John O. Rawlings, Sastry G. Pantula, David A. Dickey
    *Applied Regression Analysis: A Research Tool*, Second Edition
    Springer, New York, NY, 1998
    
  """
  # Some consistency checks
  if not dof2<dof1:
    raise(PE.PyAValError("The number of degrees of freedom in the new fit, dof2, should be less" + \
          " than the number of degrees of freedom in the original fit, dof1."))
  if not chi2 < chi1:
    raise(PE.PyAValError("The new fit (2) should have lower chi-squared than the original fit (1)."))

  # Compute the test statistic and the probability
  if compute_ratio_of_variance:
    F = (chi1/dof1) / (chi2/dof2)
    prob = f.cdf(F, dof1, dof2)
  else:
    F = ( (chi1-chi2)/(dof1-dof2) )   /  (chi2/dof2)
    prob = f.cdf(F, dof1-dof2, dof2)

  # Convert probability to Gaussian 'sigmas'
  def err(n,prob):
    return erf(n/numpy.sqrt(2))-prob
  sigmalvl = float(fsolve(err, 1.0, args=(prob,)))
  
  return {"F statistic": F, "p-value": 1.-prob, "Gaussian sigma level": sigmalvl}
