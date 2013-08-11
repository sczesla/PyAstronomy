import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE

class NelderMead:
  """
    Downhill-Simplex algorithm for minimization.
    
    This implementation is based on the publication:
    Nelder and Mead, The Computer Journal 7, 308-313, 1965 (NM 1965)
    
    :Halting criterion:
    The default stop criterion is the one used by NM 1965. In
    particular, the value
    
    .. math:: \\sqrt{\\sum (\\bar{y}-y_i) / n}
    
    is calculated. If it falls below the limit defined by the
    attribute `nmCritLim`, the iteration stops.
   
    
    Attributes
    ----------
    alpha, beta, gamma : float
        The reflection-, expansion-, and contraction-coefficients.
        The default values (after NM 1965) are 1.0, 0.5, and 2.0.
        These coefficients control the modification of the simplex.
    initialStepWidthFac : float
        This factor determines how the initial simplex is calculated.
        The first simplex point is the starting value, the others are
        constructed by adding a fraction defined by this factor to
        the starting value. The default is 0.05.
    nmCritLim : float
        Critical value for the NM 1965 stopping criterion. The
        default is 1e-8.
    _maxIter : int
        The maximum number of iterations. The default is 10000.
        
  """
  
  def __init__(self):
    # Best guesses after NM 1965
    self.alpha, self.beta, self.gamma = (1.0, 0.5, 2.0)
    # Step-width factor for initial simplex (if not specified otherwise)
    self.initialStepWidthFac = 0.05
    # Critical limit for the NM 1965 stopping criterion
    self.nmCritLim = 1e-8
    # Stopping criterion
    self._stopCrit = self._stopNM1965
    # Maximum number of iterations
    self._maxIter = 10000
  
  def _initSimplex(self, m, initDelta):
    """
      Define the initial simplex.
    """
    # The simplex: n+1 points in n-dimensional space
    self._simplex = np.zeros( (self._n+1, self._n) )
    # The function values
    self._yi = np.zeros(self._n+1)
    # Define first point (order of parameters guaranteed)
    self._simplex[0,::] = np.array(m.pars.getFreeParams())
    self._yi[0] = m.miniFunc(self._simplex[0,::])
    for i, p in enumerate(self._fpns):
      self._simplex[i+1,::] = self._simplex[0,::].copy()
      if (initDelta is not None) and (p in initDelta):
        self._simplex[i+1, i] += initDelta[p]
      else:
        if self._simplex[0,i] == 0.0:
          PE.PyAValError("The start value of the parameter '" + p + "' is 0.0. This will lead" + \
                         "to an appropriate simplex.", \
                          solution="Specify an initial step via the `initDelta` parameter.\n" + \
                          "The size of the step should reflect the 'scale' of the problem.")
        self._simplex[i+1, i] += self._simplex[0,i] * self.initialStepWidthFac
      self._yi[i+1] = m.miniFunc(self._simplex[i+1,::])
    
  def _step(self, m):
    """
      Proceed one step of the simplex algorithm.
    """
    # Indices with lowest and highest objective function value
    h = np.argmax(self._yi)
    l = np.argmin(self._yi)
    # Barycenter (excluding the worst point)
    pb = (np.sum(self._simplex, 0) - self._simplex[h,::]) / self._n
    # New suggestion
    ps = (1. + self.alpha) * pb - self.alpha*self._simplex[h,::]
    ys = m.miniFunc(ps)
    if ys < self._yi[l]:
      # In this case, calculate pss
      pss = self.gamma*ps + (1. - self.gamma) * pb
      yss = m.miniFunc(pss)
      if yss < self._yi[l]:
        # Replace by pss
        self._simplex[h,::] = pss
        self._yi[h] = yss
        return
      # Replace h by ps
      self._simplex[h,::] = ps
      self._yi[h] = ys
      return
    # We have not produced a new minimum, and we know
    # that ys >= y[l]
    if ys < self._yi[h]:
      self._simplex[h,::] = ps
      self._yi[h] = ys
      return
    else:
      # There is a new maximum
      pssb = self.beta*self._simplex[h,::] + (1.0 - self.beta)*pb
      yssb = m.miniFunc(pssb)
      if yssb > self._yi[h]:
        # Contract
        pl = self._simplex[l,::].copy()
        for i in xrange(self._n+1):
          if i == l: continue
          self._simplex[i,::] = (self._simplex[i,::] + pl)/2.0
          self._yi[i] = m.miniFunc(self._simplex[i,::])
        return
      else:
        self._simplex[h,::] = pssb
        self._yi[h] = yssb
        return
    raise(PE.PyAOrderError("Nothing done in Nelder-Mead step. If you see this, there is a coding error."))
  
  def _stopNM1965(self):
    """
      Apply the NM 1965 stopping criterion.
      
      Returns
      -------
      Reached : boolean
          True, if the stopping criterion has been reached.
    """
    ym = np.mean(self._yi)
    return (np.sqrt(np.sum((self._yi-ym)**2))/self._n) < self.nmCritLim
  
  def fit(self, m, ds, objf="chisqr", initDelta=None):
    # Number of free parameters
    self._n = m.numberOfFreeParams()
    # Set objective function
    m.setObjectiveFunction(objf)
    # Assign data object
    m._fufDS = ds
    # Names of free parameters (order guaranteed)
    self._fpns = m.freeParamNames()
    # Initial simplex
    self._initSimplex(m, initDelta)
    
    self.iterCount = 0
    while (not self._stopCrit()) and (self.iterCount < self._maxIter):
      self._step(m)
      self.iterCount += 1
    
    # Find the optimum parameter set
    l = np.argmin(self._yi)
    m.pars.setFreeParams(self._simplex[l,::])
    # Evaluate model so that model attribute holds the best match
    m.evaluate(ds.x)
    
    if self.iterCount == self._maxIter:
      print "Maxiter reached"
    