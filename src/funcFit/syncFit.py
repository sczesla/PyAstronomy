from __future__ import print_function, division
import copy
import re
import numpy
from PyAstronomy.pyaC import pyaErrors as PE
from .onedfit import MiniFunc
from .params import equal
from .onedfit import _PyMCSampler, _OndeDFitParBase
from .nameIdentBase import ModelNameIdentBase
from PyAstronomy import pyaC 
import six
import six.moves as smo

from PyAstronomy.funcFit import _scoImport, _pymcImport
if _scoImport:
  import scipy.optimize as sco
if _pymcImport:
  import pymc


class MiniFuncSync:
  """
    This decorator can be applied to use
    self-defined objective functions.
    
    Applied to an objective function, it adds the functionality
    needed to evaluate the model given a certain parameter vector,
    so that the user does only have to take care
    about the quantity to be minimized.
    
    Parameters
    ----------
    odf : fitting object
        The fitting object that is supposed to use the self-defined
        objective function.
  """
  
  def __init__(self, odf):
    """
      Parameter:
        - `odf` - An instance of a fitting objects such as for example *GaussFit1d*.
    """
    # Save a REFERENCE to the fitting object
    self.odf = odf
  
  def __call__(self, f):
    """
      Parameter:
       - `f` - The user-defined objective function.
    """
    def miniFunc(P):
      # Update the parameter values in the 'Params' class instance.
      self.odf.pars.setFreeParams(P)
      # Update self.model to hold the evaluated function.
      self.odf.updateModel()
      
      val = f(self.odf, P)
      # Assign penalty
      val += self.odf.pars.getPenalty(penaltyFact=self.odf.penaltyFactor)[0]
      val += self.odf.pars.applyConditionalRestrictions()
     
      return val
    return miniFunc


class SyncFitContainer(_PyMCSampler, _OndeDFitParBase):
  
  def addComponent(self, newCompo):
    """
      Add a new component to the container.
      
      Parameters:
        - `newCompo` - A funcFit model.
      
      Returns:
        The component identifier.
    """
    # Copy component into internal model container
    # (use 1 for the first such model)
    n = len(self._compos) + 1
    self._compos[n] = copy.deepcopy(newCompo)
    # Rename the parameters (attach '_sn' with n being the component number)
    for c in self._compos[n]._compoWalk():
      if c._isComposed(): continue
      c.setRootName(c.naming.getRoot()+"[s"+str(n)+"]", rename=True)
    # Update the pars member
    if self.pars is not None:
      self.pars = self.pars + self._compos[n].pars
    else:
      self.pars = copy.deepcopy(self._compos[n].pars)
    # Assign reference to updated parameter set (in ALL models)
    for k in six.iterkeys(self._compos):
      self._compos[k].pars = self.pars
      self._compos[k]._newPars(self.pars)
    # Define the propMap (identical mapping)
    self.propMap = dict(zip(self.parameters().keys(), self.parameters().keys()))
    # Returns the component identifier
    return n
    
  def components(self):
    """
      Returns:
        A list holding the component names.
    """
    return list(self._compos.keys())
      
  def evaluate(self, axes, component=None):
    """
      Parameters:
        - `axes` - If `component` is not given, a dictionary holding the x-axis for each \
                   component name. Otherwise, the x-axis for the specified component.
        - `component` - string, optional, The name of the component to be evaluated.
      
      The evaluated model(s) is saved in the `models` dictionary.
    """
    if component is not None:
      self.models[component] = self._compos[component].evaluate(axes)
    else:
      for k, v in six.iteritems(axes):
        self.models[k] = self._compos[k].evaluate(v)
  
  def updateModel(self):
    """
      Evaluate all components. Updates the values in the `models` dictionary.
    """
    for c in six.iterkeys(self._compos):
      self.evaluate(self.data[c][0], component=c)
  
  def __chiSqr(self):
    @MiniFuncSync(self)
    def miniChiSqr(odf, P):
      # Calculate chi^2 and apply penalty if boundaries are violated.
      chi = 0.0
      for k in six.iterkeys(self._compos):
        chi += numpy.sum(((self.data[k][1] - self.models[k])/self.yerr[k])**2)
      return chi
    return miniChiSqr
  
  def __sqrDiff(self):
    @MiniFuncSync(self)
    def minisqr(odf, P):
      # Calculate squared difference
      sqr = 0.0
      for k in six.iterkeys(self._compos):
        sqr += numpy.sum((self.data[k][1] - self.models[k])**2)
      return sqr
    return minisqr

  def __cash79(self):
      @MiniFuncSync(self)
      def miniCash79(odf, P):
        # Calculate Cash statistics according to Cash 1979 (ApJ 228, 939)
        cc = 0
        for k in six.iterkeys(self._compos):
          cc += -2.0 * numpy.sum(self.data[k][1] * numpy.log(self.models[k]) - self.models[k])
        return cc
      return miniCash79
  
  def __chiSqrRobust(self):
    @MiniFuncSync(self)
    def miniChiSqr(odf, P):
      # Calculate chi^2 and apply penalty if boundaries are violated.
      chi = 0.0
      for k in six.iterkeys(self._compos):
        chi += numpy.nansum(((self.data[k][1] - self.models[k])/self.yerr[k])**2)
      return chi
    return miniChiSqr
  
  def __sqrDiffRobust(self):
    @MiniFuncSync(self)
    def minisqr(odf, P):
      # Calculate squared difference
      sqr = 0.0
      for k in six.iterkeys(self._compos):
        sqr += numpy.nansum((self.data[k][1] - self.models[k])**2)
      return sqr
    return minisqr

  def __cash79Robust(self):
      @MiniFuncSync(self)
      def miniCash79(odf, P):
        # Calculate Cash statistics according to Cash 1979 (ApJ 228, 939)
        cc = 0
        for k in six.iterkeys(self._compos):
          cc += -2.0 * numpy.nansum(self.data[k][1] * numpy.log(self.models[k]) - self.models[k])
        return cc
      return miniCash79
  
  def treatAsEqual(self, parameter):
    """
      Treat parameters as equal.
      
      `parameter` - string or list of string,
                    If a string is given, all parameters with this "base name" (i.e., neglecting
                    everything after an underscore) will be treated as equal.
                    Otherwise the specified parameters will be treated as equal.
      
      This method uses the *relations* known from *OneDFit* to treat parameters as equal.
      Dependent variables are thawed before the relation is applied, if they are not already
      free.
    """
    if isinstance(parameter, six.string_types):
      ps = []
      for p in six.iterkeys(self.parameters()):
        r = re.match(parameter+"_.*", p)
        if r is not None:
          ps.append(p)
      ps = sorted(ps)    
    else:
      ps = parameter
    for i in smo.range(1,len(ps)):
      if not ps[i] in self.freeParamNames():
        self.thaw(ps[i])
      self.relate(ps[i], [ps[0]], equal)
  
  def parameterSummary(self, toScreen=True, prefix=""):
    lines = []
    for k, v in six.iteritems(self._compos):
      lines.append(prefix)
      lines.append(prefix + "Parameters for syncFit Component: " + str(k))
      lines.append(prefix + "=" * len(lines[-1]))
      lines.extend(v.parameterSummary(toScreen=False, prefix=prefix))
    if toScreen:
      for l in lines:
        print(l)
    return lines
  def setObjectiveFunction(self, miniFunc="chisqr"):
    """
      Define the objective function.
      
      This function sets the `miniFunc` attribute, which is used
      to calculate the quantity to be minimized.
      
      Parameters
      ----------
      miniFunc : str {chisqr, cash79, sqrdiff} or callable
          The objective function. If "chisqr", chi-square will be
          minimzed. If "cash 79", the Cash statistics 
          (Cash 1979, ApJ 228, 939, Eq. 5) will be used.
          If "sqrdiff" is specified, 
          Otherwise, a user-defined function is assumed.
    """
    # Determine function to be minimized
    if miniFunc == "chisqr":
      self.miniFunc = self.__chiSqr()
      return
    elif miniFunc == "cash79":
      self.miniFunc = self.__cash79()
      return
    elif miniFunc == "sqrdiff":
      self.miniFunc = self.__sqrDiff()
      return
    elif miniFunc == "chisqrRobust":
      self.miniFunc = self.__chiSqrRobust()
      return
    elif miniFunc == "cash79Robust":
      self.miniFunc = self.__cash79Robust()
      return
    elif miniFunc == "sqrdiffRobust":
      self.miniFunc = self.__sqrDiffRobust()
      return
    else:
      if not hasattr(miniFunc, '__call__'):
        raise(PE.PyAValError("`miniFunc` is neither None, a valid string, or a function.",
                             where="OneDFit::fit",
                             solution="Use, e.g., 'chisqr' or another valid choice from the documentation."))
      
      # A function has been specified
      self.miniFunc = miniFunc
      return
  
  def fit(self, data, yerr=None, X0 = None, minAlgo=None, miniFunc=None, *fminPars, **fminArgs):
    """
      Carries out the fit.
      
      In principle, any fit algorithm can be used. If none is specified, the default is \
      scipy.optimize.fmin (Nelder-Mead Simplex). Another choice could for instance be \
      scipy.optimize.fmin_powell. After the fit, the return value of the fitting method \
      is stored in the class property `fitResult` and the `model` property is set to \
      the best fit.
      
      Parameters:
        - `data` - Dictionary of the form: {c:[x, y], ...}. Here `c` is the component number (starts
                   with one, and `x`, `y` are the x-axis and associated values. 
        - `yerr` - array, optional: Error of data values. A dictionary of the from: {c:yerr}, where
                   `c` is the component and yerr the array of error bars.
        - `X0`   - list, optional: The initial guess. If not provided, it will be assumed \
                 that self.pars already contains the initial guess.
        - `minAlgo` - callable, The minimization algorithm. Default is scipy.optimize.fmin; \
                      other algorithms from scipy may be chosen. Alternatively, any callable \
                      object taking the function to minimize as the first, the vector of starting \
                      values as the second, and a full_output flag as arguments can be used.
        - `fminArgs` - Keywords which are passed to the minimization method (default is \
                       of scipy.optimize.fmin) (e.g., `xtol` or `ftol`).
        - `fminPars` - Non-keyword arguments passed to the  minimization method (e.g., fprime in \
                       scipy.optimize.fmin_ncg).
    """
    # Assign internal data properties.
    
    if data is not None:
      self.data = data
    elif self.data is not None:
        data = self.data
    else:
        raise(PE.PyAValError("You must provide data to fit.", solution="Call fit with data."))
        
    
    if yerr is not None:
      self.yerr = yerr
    # Choose minimization algorithm
    if minAlgo is None:
      # If not specified use default.
      # Check whether it is available...
      global _scoImport
      if not _scoImport:
        raise(PE.PyARequiredImport("SciPy.optimize could not be imported.", solution="Install SciPy (see www.scipy.org/)."))
      self.minAlgo = sco.fmin
    else:
      self.minAlgo = minAlgo
    # Determine function to be minimized
    if (miniFunc is None) and (yerr is not None):
        miniFunc = "chisqr"
    elif (miniFunc is None) and (yerr is None):
        miniFunc = "sqrdiff"
    self.setObjectiveFunction(miniFunc)
    # Assign initial guess if necessary
    if X0 is not None:
      self.pars.setFreeParams(X0)
    # Save fminPars and fminArgs to internal variables
    self.fminArgs = fminArgs
    self.fminPars = fminPars
    # Carry out fit
    self.fitResult = self.minAlgo(self.miniFunc, self.pars.getFreeParams(), *self.fminPars, \
                             full_output=True, **self.fminArgs)
    self.pars.setFreeParams(self.fitResult[0])
    self.updateModel()    
    self._stepparEnabled = True

  def fitMCMC(self, data, X0, Lims, Steps, yerr=None, pymcPars=None, pyy=None, \
              potentials=None, dbfile="mcmcSample.tmp", dbArgs=None, adaptiveMetropolis=False,
              **sampleArgs):
    """
      Carry out MCMC fit/error estimation.
      
      This member is designed to provide a flexible but *easy to use* interface to \
      the capabilities of pymc. In the simplest case, it assumes a Poisson or Gaussian \
      distribution of data points and uses continuous, uniform variables (all free fitting \
      variables) with starting values defined by `X0`, Limits given by `Lims`, and
      step sizes given by `Steps` to sample from the posterior.
      
      .. note::
         The result (the Marchov-Chain/PyMC MCMC-object) will \
         be contained in the self.MCMC property; the output dictionary of MCMC.stats() \
         (Mean, HPD interval etc.) is saved to self.basicStats.
      
      Parameters:
        - `x` - An array providing the x-values of the data points.
        - `y` - An array providing the y-values of the data points. \
                Note that for the MCMC process, it is essential to know the underlying \
                distribution of the data points. *fitMCMC* assumes Poisson distributed data \
                of `yerr` is not specified and Gaussian data if it is specified. If other \
                distributions shall be used, the `pyy` parameter must contain a \
                *pymc* random variable specifying it.
        - `yerr` - array, optional,
                   Error of data values. A dictionary of the from: {c:yerr}, where
                   `c` is the component and yerr the array of error bars. If specified a Gaussian
                   distribution will be assumed for the data points, otherwise a Poisson distribution
                   is assumed.
        - `pyy` - *optional*,
                  Can be used to handle a PyMC variable containing the \
                  data. This can be useful if the distribution is neither Poisson nor \
                  Gaussian as otherwise assumed by this function.
        - `X0` - A dictionary holding {"parName":value, ...} specifying the start values. \
                 Note that parameters treated in pymcPars must not be part of this.
        - `Lims` - A dictionary of type {"ParName:[min,max], ...} specifying the lower \
                   and upper limit of a variable. \
                   Note that parameters treated in pymcPars must not be part of this.
        - `Steps` - A dictionary providing the step sizes for the MCMC \
                    sampler.
        - `pymcPars` - *optional*,
                       This variable is supposed to hold a dictionary \
                       of the form {"parName":PyMC-Variable, ...}. `pymcPars` can be \
                       used to specify a nonuniform distribution for a parameter.
        - `potentials` - *optional*,
                         Can be used to provide a list of PyMC potentials,
                         which may be needed to provide priors.
        - `dbfile` - The name of the output file, which is to hold the MCMC chain \
                     produced during sampling.
        - `**sampleArgs` - *optional*,
                           Here additional keywords can be specified, which \
                           will be handed to the *isample* member of PyMC. Most notably \
                           it is useful to specify **iter**, **burn**, and **thin**. For other \
                           possibilities see PyMC documentation.
    """
    global _pymcImport
    if not _pymcImport:
      raise(PE.PyARequiredImport("pymc package could not be imported.", solution="Install pymc (see http://code.google.com/p/pymc/"))
    # Assign mutable default parameters
    if pymcPars is None:
      pymcPars = {}
    if dbArgs is None:
      dbArgs = {}
    if potentials is None:
      potentials = []
    self.data = data
    self.yerr = yerr
    # Build up "concatenated" y-axis and yerr axis
    self.ycon = None
    self.yerrcon = None
    for k in six.iterkeys(self._compos):
      if self.ycon is None:
        self.ycon = self.data[k][1].copy()
      else:
        self.ycon = numpy.concatenate( (self.ycon, self.data[k][1]) )
      if (self.yerr is not None):
        if (self.yerrcon is None):
          self.yerrcon = self.yerr[k].copy()
        else:
          self.yerrcon = numpy.concatenate( (self.yerrcon, self.yerr[k]) )
    # Copy the pymcPars dictionary (prevents error on multiple sampler calls)
    pymcPars = pymcPars.copy()
    # Get the names of the free parameters
    freeNames = self.freeParamNames()
    print("Free parameters: ", freeNames)
    # Check whether parameter lists are complete, define default steps
    # if necessary. 
    self._dictComplete(freeNames, X0, "start values", forget=list(pymcPars))
    self._dictComplete(freeNames, Lims, "limits", forget=list(pymcPars))
    self._dictComplete(freeNames, Steps, "steps")
    
    # Define (or complete) the pymcPars dictionary by defining uniformly distributed
    # variables in the range [lim[0], lim[1]] with starting values defined by X0.
    for par in freeNames:
      if par in pymcPars: continue
      print("Using uniform distribution for parameter: ", par)
      print("  Start value: ", X0[par], ", Limits = [", Lims[par][0], ", ", Lims[par][1], "]")
      pymcPars[par] = pymc.Uniform(par, lower=Lims[par][0], upper=Lims[par][1], value=X0[par], doc="Automatically assigned parameter.")
    
    def getConcatenatedModel():
      result = None
      for k in six.iterkeys(self._compos):
        if result is None:
          result = self.models[k]
        else:
          result = numpy.concatenate( (result, self.models[k]) )
      return result
    
    # This function is used to update the model
    def getModel(**vals):
      self.assignValue(vals)
      self.updateModel()
      return getConcatenatedModel()
      
    modelDet = pymc.Deterministic(
        eval = getModel,
        name = 'Model',
        parents = pymcPars,
        doc = 'The model',
        trace = True,
        verbose = 0,
        dtype=float,
        plot=False,
        cache_depth = 2)
    
    # Define the 'data' (y-values)
    if pyy is None:
      if yerr is None:
        print("Assuming Poisson distribution for 'y'. Use 'pyy' parameter to change this!")
        pyy = pymc.Poisson("y", mu=modelDet, value=self.ycon, observed=True)
      else:
        print("Assuming Gaussian distribution for 'y'. Use 'pyy' parameter to change this!")
        pyy = pymc.Normal("y", mu=modelDet, tau=1.0/self.yerrcon**2, value=self.ycon, observed=True)

    # Add data to the Model
    Model = [pyy]
    # Add potentials (e.g., priors)
    Model.extend(potentials)
    # Add free parameters
    for v in six.itervalues(pymcPars):
      Model.append(v)
    
    # Check database arguments
    if not "dbname" in dbArgs:
      dbArgs["dbname"] = dbfile
    dbArgs = self._checkDbArgs(dbArgs)
    
    print("Using database arguments: ", dbArgs)
    self.MCMC = pymc.MCMC(Model, **dbArgs)
    
    # Tell the MCMC class to use the MH algorithm with specified step width
    if adaptiveMetropolis:
      self.MCMC.use_step_method(pymc.AdaptiveMetropolis, list(pymcPars.values()), shrink_if_necessary=True)
    else:
      for par in six.iterkeys(pymcPars):
        self.MCMC.use_step_method(pymc.Metropolis, pymcPars[par], proposal_sd=Steps[par], proposal_distribution='Normal')
          
    if not "iter" in sampleArgs:
      sampleArgs["iter"] = 2000
    if not "burn" in sampleArgs:
      sampleArgs["burn"] = 0
    if not "thin" in sampleArgs:
      sampleArgs["thin"] = 1
    
    print("Giving the following arguments to 'isample':")
    print("  ", sampleArgs)
    
    self.MCMC.isample(**sampleArgs)
    self.basicStats = self.MCMC.stats()
    self._basicStatMCMCOutput(self.basicStats)
    
    # Setting values to ``best fit values'' (lowest deviance)
    mindex = numpy.argmin(self.MCMC.trace("deviance")[:])
    for par in six.iterkeys(pymcPars):
      self[par] = self.MCMC.trace(par)[mindex]
    self.updateModel()
    self.MCMC.db.close() 

  
  def __extractFunctionValue(self, fr):
    """
      Returns the function value (e.g., chi-square).
      
      Parameters
      ----------
      fr : list
          The fit result returned by the fit method
          used by the `fit` method.
      
      Returns
      -------
      Function value : float
          For example, chi-square.
    """
    return fr[1]

  def steppar(self, pars, ranges, extractFctVal=None, quiet=False):
    """
      Allows to step a parameter through a specified range.
      
      This function steps the specified parameters through the given
      ranges. During each steps, all free parameters, except for those
      which are stepped, are fitted. The resulting contours allow
      to estimate confidence intervals.
      
      This command uses the fitting parameters specified on a call
      to the `fit` method. In particular, the same values for `x`,
      `y`, `yerr`, `minAlgo`, `miniFunc`, `fminPars`, and `fminArgs`
      are used.
      
      .. note:: You need to have carried out a fit before you can
                use `steppar`.
      
      Parameters
      ----------
      pars : string or list of strings
          The parameter(s) which are to be stepped.
      ranges : dictionary
          A dictionary mapping parameter name to range specifier.
          The latter is a list containing [lower limit, upper limit,
          no. of steps, 'lin'/'log']. The fourth entry, which
          is optional, is a string specifying whether a constant
          linear step size ('lin') or a constant logarithmic
          step size ('log') shall be used.
      quiet : boolean, optional
          If True, output will be suppressed.
      extractFctVal : callable, optional
          A function specifying how the function value is extracted
          from the fit result. If standard settings are used, the
          default of None is adequate.
      
      Returns
      -------
      Parameter steps : list
          The return value is a list of lists. Each individual list
          contains the values of the stepped parameters as the first
          entries (same order as the input `pars` list), the
          following entry is the value of the objective function
          (e.g., chi square), and the last entry is a tuple
          containing the indices of the steps of the parameter values.
          This last entry can be useful to convert the result into
          an arrow to plot, e.g., contours. 
    """
    if not self._stepparEnabled:
      raise(PE.PyAOrderError("Before you can use steppar, you must call a function, which enables its use (e.g., `fit`).", \
            solution="Call the `fit` method first and then try again."))
    if isinstance(pars, six.string_types):
      # Make it a list
      pars = [pars]
    # Check parameter consistency
    for p in pars:
      # Check existence
      tmp = self[p]
      if not p in ranges:
        raise(PE.PyAValError("There is no range for parameter: " + p, \
                             solution="Specify a range; e.g., {'xyz':[0.5,1.9,20,'lin']}"))
    # Function to extract function value from the fit result
    if extractFctVal is None:
      self._extractFctVal = self.__extractFunctionValue
    else:
      if not hasattr(extractFctVal, "__call__"):
        raise(PE.PyAValError("`extractFctVal` needs to be callable!", \
                             solution="Specify a function here or try to use None."))
      self._extractFctVal = extractFctVal
    # Set up ranges
    rs = []
    for par in pars:
      r = ranges[par]
      if len(r) > 4:
        # Use the axis as given
        rs.append(r)
        continue
      if len(r) < 4:
        # By default, use linear spacing
        mode = 'lin'
      else:
        if not isinstance(r[3], six.string_types):
          raise(PE.PyAValError("If the range has 4 entries, the fourth must be a string specifying the mode.", \
                               solution="Use either 'lin' or 'log' as the fourth entry."))
        mode = r[3]
      if mode == 'lin':
        rs.append(numpy.linspace(r[0], r[1], r[2]))
      elif mode == 'log':
        # Calculate factor
        s = numpy.power((r[1]/r[0]), 1.0/r[2])
        rs.append( r[0] * numpy.power(s, numpy.arange(r[2])) )
      else:
        raise(PE.PyAValError("Unknown mode: " + str(mode), \
                             solution="Use either 'lin' or 'log'."))
    # Save state of object
    saveObj = self.saveState()
    saveFitResult = self.fitResult
    saveModels = {}
    for k in six.iterkeys(self._compos):
      saveModels[k] = self.models[k].copy()
    # Freeze parameters, which are affected
    self.freeze(pars)
    # Store result
    result = []
    # Loop over the axes
    nli = pyaC.NestedLoop(list(map(len, rs)))
    for index in nli:
      for i, p in enumerate(pars):
        self[p] = rs[i][index[i]]
      # Fit using previous setting
      # Note that mAA is dispensable, because self.minAlgo will be a callable.
      self.fit(None, None, minAlgo=self.minAlgo, miniFunc=self.miniFunc, \
               *self.fminPars, **self.fminArgs)
      # Build up result
      ppr = []
      for par in pars:
        ppr.append(self[par])
      try:
        ppr.append(self._extractFctVal(self.fitResult))
      except Exception as e:
        PE.warn(PE.PyAValError("The call to the `extractFctVal` function failed. Using full output." + \
                               "\n  Original message: " + str(e)))
        ppr.append(self.fitResult)
      if not quiet:
        print("Result from last iteration:")
        print("  ", ppr)
      ppr.append(index)
      result.append(ppr)
    # Restore old state of object
    self.restoreState(saveObj)
    self.fitResult = saveFitResult
    for k in six.iterkeys(self._compos):
      self.models[k] = saveModels[k]
    return result
  
  def __init__(self):
    """
      Simultaneous model fitting.
      
      As an example, take a simultaneous measurement of a photometric planetary transit and
      the Rossiter-McLaughlin effect. Surely, both should be described by a subset of common
      parameters like the size of the planet and the large semi-major axis, but the
      models/measurements
      refer to quite different regimes: brightness and radial-velocity shift. This class
      can be used to carry out a fit of both simultaneously. 
      
      Attributes
      ----------
      pars : Instance of Params
          Manages the model parameters.
      models : dictionary
          A dictionary of the form component-number model; saves the evaluated
          models.
      penaltyFactor : float
          Factor used to scale the penalty imposed if parameter
          restrictions are violated.
      _compos : dictionary
          A dictionary of the form component-number model-component. The
          component number uniquely identifies every model component.
    
    """
    self._compos = {}
    self.models = {}
    self.pars = None
    self.penaltyFactor = 1e20
    self.naming = ModelNameIdentBase()
    self._stepparEnabled = False
    