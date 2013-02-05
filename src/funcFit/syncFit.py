import copy
import re
import numpy
from PyAstronomy.pyaC import pyaErrors as PE
from onedfit import MiniFunc
from params import equal
from onedfit import _PyMCSampler, _OndeDFitParBase
from nameIdentBase import ModelNameIdentBase

from PyAstronomy.funcFit import _scoImport, _pymcImport
if _scoImport:
  import scipy.optimize as sco
if _pymcImport:
  import pymc


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
    for k in self._compos.iterkeys():
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
    return self._compos.keys()
      
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
      for k, v in axes.iteritems():
        self.models[k] = self._compos[k].evaluate(v)
  
  def updateModel(self):
    """
      Evaluate all components. Updates the values in the `models` dictionary.
    """
    for c in self._compos.iterkeys():
      self.evaluate(self.data[c][0], component=c)
  
  def __chiSqr(self):
    @MiniFunc(self)
    def miniChiSqr(odf, P):
      # Calculate chi^2 and apply penalty if boundaries are violated.
      chi = 0.0
      for k in self._compos.iterkeys():
        chi += numpy.sum(((self.data[k][1] - self.models[k])/self.yerr[k])**2)
      return chi
    return miniChiSqr
  
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
    if isinstance(parameter, basestring):
      ps = []
      for p in self.parameters().iterkeys():
        r = re.match(parameter+"_.*", p)
        if r is not None:
          ps.append(p)
    else:
      ps = parameter
    for i in xrange(1,len(ps)):
      if not ps[i] in self.freeParamNames():
        self.thaw(ps[i])
      self.relate(ps[i], [ps[0]], equal)
  
  def parameterSummary(self, toScreen=True, prefix=""):
    lines = []
    for k, v in self._compos.iteritems():
      lines.append(prefix)
      lines.append(prefix + "Parameters for syncFit Component: " + str(k))
      lines.append(prefix + "=" * len(lines[-1]))
      lines.extend(v.parameterSummary(toScreen=False, prefix=prefix))
    if toScreen:
      for l in lines:
        print l
    return lines
  
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
    self.data = data
    if yerr is not None:
      self.yerr = yerr
    else:
      self.yerr = {}
      for k in self._compos.iterkeys():
        self.yerr[k] = numpy.ones(len(self.data[k][0]))
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
    if miniFunc is None:
      self.miniFunc = self.__chiSqr()
    else:
      self.miniFunc = miniFunc
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

  def fitMCMC(self, data, X0, Lims, Steps, yerr=None, pymcPars={}, pyy=None, \
              potentials=[], dbfile="mcmcSample.tmp", dbArgs={}, adaptiveMetropolis=False,
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
    self.data = data
    self.yerr = yerr
    # Build up "concatenated" y-axis and yerr axis
    self.ycon = None
    self.yerrcon = None
    for k in self._compos.iterkeys():
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
    print "Free parameters: ", freeNames
    # Check whether parameter lists are complete, define default steps
    # if necessary. 
    self._dictComplete(freeNames, X0, "start values", forget=pymcPars.keys())
    self._dictComplete(freeNames, Lims, "limits", forget=pymcPars.keys())
    self._dictComplete(freeNames, Steps, "steps")
    
    # Define (or complete) the pymcPars dictionary by defining uniformly distributed
    # variables in the range [lim[0], lim[1]] with starting values defined by X0.
    for par in freeNames:
      if par in pymcPars: continue
      print "Using uniform distribution for parameter: ", par
      print "  Start value: ", X0[par], ", Limits = [", Lims[par][0], ", ", Lims[par][1], "]" 
      pymcPars[par] = pymc.Uniform(par, lower=Lims[par][0], upper=Lims[par][1], value=X0[par], doc="Automatically assigned parameter.")
    
    def getConcatenatedModel():
      result = None
      for k in self._compos.iterkeys():
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
        print "Assuming Poisson distribution for 'y'. Use 'pyy' parameter to change this!"
        pyy = pymc.Poisson("y", mu=modelDet, value=self.ycon, observed=True)
      else:
        print "Assuming Gaussian distribution for 'y'. Use 'pyy' parameter to change this!"
        pyy = pymc.Normal("y", mu=modelDet, tau=1.0/self.yerrcon**2, value=self.ycon, observed=True)

    # Add data to the Model
    Model = [pyy]
    # Add potentials (e.g., priors)
    Model.extend(potentials)
    # Add free parameters
    for v in pymcPars.itervalues():
      Model.append(v)
    
    # Check database arguments
    if not "dbname" in dbArgs:
      dbArgs["dbname"] = dbfile
    dbArgs = self._checkDbArgs(dbArgs)
    
    print "Using database arguments: ", dbArgs
    self.MCMC = pymc.MCMC(Model, **dbArgs)
    
    # Tell the MCMC class to use the MH algorithm with specified step width
    if adaptiveMetropolis:
      self.MCMC.use_step_method(pymc.AdaptiveMetropolis, pymcPars.values(), shrink_if_necessary=True)
    else:
      for par in pymcPars.keys():
        self.MCMC.use_step_method(pymc.Metropolis, pymcPars[par], proposal_sd=Steps[par], proposal_distribution='Normal')
          
    if not "iter" in sampleArgs:
      sampleArgs["iter"] = 2000
    if not "burn" in sampleArgs:
      sampleArgs["burn"] = 0
    if not "thin" in sampleArgs:
      sampleArgs["thin"] = 1
    
    print "Giving the following arguments to 'isample':"
    print "  ", sampleArgs
    
    self.MCMC.isample(**sampleArgs)
    self.basicStats = self.MCMC.stats()
    self._basicStatMCMCOutput(self.basicStats)
    
    # Setting values to ``best fit values'' (lowest deviance)
    mindex = numpy.argmin(self.MCMC.trace("deviance")[:])
    for par in pymcPars.iterkeys():
      self[par] = self.MCMC.trace(par)[mindex]
    self.updateModel()
    self.MCMC.db.close() 
  
  def __init__(self):
    """
      This class can be used in the case that two related models defined on different axes
      are to be fitted simultaneously.
      
      As an example, take a simultaneous measurement of a photometric planetary transit and
      the Rossiter-McLaughlin effect. Surely, both should be described by a subset of common
      parameters like size of the planet and large semi-major axis, but the models/measurements
      refer to quite different regimes: brightness and radial-velocity shift. This class
      can be used to carry out a fit of both simultaneously. 
      
      :Class properties:
        - `_compos` - A dictionary of the form {component-number:model-component}. The \
                      component number identifies every model component.
        - `models` - A dictionary of the form {component-number:model}; saves the evaluated \
                     model for every model component. This property is updated on a call \
                     to *updateModel*.
        - `pars` - Instance of the *Params* class to manage the model parameters.
        - `penaltyFactor` - float, Factor used to scale the penalty imposed if parameter \
                            restrictions are violated.  
    """
    self._compos = {}
    self.models = {}
    self.pars = None
    self.penaltyFactor = 1e20
    self.naming = ModelNameIdentBase()