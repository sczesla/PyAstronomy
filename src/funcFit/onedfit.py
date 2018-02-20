# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
from .params import Params
import re
import copy
import types
from PyAstronomy.pyaC import pyaErrors as PE
from .nameIdentBase import ModelNameIdentBase
from PyAstronomy import pyaC
from time import time as timestamp
from .fufDS import FufDS
from .extFitter import NelderMead
import six
import six.moves as smo

from PyAstronomy.funcFit import _pymcImport, _scoImport, ic

if _pymcImport:
    import pymc
if _scoImport:
    import scipy.optimize as sco
if ic.check["emcee"]:
    import emcee
if ic.check["progressbar"]:
    import progressbar


def addEval(self, x):
    return (self.leftCompo.evaluate(x) + self.rightCompo.evaluate(x))


def subEval(self, x):
    return (self.leftCompo.evaluate(x) - self.rightCompo.evaluate(x))


def divEval(self, x):
    return (self.leftCompo.evaluate(x) / self.rightCompo.evaluate(x))


def mulEval(self, x):
    return (self.leftCompo.evaluate(x) * self.rightCompo.evaluate(x))


def powEval(self, x):
    return (self.leftCompo.evaluate(x) ** self.rightCompo.evaluate(x))


class MiniFunc:
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
            # Obtain penalties
            pfac, pdict = self.odf.pars.getPenalty(
                penaltyFact=self.odf.penaltyFactor)
            # Assign penalty
            val = pfac
            # Apply conditional restrictions
            val += self.odf.pars.applyConditionalRestrictions()
            # Assign x, y, and yerr attributes. This is a not-so-beautiful
            # way of not breaking the API and having the funcFit data
            # storage object.
            self.odf.x = self.odf._fufDS.x
            self.odf.y = self.odf._fufDS.y
            self.odf.yerr = self.odf._fufDS.yerr

            try:
                # Update self.model to hold the evaluated function.
                self.odf.updateModel()
            except Exception as e:
                if val > 0.0:
                    # Allow model evaluation to fail if parameters wandered outside the
                    # valid range (returning immediately maybe an option, but this largely maintains
                    # behavior, preventing annoying errors)
                    return val
                # If no penalty is present, raise exception
                raise(PE.PyAAlgorithmFailure("Could not evaluate model for parameters: " + str(self.odf.parameters()), \
                                             where="ondeDFit", \
                                             tbfe=e, \
                                             solution=["Try to define 'restrictions' via setRestriction if parameter values are invalid.", \
                                                       "Adjust implementation of model to prevent error."]))
            # Add value of actual objective function
            val += f(self.odf, P)
            return val
        return miniFunc


class _PyMCSampler:
    """
    This class encapsulates a number of methods helping to set up the
    PyMC sampler.
    """

    def _dictComplete(self, l, d, whichDict, forget=[]):
        """
        Checks whether the list `l` contains all keys present in dictionary `d`.

        Parameters:
          `l` - A list of strings
          `d` - A dictionary with string keys.
          `whichDict` - string, Which dictionary (e.g., start values) is under consideration?
          `forget` - A list of string specifying keys, which can be omitted in d.
        """
        if len(l) != (len(d) + len(forget)):
            message = "Dictionary of " + whichDict + \
                " has not the correct number of entries - "
            if len(l) < (len(d) + len(forget)):
                message += " there are too many.\n"
            else:
                message += " there are too few.\n"
            message += "  Needed entries: " + ', '.join(l) + "\n"
            given = list(d)
            given.extend(forget)
            message += "  Given entries: " + ', '.join(given)
            raise(PE.PyAValError(message, where="fitMCMC",
                                 solution="Adjust input dictionary."))
        for elem in l:
            if elem in forget:
                continue
            if elem not in d:
                message = "Error in " + whichDict + " dictionary.\n"
                message += "  Key " + elem + " is missing!"
                raise(PE.PyAValError(message, where="fitMCMC",
                                     solution="Adjust input dictionary."))

    def _basicStatMCMCOutput(self, bs):
        print("Basic statistics of MCMC analysis: ")
        print("-----------------------------------------------------")
        for k in six.iterkeys(bs):
            print("Parameter: ", k)
            if bs[k] is None:
                print("  No output available!")
                continue
            for k2, v2 in six.iteritems(bs[k]):
                if k2 != "quantiles":
                    print("  " + k2 + " : " + str(v2))
                else:
                    s = "  quantiles : "
                    for qk, qv in six.iteritems(bs[k]["quantiles"]):
                        print(s + str(qk) + " : " + str(qv))
                        s = "              "
            print("-----------------------------------------------------")

    def _checkDbArgs(self, dbArgs):
        """
        Check whether database arguments are given. If not, use default parameters.
        """
        # If no db is specified, use pickle
        if not "db" in dbArgs:
            dbArgs["db"] = "pickle"

        if dbArgs["db"] == "pickle":
            if not "dbname" in dbArgs:
                dbArgs["dbname"] = "tmp.pickle"
            return dbArgs
        elif dbArgs["db"] == "hdf5":
            if not "dbname" in dbArgs:
                dbArgs["dbname"] = "tmp.zlib"
            if not "dbmode" in dbArgs:
                dbArgs["dbmode"] = "w"
            if not "dbcomplevel" in dbArgs:
                dbArgs["dbcomplevel"] = 9
            if not "dbcomplib" in dbArgs:
                dbArgs["dbcomplib"] = "zlib"
            return dbArgs
        else:
            raise(PE.PyAValError("Database: " + str(dbArgs["db"] + " currently not supported.",
                                                    solution="Use another database (e.g., pickle or hdf5).")))
        return dbArgs

    def MCMCautoParameters(self, ranges, picky=True, stepsize=1e-2, setRestrictionsFromPriors=False):
        """
        Convenience function to generate parameters for MCMC fit.

        This function prepares the `X0`, `Lims`, and `Steps` dictionaries
        needed to call :py:func:`fitMCMC`.

        For `X0`, the current parameter values are used. `Lims` is constructed
        using the `ranges` parameter, and `Steps` is defined on the basis
        of the `stepsize` and `ranges` parameters.

        The `picky` parameter determines how this functions handles behaves
        if it encounters parameters in the `range` dictionary, which were either
        not thawed or have been thawed but have not been specified in `ranges`.
        If `picky` is True (the default), this function throws an error if
        `ranges` does not cover all and only the free parameters. If picky is
        False, the function will automatically thaw all parameters specified
        through `ranges` and freeze the rest.

        .. warning:: There is NO guarantee that the sampling parameters
                     (start values, limits for the uniform priors, and initial
                     step sizes for the sampler) are reasonable. You need to
                     check the results.

        Parameters
        ----------
        ranges : dictionary
          Holds the fit ranges for the individual parameters. If single values 
          are given, the sampling range (uniform prior) will be arranged
          symmetrically around the parameter's current value.
          It is also possible to specify a range
          directly using, e.g., "A1":[0,100].
        stepsize : float, optional
          Defines the step size as a fraction of the fit range
          given in `ranges`.
        picky : boolean, optional
          If True (default), the list of free parameters has to match exactly
          the list of parameters specified in `ranges`. If False, the list
          of free parameters will be adapted to those given in `ranges`.
        setRestrictionsFromPriors : boolean, optional 
          Default: False. If True, parameter restrictions are applied according to
          the ranges of the uniform priors.

        Returns
        -------
        X0 : dictionary
            Maps parameter name to start value.
        lims : dictionary
            Maps parameter name to [lower, upper] limit.
        steps : dictionary
            Maps parameter name to step size. 

        Examples
        --------

        ::

          from PyAstronomy import funcFit as fuf
          import numpy as np
          import matplotlib.pylab as plt

          x = np.linspace(0,30,1000)
          gauss = fuf.GaussFit1d()
          gauss["A"] = 1
          gauss["mu"] = 23.
          gauss["sig"] = 0.5
          yerr = np.random.normal(0., 0.05, len(x))
          y = gauss.evaluate(x) + yerr
          # This step is not necessary if <picky>=False in MCMCautoParameters.
          gauss.thaw(["A","mu","sig"])
          X0, lims, steps = gauss.MCMCautoParameters({"A":[0,10],"mu":3, "sig":[0.1,1.0]})
          gauss.fitMCMC(x, y, X0, lims, steps, yerr=yerr, iter=1000)

          plt.plot(x, y, 'k+')
          plt.plot(x, gauss.evaluate(x), 'r--')
          plt.show()

        """

        X0 = {}
        steps = {}
        lims = {}

        # parameters which are free but not given
        missingParams = []
        # parameters which are not free but given
        extraParams = []

        # Check if ranges are given for all free parameters
        for p in self.freeParameters():
            if p not in ranges:
                if picky == False:
                    self.freeze(p)
                missingParams.append(p)
        # Check if all given ranges pertain to free parameters
        for p in ranges:
            if p not in self.freeParameters():
                if picky == False:
                    self.thaw(p)
                extraParams.append(p)
        # Raise exception if missing or extra parameters exits unless picky=False
        if picky:
            if len(missingParams) > 0:
                raise(PE.PyAValError("Not enough parameters: " + str(missingParams), where="_PyMCSampler::MCMCautoParameters",
                                     why="Ranges for these free fit parameters are missing."))
            if len(extraParams) > 0:
                raise(PE.PyAValError("Too many parameters: " + str(extraParams), where="_PyMCSampler::MCMCautoParameters",
                                     why="Ranges were given for non-free parameters."))

        if stepsize >= 1.0:
            raise(PE.PyAValError("The step size is larger than one. Must be a fraction of the range, i.e., smaller than 1.",
                                 solution="Redefine the step size. Use, e.g., 0.01."))

        # Loop through given parameters
        for p in ranges:
            X0[p] = self[p]
            # range is [a, b] or (a, b) then use this range, otherwise symmetric around
            #    the current value
            if isinstance(ranges[p], (list, tuple)):
                lims[p] = [ranges[p][0], ranges[p][1]]
            else:
                lims[p] = [X0[p] - ranges[p] / 2., X0[p] + ranges[p] / 2.]
            steps[p] = (max(lims[p]) - min(lims[p])) * stepsize
            if setRestrictionsFromPriors:
                self.setRestriction({p: lims[p]})

        return X0, lims, steps

    def autoFitMCMC(self, x, y, ranges, picky=True, stepsize=1e-3, yerr=None,
                    pymcPars={}, pyy=None, potentials=[], dbfile="mcmcSample.tmp",
                    dbArgs={}, **sampleArgs):
        """
        Convenience function to using auto-generated sampling parameters in MCMC.

        This function is essentially a wrapper around :py:func:`fitMCMC`.
        It allows you to use the 
        start values, step sizes, and limits constructed using
        :py:func:`MCMCautoParameters`.
        This method takes all parameters but `X0`, `Lims`, and `Steps`,
        which are expected by :py:func:`fitMCMC`. Additionally,
        the following parameters are available.

        .. warning:: There is NO guarantee that the sampling parameters
                     (start values, limits for the uniform priors, and initial
                     step sizes for the sampler) are reasonable. You need to
                     check the results.

        Parameters
        ----------
        ranges : dictionary
          Holds the fit ranges for the individual parameters. If a single values 
          is given, the fit range will be arranged symmetrically around the 
          parameter's current value. It is also possible to specify a range
          directly using, e.g., "A1":[0,100].
        stepsize : float, optional
          Defines the step size as a fraction of the fit range
          given in `ranges`.
        picky : boolean, optional
          If True (default), the list of free parameters has to match exactly
          the list of parameters specified in `ranges`. If False, the list
          of free parameters will be adapted to those given in `ranges`.

        Examples
        --------
        ::

          from PyAstronomy import funcFit as fuf
          import numpy as np
          import matplotlib.pylab as plt

          x = np.linspace(0,30,1000)
          gauss = fuf.GaussFit1d()
          gauss["A"] = 1
          gauss["mu"] = 23.
          gauss["sig"] = 0.5
          yerr = np.random.normal(0., 0.05, len(x))
          y = gauss.evaluate(x) + yerr
          ranges = {"A":[0,5],"mu":10, "sig":[0.1,1.0]}
          gauss.autoFitMCMC(x, y, ranges, yerr=yerr, iter=1000, picky=False)

          plt.plot(x, y, 'k+')
          plt.plot(x, gauss.evaluate(x), 'r--')
          plt.show()
        """
        X0, lims, steps = self.MCMCautoParameters(
            ranges, picky=picky, stepsize=stepsize)
        self.fitMCMC(x, y, X0, lims, steps, yerr=yerr, pymcPars=pymcPars, pyy=pyy,
                     potentials=potentials, dbfile=dbfile, dbArgs=dbArgs, **sampleArgs)

    def fitMCMC(self, x, y, X0, Lims, Steps, yerr=None, pymcPars={}, pyy=None,
                potentials=[], dbfile="mcmcSample.tmp", dbArgs={}, **sampleArgs):
        raise(PE.PyANotImplemented("This method has to be implemented"))


class _OndeDFitParBase:

    def __init__(self, parList, **kwargs):
        self.pars = Params(parList)
        # Set up the propMap (used to create a mapping between "property" (e.g., mu)
        # and variable name (e.g., "mu_1"))
        self.propMap = {}
        for p in parList:
            self.propMap[p] = p
        # Enable Root naming and component counters
        self.naming = ModelNameIdentBase(**kwargs)

    def __getitem__(self, specifier, **kwargs):
        name = self.naming.convertSpecifier(specifier)[0]
        if not name in self.propMap:
            if name in self.parameters():
                return self.pars.__getitem__(name)
            else:
                raise(PE.PyAValError("No such parameter: " + str(name), where="OneDFit::__getitem__",
                                     why="Did you perhaps misspell the parameter name?"))
        return self.pars.__getitem__(self.propMap[name])

    def __setitem__(self, specifier, value):
        name = self.naming.convertSpecifier(specifier)[0]
        if not name in self.propMap:
            if name in self.parameters():
                self.pars.assignValue({name: value})
                return
            raise(PE.PyAValError("No such parameter: " + str(name), where="OneDFit::__setitem__",
                                 why="Did you perhaps misspell the parameter name?"))
        self.pars.assignValue({self.propMap[name]: value})

    def hasVariable(self, specifier):
        """
        Determine whether the variable exists.

        Parameters
        ----------
        specifier : string or tuple
            Defines the name of the variable to be checked by
            string or specifier.

        Returns
        -------
        Flag : boolean
            True if the variable exists.
        """
        varName = self.naming.specifierToName(specifier)
        return self.pars.hasParam(varName)

    def assignValue(self, specval):
        """
        Assign new values to variables.

        Parameters
        ----------
        specval : dictionary
            Dictionary mapping variable names (wither given as string
            or specifier tuple) to value.
        """
        return self.pars.assignValue(self.naming.specifierToName(specval))

    def assignValues(self, specval):
        return self.assignValue(specval)
    assignValues.__doc__ = assignValue.__doc__

    def thaw(self, specifiers):
        """
        Consider variables fixed.

        Parameters
        ----------
        specifiers : list of strings or tuples
            The names of the variables to be fixed. Either given as
            string or specifier tuple.
        """
        self.pars.thaw(self.naming.specifierToName(specifiers))

    def freeze(self, specifiers):
        """
        Consider variables free to float.

        Parameters
        ----------
        specifiers : list of strings or tuples
            The names of the variables to be thawed. Either given as
            string or specifier tuple.
        """
        self.pars.freeze(self.naming.specifierToName(specifiers))

    def setRestriction(self, restricts):
        """
        Define restrictions.

        Parameters
        ----------
        restricts : dictionary
            A dictionary associating a variable (given as string or
            specifier tuple) with a restriction of the form:
            [lower, upper]. Use "None" where no restrictions shall
            be applied.
        """
        self.pars.setRestriction(self.naming.specifierToName(restricts))

    def untie(self, parName, forceFree=False):
        """
        Remove all relations of parameter parName, i.e., the parameter is not dependend 
        on other parameters. The parameter parName is set to "freeze".

        Parameters
        ----------
        parName : string
            The name of the dependent variable which should become "unrelated".
        forceFree : boolean
            Set parName to "free" instead of "frozen" if set to True.
        """
        self.pars.untie(self.naming.specifierToName(
            parName), forceFree=forceFree)

    def relate(self, dependentVar, independentVars, func=None, **kwargs):
        """
        Define a relation.

        Parameters
        ----------
        dependentVar : string or tuple
            The dependent variable given by string or specifier tuple.
        independentVars : string or list of strings
            The independent variables. You may also use specifier tuples
            to address the variables.
        func : callable
            The functional form of the relation. Must take the
            independent variables as arguments and returns the
            value of the dependent variable. If None is given,
            equality will be assumed.
        """
        self.pars.relate(self.naming.specifierToName(
            dependentVar), self.naming.specifierToName(independentVars), func, **kwargs)

    def getRelationsOf(self, specifier):
        """
        Return relations of a variable.

        Parameters
        ----------
        specifier : string or tuple
            Variable name or specifier tuple.

        Returns
        -------
        Relations : list of relations
            Those relations in which `specifier` is the dependent variable, i.e.,
            specifier = f(x,y).
        """
        return self.pars.getRelationsOf(self.naming.specifierToName(specifier))

    def freeParameters(self):
        return self.pars.freeParameters()
    freeParameters.__doc__ = Params.freeParameters.__doc__

    def frozenParameters(self):
        return self.pars.frozenParameters()
    frozenParameters.__doc__ = Params.frozenParameters.__doc__

    def freeParamNames(self):
        return self.pars.freeParamNames()
    freeParamNames.__doc__ = Params.freeParamNames.__doc__

    def numberOfFreeParams(self):
        return self.pars.numberOfFreeParams()
    numberOfFreeParams.__doc__ = Params.numberOfFreeParams.__doc__

    def getRestrictions(self):
        return self.pars.getRestrictions()
    getRestrictions.__doc__ = Params.getRestrictions.__doc__

    def delRestriction(self, parName):
        return self.pars.delRestriction(parName)
    delRestriction.__doc__ = Params.delRestriction.__doc__

    def availableParameters(self):
        return self.pars.availableParameters()
    availableParameters.__doc__ = Params.availableParameters.__doc__

    def parameters(self):
        return self.pars.parameters()
    parameters.__doc__ = Params.parameters.__doc__

    def saveState(self, *args, **kwargs):
        return self.pars.saveState(*args, **kwargs)
    saveState.__doc__ = Params.saveState.__doc__

    def restoreState(self, resource):
        self.pars.restoreState(resource)
    restoreState.__doc__ = Params.restoreState.__doc__

    def addConditionalRestriction(self, *args):
        return self.pars.addConditionalRestriction(*args)
    addConditionalRestriction.__doc__ = Params.addConditionalRestriction.__doc__

    def removeConditionalRestriction(self, *args):
        return self.pars.removeConditionalRestriction(*args)
    removeConditionalRestriction.__doc__ = Params.removeConditionalRestriction.__doc__

    def showConditionalRestrictions(self, **kwargs):
        return self.pars.showConditionalRestrictions(**kwargs)
    showConditionalRestrictions.__doc__ = Params.showConditionalRestrictions.__doc__


class IFitterBase:
    """
    Base class for internal fitter.
    """

    def __init__(self):
        self._objfval = None

    def fit(self, minifunc, x0):
        """
        Carry out the minimization.

        Parameters
        ----------
        minifunc : callable
            Objective function.
        x0 : list
            Starting values

        Returns
        -------
        Result : list or tuple
            First item is a list of the best-fit values and second
            item is the value of the objective function.
        """
        raise(PE.PyANotImplemented(
            "IFitterBase: fit method needs to be implemented"))

    def _digestkwargs(self, kwargs):
        """
        Check validity of keywords.

        Tests whether all given keywords are included in the list
        of "allowed" keywords (attribute _allowedKWs).

        Parameters
        ----------
        kwargs : dictionary
            Keyword arguments handed to the minimizer.
        """
        if not hasattr(self, "_allowedKWs"):
            raise(PE.PyANotImplemented("No _allowedKWs attribute. Must be set to use '_digestkwargs'.",
                                       solution="Specify _allowedKWs, e.g., in constructor."))
        for k in six.iterkeys(kwargs):
            if not k in self._allowedKWs:
                raise(PE.PyAValError("Keyword " + k + " not allowed in algorithm: " + self.name,
                                     solution="Allowed keywords: " +
                                     ', '.join(self._allowedKWs),
                                     where=self.name))

    def getObjFuncValue(self):
        """
        Access value of objective function.

        Returns
        -------
        value : float
            Value of objective function.
        """
        return self._objfval

    def __call__(self, *args, **kwargs):
        """
        Wrapper around the actual fit method.
        """
        return self.fit(*args, **kwargs)


class ScipyFMIN(IFitterBase):
    """
    Wrapper around scipy.optimize.fmin.
    """

    def __init__(self, *args, **kwargs):
        IFitterBase.__init__(self)
        self._allowedKWs = ["xtol", "ftol", "maxiter",
                            "maxfun", "disp", "retall", "callback"]
        self.name = "scipy.optimize.fmin"

    def fit(self, miniFunc, x0, *fminpars, **fminargs):
        """
        Wrapper around scipy.optimize.fmin.
        """
        self._digestkwargs(fminargs)
        self._result = sco.fmin(miniFunc, x0, *fminpars,
                                full_output=True, **fminargs)
        self._objfval = self._result[1]
        return self._result[0], self._result[1]


class FuFNM(IFitterBase):

    def __init__(self, *args, **kwargs):
        IFitterBase.__init__(self)
        self._allowedKWs = ["initDelta", "maxIter", "callback", "nmCritLim"]
        self.name = "funcFit NM65"
        self._obj = args[0]
        self._nm = NelderMead(**kwargs)

    def fit(self, miniFunc, x0, *fminpars, **fminargs):
        """
        Wrapper around funcFit's NelderMead implementation.

        See the implementation of the `fit` method of :py:func:`NelderMead`
        for the available keyword arguments (fminargs).
        """
        self._digestkwargs(fminargs)
        self._bestFit = self._nm.fit(
            self._obj, self._obj._fufDS, self._obj.miniFunc, **fminargs)
        self._bestFitVals = self._obj.pars.getFreeParams()
        self._objfval = self._obj.miniFunc(self._bestFitVals)
        return self._bestFitVals, self._objfval


FuFNM.__doc__ = NelderMead.__doc__


class FuFPrior(object):
    """
    A number of priors.

    Properly initialized, an object of type "FuFPrior"
    is callable. On call, it expects a dictionary with
    the names and values of the parameters as first
    argument and the name of the parameter under consideration
    as second argument. The return value should be the
    (natural) logarithm of the associated prior probability
    distribution.

    Parameters
    ----------
    lnp : string, {uniform, jeffreyPS, gaussian, limuniform}
        uniform : improper uniform prior. limuniform: proper
        uniform prior. 'lower' and 'upper' define the lower and
        upper bounds of the interval.
        jeffreyPS: Jeffreys
        prior for a Poisson scaling parameter. gaussian: A
        Gaussian prior. The keyswords 'mu' and 'sig' must be
        specified to define the mean and standard deviation of
        the Gaussian.
    """

    def _uniform(self, **kwargs):
        def uniform(ps, n, **rest):
            return 0.0
        return uniform

    def _uniformLimit(self, **kwargs):
        if kwargs["upper"] < kwargs["lower"]:
            raise(PE.PyAValError("upper needs to be larger than lower",
                                 where="FuFPrior (limited uniform distribution)",
                                 solution="Adapt upper and lower."))
        p = np.log(1.0 / (kwargs["upper"] - kwargs["lower"]))

        def unilimit(ps, n, **rest):
            if (ps[n] >= kwargs["lower"]) and (ps[n] <= kwargs["upper"]):
                return p
            else:
                return -np.Inf
        return unilimit

    def _jeffreyPoissonScale(self, **kwargs):
        def jps(ps, n, **rest):
            return -0.5 * np.log(ps[n])
        return jps

    def _gaussian(self, **kwargs):
        r = -0.5 * np.log(2.0 * np.pi * kwargs["sig"]**2)

        def gaussianPrior(ps, n, **rest):
            return r - (ps[n] - kwargs["mu"])**2 / (2.0 * kwargs["sig"]**2)
        return gaussianPrior

    def _callDelegator(self, *args, **kwargs):
        """ Overwritten by the method to represent __call__ """
        raise(PE.PyANotImplemented("_callDelegator is not implemented."))

    def __call__(self, *args, **kwargs):
        return self._callDelegator(*args, **kwargs)

    def __init__(self, lnp, **kwargs):
        if isinstance(lnp, six.string_types):
            if lnp == "uniform":
                self._callDelegator = self._uniform(**kwargs)
            elif lnp == "limuniform":
                self._callDelegator = self._uniformLimit(**kwargs)
            elif lnp == "jeffreyPS":
                self._callDelegator = self._jeffreyPoissonScale(**kwargs)
            elif lnp == "gaussian":
                self._callDelegator = self._gaussian(**kwargs)
            else:
                raise(PE.PyAValError("No prior defined for " + str(lnp),
                                     where="FuFPrior",
                                     solution="Use either of {uniform, limuniform, jeffreyPS, gaussian}"))


class OneDFit(_OndeDFitParBase, _PyMCSampler):
    """
    The base class for fitting objects.

    Parameters
    ----------
    parList : list of strings
        Contains the names of the properties
        defining the model. By default, variables of the same name
        are used to represent them.

    Attributes
    ----------
    model : array
        Used by the `updateModel` method to store the
        evaluated model for current parameter settings.
        Holds the best-fit model after a call to a fit method.
    penaltyFactor : float
        The penalty factor used to apply penalties for
        enforcing restrictions (default = 10**20).

    Notes
    -----

    The purpose of the class

      The purpose of this class is to provide a convenient interface
      to various fitting algorithms.
      It provides the functionality, which allows for parameter
      fitting, but does not implement a particular model.
      The class can be used to fit any kind of model, which has
      to be implemented in a class, which inherits from
      the *OneDFit* class.

    Management of fitting parameters

      The fitting parameters are managed by a *Params* class
      instance, which provides a wealth of possibilities to
      influence the behavior of the parameters during the
      fitting process. This includes deciding whether a particular
      parameter is supposed to be a free fitting parameter,
      applying restrictions to limit the valid range for a
      parameters, or the introduction of functional dependencies among
      different parameters.

    Properties versus variable names

      Each model is described by a number of *properties*, such
      as, for example, mass and radius. These may be represented
      by arbitrarily named variables. Normally, it is convenient to
      name the variables according to the properties they describe,
      which is the default behavior. However, in some cases, for example
      if a model consists of two equal subcomponents, such a naming
      scheme leads to nonunique variable names, which has to be avoided.
      Now it is necessary to distinguish between the *property* and
      the describing variable. This class uses the `propMap`
      dictionary, which maps property name to
      variable name to manage these situations.

    Combining fitting objects

      Often, it can be convenient to combine a number of simple models
      to form a new, more complex one. The *OneDFit* class allows to
      combine objects using the arithmetic operators +-\*/,
      and the power (\*\*) operator.

    Naming scheme for models

      For simple models it is convenient to use a one-to-one mapping
      between property and variable name. It may, however, become
      necessary to deviate from this scheme, for example, to keep
      variable names unique. This class supports the following naming scheme:
      Each model has a "root name", which is supposed to be a concise
      string describing the model (for instance, "Gaussian"). The root
      name is attached to the property name using an underscore.
      If a complex model consists of more than one component
      with the same root name, a component counter, enclosed in
      parenthesis, is attached to the variable name.
      A variable name could, for example, look like: "mu_Gaussian(1)".

    Methods to be implemented in a model class

      A valid model class inheriting this interface class must provide
      the following methods

     - **__init__()** - The constructor.

                        Defines the set of properties describing the model.
     - **evaluate(x)** - An *evaluate* method.

                         This method takes a single
                         argument, x, which is an array of points at which
                         the model is to be evaluated. To access the
                         current model parameters, this method should use
                         the set/getitem methods. The return value is
                         an array holding the model evaluated at the
                         points given by `x`.
    """

    def __init__(self, parList, **kwargs):
        _OndeDFitParBase.__init__(self, parList, **kwargs)
        # Left and right compo(nent) are necessary for combining models
        self.leftCompo = None
        self.rightCompo = None
        self.penaltyFactor = 1e20
        self.model = None
        self._fufDS = None
        self.fitResult = None
        # Determines whether steppar can be used
        self._stepparEnabled = False

    def _compoWalk(self):
        """
        TBD
        """
        def walk(c, refs):
            refs.append(c)
            if c.leftCompo is not None:
                walk(c.leftCompo, refs)
            if c.rightCompo is not None:
                walk(c.rightCompo, refs)

        refs = []
        walk(self, refs)
        for c in refs:
            yield c

    def renameVariable(self, oldName, newName):
        """
        Change name of variable.

        Parameters
        ----------
        oldName : string
            Current variable name.
        newName : string
            New variable name.

        Notes
        -----
        Variable names and properties are not the same.
        """
        # First, walk down the left and right components (for combined models)
        # and change the variable names.
        if self.leftCompo is not None:
            try:
                self.leftCompo.renameVariable(oldName, newName)
            except PE.PyAValError:
                pass
        if self.rightCompo is not None:
            try:
                self.rightCompo.renameVariable(oldName, newName)
            except PE.PyAValError:
                pass
        # Now do the same for the "top" component
        if newName == oldName:
            # Ignore identical transformations
            return
        if newName in list(self.propMap.values()):
            raise(PE.PyANameClash("A variable named " + newName +
                                  " does already exist.", where="OneDFit::renameVariable"))
        if newName in self.propMap:
            if self.propMap[newName] != oldName:
                raise(PE.PyANameClash("You may not assign a name to a variable, which corresponds to the name of another property.",
                                      where="OneDFit::renameVariable"))
        if not oldName in list(self.propMap.values()):
            raise(PE.PyAValError("A variable named " + oldName +
                                 " does not exist.", where="OneDFit::renameVariable"))
        for k in six.iterkeys(self.propMap):
            if self.propMap[k] == oldName:
                self.propMap[k] = newName
                break
        # Tell the parameter class about the renaming (if this did not already happen)
        if not newName in six.iterkeys(self.pars.parameters()):
            self.pars.renameParameter(oldName, newName)

    def _isComposed(self):
        """
        Determines whether current model is "composed".

        A model is composed, if there are left and right components.

        Returns True if model is composed and False otherwise.
        """
        return ((self.leftCompo is not None) and (self.rightCompo is not None))

    def __combineRemapping(self, left, right):
        """
        This member is essentially a renaming machine. When combining models
        it can easily happen that two variables share the same name. If the
        models are combined, unique variable names are needed. This method
        uses the "root name" and "component counter" to assign new, unique
        names to the variables. 

        Parameters:
          - `left`, `right` - Two fitting objects (derived from OneDFit).
        """

        def extendCoDat(coDat, c):
            ident = c.naming.getRoot()
            if not ident in coDat:
                coDat[ident] = [c]
            else:
                coDat[ident].append(c)
            return coDat

        # Build up a dictionary assigning root to a list of corresponding components
        coDat = {}
        for c in left._compoWalk():
            extendCoDat(coDat, c)
        for c in right._compoWalk():
            extendCoDat(coDat, c)

        for k in six.iterkeys(coDat):
            # Loop over all available root names
            if len(coDat[k]) == 1:
                # Only a single component with that root, no problem, no renaming
                continue
            elif len(coDat[k]) >= 2:
                # Get all component counters pertaining to that root
                cs = np.array([c.naming.getComponentCounter()
                               for c in coDat[k]])
                # What is the current maximum
                highest = np.max(cs)
                # Check whether one or more counters are zero;
                # if so, increase the number and set it to "highest + 1".
                for c in coDat[k]:
                    if c.naming.getComponentCounter() == 0:
                        c.naming.setComponentCounter(highest + 1)
                        highest += 1
                # Now re-check whether counters are unique
                while True:
                    cs = np.array([c.naming.getComponentCounter()
                                   for c in coDat[k]])
                    u = np.unique(cs)
                    if len(u) == len(cs):
                        # All counters are unique, work is finished
                        break
                    # Set lowest counter
                    again = False
                    for c in coDat[k]:
                        # Loop over components
                        if len(np.where(cs == c.naming.getComponentCounter())[0]) > 1:
                            # If there are nonunique counters, make them unique and retry
                            c.naming.setComponentCounter(highest + 1)
                            highest += 1
                            again = True
                            break
                    if not again:
                        break

        def renameIfNecessary(c, baseComponent):
            # Rename all variables according to standard naming rules
            for prop in six.iterkeys(c.propMap):
                # Compose new name using "standard conventions"
                newName = c.naming.composeVariableName(prop)
                if not (newName == c.propMap[prop]):
                    # Name has changed
                    baseComponent.renameVariable(c.propMap[prop], newName)

        # Rename all variables mentioned in NONCOMPOSED models
        for c in left._compoWalk():
            if c._isComposed():
                continue
            renameIfNecessary(c, left)
        for c in right._compoWalk():
            if c._isComposed():
                continue
            renameIfNecessary(c, right)

    def _newPars(self, npars):
        """
        Change pars reference for left/right Compo recursively.
        """
        if self.leftCompo is not None:
            self.leftCompo.pars = npars
            self.leftCompo._newPars(npars)
        if self.rightCompo is not None:
            self.rightCompo.pars = npars
            self.rightCompo._newPars(npars)

    def __combineFittingObjects(self, right):
        """
        Creates and returns a fitting object combining the properties/variables of
        the current objects and that given by `right`.

        Parameters:
          - `right` - A fitting object (derived from OneDFit).
        """
        # Obtain deep copies of the involved fitting objects
        # The actual reason for doing this, is that the evaluate
        # methods have to be conserved (copied).
        left = copy.deepcopy(self)
        right = copy.deepcopy(right)
        # Find new parameter names to avoid parameter clashes
        # (parameter names must be unique)
        self.__combineRemapping(left, right)
        # Combine the Param instances
        npars = left.pars + right.pars
        # Obtain a new Fitting object
        result = OneDFit(list(npars.parameters().keys()))
        # Assign new combined Param instance
        result.pars = npars
        # Assign new parameter class to ``old'' individual fitting
        # objects.
        left.pars = result.pars
        right.pars = result.pars

        # Assign left and right component (for evaluation)
        result.leftCompo = left
        result.rightCompo = right
        # Recursively replace pars property by new pars reference
        result._newPars(npars)
        result.setRootName("combinedModel")
        return result

    def __add__(self, right):
        result = self.__combineFittingObjects(right)
        # For a discussion about bound instance methods see
        # https://stackoverflow.com/questions/972/adding-a-method-to-an-existing-object-instance
        result.evaluate = types.MethodType(addEval, result)
        # Save the 'operator' relating left and right component
        result._operator = '+'
        return result

    def __sub__(self, right):
        result = self.__combineFittingObjects(right)
        result.evaluate = types.MethodType(subEval, result)
        # Save the 'operator' relating left and right component
        result._operator = '-'
        return result

    def __mul__(self, right):
        result = self.__combineFittingObjects(right)
        result.evaluate = types.MethodType(mulEval, result)
        # Save the 'operator' relating left and right component
        result._operator = '*'
        return result

    def __div__(self, right):
        result = self.__combineFittingObjects(right)
        result.evaluate = types.MethodType(divEval, result)
        # Save the 'operator' relating left and right component
        result._operator = '/'
        return result

    def __pow__(self, right):
        result = self.__combineFittingObjects(right)
        result.evaluate = types.MethodType(powEval, result)
        # Save the 'operator' relating left and right component
        result._operator = '**'
        return result

    def __sqrDiff(self):
        @MiniFunc(self)
        def miniSqrDiff(odf, P):
            # Calculate squared difference
            chi = np.sum((self._fufDS.y - self.model)**2)
            return chi
        return miniSqrDiff

    def __chiSqr(self):
        @MiniFunc(self)
        def miniChiSqr(odf, P):
            # Calculate chi^2 and apply penalty if boundaries are violated.
            chi = np.sum(((self._fufDS.y - self.model) / self._fufDS.yerr)**2)
            return chi
        return miniChiSqr

    def __cash79(self):
        @MiniFunc(self)
        def miniCash79(odf, P):
            # Calculate Cash statistics according to Cash 1979 (ApJ 228, 939)
            c = -2.0 * np.sum(self._fufDS.y * np.log(self.model) - self.model)
            return c
        return miniCash79

    def setRootName(self, root, rename=False):
        """
        Define the root name of the model.

        Parameters
        ----------
        root : string
            A concise description of the model.
        rename : bool, optional, default=False
            If true, all model variables will be renaming using the root.
        """
        if not hasattr(self, "naming"):
            raise(PE.PyAUnclassifiedError("No attribute 'naming' found.",
                                          solution="Did you call the constructor of OneDFit before you tried to access this? If not, you should do that."))
        self.naming.setRootName(root)
        if rename:
            for par in six.iterkeys(self.propMap):
                oldName = self.propMap[par]
                newName = self.naming.composeVariableName(par)
                self.renameVariable(oldName, newName)

    def updateModel(self):
        """
        Recalculate the model using current settings.

        Notes
        -----
        Updates the `model` attribute of the class by
        calling `evaluate` using the `x` attribute, which
        is, e.g., assigned on call to a fit method.
        """
        if self._fufDS is None:
            return
        self.model = self.evaluate(self._fufDS.x)

    def setPenaltyFactor(self, penalFac):
        """
        Change the penalty factor.

        Parameters
        ----------
        penalFac : float
            The penalty factor (default is 1e20).

        Notes
        -----
        May also be done by accessing the `penalty` property
        directly.
        """
        self.penaltyFactor = penalFac

    def description(self, parenthesis=False):
        """
        Returns a description of the model based on the names of the individual components.

        Parameters
        ----------
        parenthesis : boolean, optional
            If True, the entire expression/description will be enclosed in
            parenthesis. The default is False.

        Returns
        -------
        Description : string
            Description of the model.
        """
        c = self
        if c._isComposed():
            if parenthesis:
                return "(" + c.leftCompo.description(True) + " " + c._operator + " " + c.rightCompo.description(True) + ")"
            else:
                return c.leftCompo.description(True) + " " + c._operator + " " + c.rightCompo.description(True)
        else:
            if c.naming.getComponentCounter() > 0:
                return c.naming.getRoot() + "(No. " + str(c.naming.getComponentCounter()) + ")"
            else:
                return c.naming.getRoot()

    def parameterSummary(self, toScreen=True, prefix="", sorting='none'):
        """
        Writes a summary of the parameters in text form.

        Parameters
        ----------
        toScreen : bool, optional, default = True
            If True, the output is written to the screen.
        prefix : string, optional, default = ""
            A prefix for each line (e.g., '#').
        sorting : string, optional, {'none', 'ps'}
            Determines the order in which the parameters are
            printed out. If 'none' is given (default), no
            particular order is imposed. If 'ps', Python's
            `sorting` routine is used to impose an order.

        Returns
        -------
        A list of strings containing the text.
        """
        lines = []

        # Collect roots and associated counters and
        # (root, counter) and associated component
        rootNums = {}
        rootNoComp = {}
        for c in self._compoWalk():
            if c._isComposed():
                continue
            # It is a noncomposed component
            if not c.naming.getRoot() in rootNums:
                rootNums[c.naming.getRoot()] = [c.naming.getComponentCounter()]
            else:
                rootNums[c.naming.getRoot()].append(
                    c.naming.getComponentCounter())
            rootNoComp[(c.naming.getRoot(), c.naming.getComponentCounter())] = c

        # Find maximal length of property, root, identifier (= root (counter)),
        # and variable. This is done to provide well formatted output
        maxPropLen = 0
        maxRootLen = 0
        maxIdentLen = 0
        maxVarLen = 0
        for k, c in six.iteritems(rootNoComp):
            maxRootLen = max(maxRootLen, len(k[0]))
            maxIdentLen = max(maxIdentLen, len(c.naming.identifier()))
            for prop, var in six.iteritems(c.propMap):
                maxPropLen = max(maxPropLen, len(prop))
                maxVarLen = max(maxVarLen, len(var))

        # Walk through all root and counter to provide output
        for root in sorted(rootNums.keys()):
            for counter in sorted(rootNums[root]):
                c = rootNoComp[(root, counter)]
                cono = "Parameters for Component: "
                if root == "":
                    cono += "unnamed"
                else:
                    cono += root
                if counter != 0:
                    cono += " (No. " + str(counter) + ")"
                lines.append("-" * len(cono))
                lines.append(cono)
                lines.append("-" * len(cono))

                # Sorting of the parameter names (component-wise)
                cprops = list(c.propMap.keys())
                if sorting == "none":
                    pass
                elif sorting == "ps":
                    cprops = sorted(cprops)
                else:
                    raise(PE.PyAValError("Unknown sorting mode: " + str(sorting),
                                         where="OneDFit::parameterSummary",
                                         solution="Use 'None' or 'ps'"))

                for prop in cprops:
                    var = c.propMap[prop]
                    # Loop over properties and variables
                    restricted = (self.getRestrictions()[var][0] is not None) or (
                        self.getRestrictions()[var][1] is not None)
                    lines.append(("Parameter: %" + str(maxPropLen) + "s  %" + str(maxIdentLen) + "s, [%" + str(maxVarLen) + "s], value: % 12g, free: %5s, restricted: %5s, related: %5s") %
                                 (prop, c.naming.identifier(), var, self[var], var in self.freeParameters(), restricted, (len(self.getRelationsOf(var)) != 0)))
                    app = " " * 4
                    if restricted:
                        res = self.getRestrictions()[var]
                        app += "Restriction: ["
                        if res[0] is not None:
                            app += "% g, " % res[0]
                        else:
                            app += "None, "
                        if res[1] is not None:
                            app += "% g]" % res[1]
                        else:
                            app += " None]"
                    relas = self.getRelationsOf(var)
                    if len(relas) != 0:
                        if restricted:
                            app += ", "
                        for r in relas:
                            app += " Relation: " + r[0] + " = f("
                            for k in range(len(r[2]) - 1):
                                app += r[2][k] + ", "
                            app += r[2][len(r[2]) - 1] + ")"

                    if re.match("\s*[^\s]+", app):
                        lines.append(app)

        for i in smo.range(len(lines)):
            lines[i] = prefix + lines[i]

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
        else:
            if not hasattr(miniFunc, '__call__'):
                raise(PE.PyAValError("`miniFunc` is neither None, a valid string, or a function.",
                                     where="OneDFit::fit",
                                     solution="Use, e.g., 'chisqr' or another valid choice from the documentation."))

            # A function has been specified
            self.miniFunc = miniFunc
            return

    def fitMCMC(self, x, y, X0, Lims, Steps, yerr=None, pymcPars=None, pyy=None,
                potentials=None, dbfile="mcmcSample.tmp", quiet=False, dbArgs=None,
                adaptiveMetropolis=False, **sampleArgs):
        """
        Carry out MCMC fit/error estimation.

        This member is designed to provide a flexible but *easy to use*
        interface to the capabilities of pymc. In the simplest case,
        it assumes a Poisson or Gaussian distribution of data points
        and uses continuous, uniform variables (all free fitting
        variables) with starting values defined by `X0`, Limits
        given by `Lims`, and step sizes given by `Steps`
        to sample from the posterior.

        .. note::
           The result (the Marchov-Chain/PyMC MCMC-object) will
           be contained in the self.MCMC property; the output dictionary
           of MCMC.stats()
           (Mean, HPD interval etc.) is saved to self.basicStats.

        Parameters
        ----------
        x, y : array 
            Array providing the abscissa and ordinate values of the
            data points.
            For the MCMC process, it is essential to know
            the underlying distribution of the data points.
            `fitMCMC` assumes Poisson distributed data if
            `yerr` is not specified and Gaussian data if it
            is specified. If other distributions shall be used,
            the `pyy` parameter must contain a
            *PyMC* random variable specifying it.
        yerr : array, optional
            The error of the data points. If not specified, Poisson
            distributed data will be assumed.
        pyy : PyMC stochastic, optional,
            A PyMC variable containing the \
            data. This can be useful if the distribution is
            neither Poisson nor
            Gaussian as otherwise assumed by this method.
        X0 : dict
             A dictionary holding {"parName":value, ...} specifying
             the start values. Note that parameters treated in
             `pymcPars` must not be part of this.
        Lims : dict,
            A dictionary of type {"ParName:[min,max], ...}
            specifying the lower and upper limit of a variable.
            Note that parameters treated in pymcPars must not
            be part of this.
        Steps : dict
            A dictionary providing the step sizes for the MCMC sampler.
        pymcPars : dict, optional,
            Holds a dictionary of the form {"parName":PyMC-Variable, ...}.
            `pymcPars` can be used to specify a nonuniform
            distribution for a parameter.
        potentials : list of PyMC potentials, optional,
            Can be used to provide a list of PyMC potentials,
            which may be needed to provide priors.
        dbfile : string, optional,
            The name of the output file, which is to hold the MCMC chain
            produced during sampling (default="mcmcSample.tmp").
        quiet : bool, optional,
            Set to True in order to suppress the text output.
        sampleArgs : dict, optional,
            Here additional keywords can be specified, which
            will be handed to the `isample` member of PyMC. Most notably
            it is useful to specify **iter**, **burn**, and **thin**.
            For other possibilities see PyMC documentation.
        dbArgs : dict, optional,
            Keyword arguments specifying the trace database to be used.
            By default a pickle database named according to the
            `dbfile` keyword.
            You may also use an hdf5 database.
        adaptiveMetropolis : bool, optional,
            Set to true in order to use the AdaptiveMetropolis step method. 
        """
        global _pymcImport
        if not _pymcImport:
            raise(PE.PyARequiredImport("pymc package could not be imported.",
                                       solution="Install pymc (see http://code.google.com/p/pymc/"))
        if (pymc.__version__).startswith('3'):
            raise(PE.PyARequiredImport("pymc is available in version '" + str(pymc.__version__) + "'. " +
                                       "Please note that funcFit does not yet support pymc3, which is a " +
                                       "complete rewrite of pymc2.",
                                       solution="Please install pymc in version 2.3."))
        # Assign mutable default parameters
        if pymcPars is None:
            pymcPars = {}
        if dbArgs is None:
            dbArgs = {}
        if potentials is None:
            potentials = []
        # Assign attributes and check x, y, and yerr.
        self._fufDS = FufDS(x, y, yerr)
        # Copy the pymcPars and dbArgs dictionaries (prevents error on multiple sampler calls)
        pymcPars = pymcPars.copy()
        dbArgs = dbArgs.copy()
        # Get the names of the free parameters
        freeNames = self.freeParamNames()
        if not quiet:
            print("Free parameters: ", freeNames)
        # Check whether parameter lists are complete, define default steps
        # if necessary.
        self._dictComplete(freeNames, X0, "start values",
                           forget=list(pymcPars.keys()))
        self._dictComplete(freeNames, Lims, "limits",
                           forget=list(pymcPars.keys()))
        self._dictComplete(freeNames, Steps, "steps")

        # Define (or complete) the pymcPars dictionary by defining uniformly distributed
        # variables in the range [lim[0], lim[1]] with starting values defined by X0.
        for par in freeNames:
            if par in pymcPars:
                continue
            if not quiet:
                print("Using uniform distribution for parameter: ", par)
                print("  Start value: ",
                      X0[par], ", Limits = [", Lims[par][0], ", ", Lims[par][1], "]")
            pymcPars[par] = pymc.Uniform(
                par, lower=Lims[par][0], upper=Lims[par][1], value=X0[par], doc="Automatically assigned parameter.")

        # This function is used to update the model
        def getModel(**vals):
            self.assignValue(vals)
            self.updateModel()
            return self.model

        modelDet = pymc.Deterministic(
            eval=getModel,
            name='Model',
            parents=pymcPars,
            doc='The model',
            trace=True,
            verbose=0,
            dtype=float,
            plot=False,
            cache_depth=2)

        # Define the 'data' (y-values)
        if pyy is None:
            if yerr is None:
                if not quiet:
                    print(
                        "Assuming Poisson distribution for 'y'. Use 'pyy' parameter to change this!")
                pyy = pymc.Poisson(
                    "y", mu=modelDet, value=self._fufDS.y, observed=True)
            else:
                if not quiet:
                    print(
                        "Assuming Gaussian distribution for 'y'. Use 'pyy' parameter to change this!")
                pyy = pymc.Normal(
                    "y", mu=modelDet, tau=1.0 / self._fufDS.yerr**2, value=self._fufDS.y, observed=True)

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

        if not quiet:
            print("Using database arguments: ", dbArgs)
        self.MCMC = pymc.MCMC(Model, **dbArgs)

        # Tell the MCMC class to use the MH algorithm with specified step width
        if adaptiveMetropolis:
            self.MCMC.use_step_method(pymc.AdaptiveMetropolis, list(
                pymcPars.values()), shrink_if_necessary=True)
        else:
            for par in six.iterkeys(pymcPars):
                self.MCMC.use_step_method(
                    pymc.Metropolis, pymcPars[par], proposal_sd=Steps[par], proposal_distribution='Normal')

        if not "iter" in sampleArgs:
            sampleArgs["iter"] = 2000
        if not "burn" in sampleArgs:
            sampleArgs["burn"] = 0
        if not "thin" in sampleArgs:
            sampleArgs["thin"] = 1

        if not quiet:
            print("Giving the following arguments to 'isample':")
            print("  ", sampleArgs)

        if not quiet:
            self.MCMC.isample(**sampleArgs)
        else:
            self.MCMC.sample(**sampleArgs)
        self.basicStats = self.MCMC.stats()
        if not quiet:
            self._basicStatMCMCOutput(self.basicStats)

        # Setting values to ``best fit values'' (lowest deviance)
        mindex = np.argmin(self.MCMC.trace("deviance")[:])
        if not quiet:
            print("Searching for lowest-deviance solution.")
            print("  Index of lowest-deviance solution: ", mindex)
        for par in six.iterkeys(pymcPars):
            # Weird way of accessing element with index "mindex";
            # prevents slicing/indexing problem in pymc-Trace
            self[par] = (self.MCMC.trace(par)[mindex:mindex + 1])[0]
            if not quiet:
                print("  Parameter: ", par, ", value: ", self[par])
        self.updateModel()
        self.MCMC.db.close()

    def fitEMCEE(self, x=None, y=None, yerr=None, nwalker=None, priors=None, pots=None, scales=None,
                 sampleArgs=None, dbfile="chain.emcee", ps=None, emcp=None, toMD=True):
        """
        MCMC sampling using emcee package.

        Sample from the posterior probability distribution using the emcee
        package. By default the likelihood is calculated as -0.5 times the
        model chi-square value.

        The emcee sampler can be accessed via the `emceeSampler` attribute,
        which may be used to continue or manipulate sampling.

        Parameters
        ----------
        nwalker : int, optional
            The number of walker to be used. By default, two times the
            number of free parameters is used.
        scales : dictionary, optional
            The scales argument can be used to control the initial distribution
            of the walkers. By default, all walkers are distributed around the
            location given by the current state of the object, i.e., the current
            parameter values. In each direction, the walker are randomly distributed
            with a Gaussian distribution, whose default standard deviation is one.
            The scales argument can be used to control the width of Gaussians used
            to distribute the walkers.
        sampleArgs : dictionary, optional
            Number controlling the sampling process. Use 'burn' (int) to specify
            the number of burn-in iterations (default is 0). Via 'iters' (int)
            the numbers of iterations after the burn-in can be specified (default 1000).
            The 'process' (int) key can be used to control the number of iterations after
            which the progress bar is updated (default is iters/100). Note that the
            'progressbar' package must be installed to get a progress bar. Otherwise
            more mundane print statements will be used.  
        priors : dictionary, optional
            For each parameter, a primary can be specified. In particular, a
            prior is a callable, which is called with two arguments: first, a
            dictionary mapping the names of the free parameters to their
            current values, and second, a string specifying the name of the
            parameter for which the prior is to apply. The return value must be
            the logarithmic prior probability (natural logarithm). A number of default
            priors are available in the form of the `FuFPrior` class. By
            default, a uniform (improper) prior is used for all parameter, for
            which no other prior was specified.
        pots : list, optional
            A list of 'potentials'. A potential is a function, which is called using
            a dictionary holding the current value for all parameters and returns
            the logarithm of the associated probability. Potentials may, e.g., be
            used to implement certain relations between parameter values not otherwise
            accounted for. 
        dbfile : string, optional
            The result of the sampling, i.e., the chain(s), the corresponding
            values of the posterior, and the names of the free parameters are
            saved to the specified file (by default 'chain.emcee' is used).
            The traces stored there can be analyzed using the 'TraceAnalysis'
            class. Set this parameter to 'None' to avoid saving the results.
        ps : tuple, optional
            A tuple holding the current position and state of the sampler. This
            tuple is returned by this method. The `ps` argument can be used
            to continue sampling from the last state. Note that no burn-in will
            ne carried out and the other arguments should be given as previously
            to continue sampling successfully. 
        emcp : dictionary, optional
            Extra arguments handed to `EnsembleSampler` object.
        toMD : boolean, optional
            If True (default), the object is set to the lowest-deviance solution
            after sampling. Otherwise, it remains in a random state.
        """

        if not ic.check["emcee"]:
            raise(PE.PyARequiredImport("Could not import the 'emcee' package.",
                                       solution="Please install 'emcee'."))

        if (not x is None) and (not y is None) and (not yerr is None):
            # Assign attributes and check x, y, and yerr.
            self._fufDS = FufDS(x, y, yerr)
        elif (not x is None) and (not y is None) and (yerr is None):
            raise(PE.PyAValError("An error on the y values is required.",
                                 where="fitEMCEE",
                                 solution="Please specify 'yerr'"))
        if self._fufDS is None:
            raise(PE.PyAValError("Please specify the data completely.",
                                 where="fitEMCEE",
                                 solution="Specify x, y, and yerr."))

        if not self._fufDS.xyyerrDefined():
            raise(PE.PyAValError("Please specify the data completely. Either of x, y, and/or yerr are missing.",
                                 where="fitEMCEE",
                                 solution="Specify x, y, and yerr."))

        # Names and values of free parameters
        fps = self.freeParameters()
        # Names of the free parameters in specific order
        fpns = self.freeParamNames()
        # Number of dimensions
        ndims = len(fps)

        if ndims == 0:
            raise(PE.PyAValError("At least one free parameter is required for sampling.",
                                 where="fitEMCEE",
                                 solution="Use 'thaw' to free same parameters."))

        if not dbfile is None:
            if re.match(".*\.emcee$", dbfile) is None:
                PE.warn(PE.PyAValError("The db filename (" + str(dbfile) + ") does not end in .emcee. TraceAnalysis will not recognize it as an emcee trace file.",
                                       solution="Use a filename of the form *.emcee"))

        # Number of walkers
        if nwalker is None:
            self.nwalker = ndims * 2
        else:
            self.nwalker = nwalker

        if self.nwalker < ndims * 2:
            raise(PE.PyAValError("The number of walkers must be at least twice the number of free parameters.",
                                 where="fitEMCEE",
                                 solution="Increase the number of walkers."))
        if self.nwalker % 2 == 1:
            raise(PE.PyAValError("The number of walkers must be even.",
                                 where="fitEMCEE",
                                 solution="Use an even number of walkers."))

        # Use default prior for those parameters not listed
        if priors is None:
            priors = {}
        for n in fpns:
            if not n in priors:
                priors[n] = FuFPrior("uniform")

        # Ensure that potentials is at least an empty list
        if pots is None:
            pots = []

        # Chi square calculator
        chisqr = self.__chiSqr()

        def likeli(names, vals):
            # The likelihood function
            likeli = -0.5 * chisqr(vals)
            return likeli

        def lnpostdf(values):
            # Parameter-Value dictionary
            ps = dict(zip(fpns, values))
            # Check prior information
            prior_sum = 0
            for name in fpns:
                prior_sum += priors[name](ps, name)
            # If log prior is negative infinity, parameters
            # are out of range, so no need to evaluate the
            # likelihood function at this step:
            pdf = prior_sum
            if pdf == -np.inf:
                return pdf
            # Likelihood
            pdf += likeli(fpns, values)
            # Add information from potentials
            for p in pots:
                pdf += p(ps)
            if np.isnan(pdf):
                raise(PE.PyAValError("Posterior value is NaN for parameters: " + str(self.parameters()) + ".", \
                                     where="fitEmcee", \
                                     solution="Possibly, a prior (e.g., 'limuniform') can be used to restrict parameter range. " + \
                                     "Note that restrictions are not automatically converted into priors."))
            return pdf

        # Set default values for sampleArgs
        if sampleArgs is None:
            sampleArgs = {}
        if not "burn" in sampleArgs:
            sampleArgs["burn"] = 0
        if not "iters" in sampleArgs:
            sampleArgs["iters"] = 1000
        if not "progress" in sampleArgs:
            sampleArgs["progress"] = sampleArgs["iters"] / 100

        if ps is None:

            if emcp is None:
                emcp = {}

            # Generate the sampler
            self.emceeSampler = emcee.EnsembleSampler(
                self.nwalker, ndims, lnpostdf, **emcp)

            if scales is None:
                scales = {}

            # Generate starting values
            pos = []
            rrs = self.pars.getRestrictions()
            for _ in smo.range(self.nwalker):
                pos.append(np.zeros(ndims))
                for i, n in enumerate(fpns):
                    if not n in scales:
                        s = 1.0
                    else:
                        s = scales[n]
                    
                    # Trial counter -- avoid values beyond restrictions
                    tc = 0
                    while True:
                        if tc == 100:
                            raise(PE.PyAAlgorithmFailure("Could not determine valid starting point for parameter: " + str(fps) + " due to restrictions", \
                                                         where="fitEmcee", \
                                                         solution=["Try to use 'scale' to limit range of trial starting values.", \
                                                                   "Change starting value before MCMC call into valid range."]))
                        propval = np.random.normal(fps[n], s)
                        if n in rrs:
                            # There is a restriction
                            if (not rrs[n][0] is None) and (propval < rrs[n][0]):
                                tc += 1
                                continue
                            if (not rrs[n][1] is None) and (propval > rrs[n][1]):
                                tc += 1
                                continue
                        break
                    pos[-1][i] = propval

            # Default value for state
            state = None

            if sampleArgs["burn"] > 0:
                # Run burn-in
                pos, prob, state = self.emceeSampler.run_mcmc(
                    pos, sampleArgs["burn"])
                # Reset the chain to remove the burn-in samples.
                self.emceeSampler.reset()

        else:
            # Assign position and state from previous run
            pos, state = ps

        if (not sampleArgs["progress"] is None) and ic.check["progressbar"]:
            widgets = ['EMCEE progress: ', progressbar.Percentage(), ' ', progressbar.Bar(marker=progressbar.RotatingMarker()),
                       ' ', progressbar.ETA()]
            pbar = progressbar.ProgressBar(
                widgets=widgets, maxval=sampleArgs["iters"]).start()

        n = 0
        for pos, prob, state in self.emceeSampler.sample(pos, rstate0=state, iterations=sampleArgs["iters"], thin=1, storechain=True):
            n += 1
            if (not sampleArgs["progress"] is None) and (n % sampleArgs["progress"] == 0):
                if ic.check["progressbar"]:
                    pbar.update(n)
                else:
                    print("EMCEE: Reached iteration ",
                          n, " of ", sampleArgs["iters"])

        # Save the chain to a file
        if not dbfile is None:
            np.savez_compressed(open(dbfile, 'wb'), chain=self.emceeSampler.chain, lnp=self.emceeSampler.lnprobability,
                                pnames=np.array(fpns, dtype=np.unicode_))

        if toMD:
            # Set to lowest-deviance (highest likelihood) solution
            indimin = np.argmax(self.emceeSampler.lnprobability)
            for i, p in enumerate(self.freeParamNames()):
                self[p] = self.emceeSampler.flatchain[indimin, i]

        return pos, state

    def _resolveMinAlgo(self, minAlgo, default=None, mAA=None):
        """
        Resolve minimization algorithm (minAlgo).

        Parameters
        ----------
        minAlgo : callable, string, or None
            If None, the default will be used. If it is
            a callable, it will be assumed to be a valid
            implementation of a minimization algorithm. If
            it is a string, the function will try to resolve
            it.
        default : callable, string, or None
            The algorithm used, if minAlgo itself is None.
        mAA : dictionary, optional
            Keyword arguments given to constructor of minAlgo.

        Returns
        -------
        minAlgo : callable
            An instance of the minimization algorithm.
        """
        if minAlgo is None:
            # If not specified use default.
            minAlgo = default

        if minAlgo is None:
            # This should not happen. Default must be specified.
            raise(PE.PyAValError("No minimization algorithm specified.",
                                 solution="Use, e.g., minAlgo='spfmin' or minAlgo='fufnm'"))

        if hasattr(minAlgo, '__call__'):
            # It is a callable. Just assume that it can be used
            return minAlgo

        # Check keyword arguments to be handed to constructor
        # of the fitter
        if mAA is None:
            maa = {}
        else:
            maa = mAA

        if isinstance(minAlgo, six.string_types):
            if minAlgo == "fufnm":
                return FuFNM(self, **maa)
            elif minAlgo == "spfmin":
                # scipy.optimize.fmin
                global _scoImport
                if not _scoImport:
                    raise(PE.PyARequiredImport("SciPy.optimize could not be imported.",
                                               solution=["Install SciPy (see www.scipy.org/).",
                                                         "Use funcFit's Nelder-Mead simplex implementation (minAlgo='fufnm')"]))
                return ScipyFMIN(**maa)
            else:
                raise(PE.PyAValError("Unknown string identifier for minimization algorithm: " + minAlgo,
                                     solution="Use, e.g., minAlgo='spfmin' or minAlgo='fufnm'"))

        else:
            raise(PE.PyAValError("Could not resolve 'minAlgo'.",
                                 where="funcFit::OneDFit",
                                 solution="Use, e.g., minAlgo='spfmin' or minAlgo='fufnm'"))

    def fit(self, x, y, yerr=None, X0=None, minAlgo=None, mAA=None, miniFunc=None, printTime=False, *fminPars, **fminArgs):
        """
        Carries out a fit.

        Uses an internal optimizer to find the best-fit parameters.
        By default, the method uses the scipy.optimize.fmin algorithm
        (Nelder-Mead Simplex) to find the best-fit parameters. After the
        fit, the parameter values are set to the best-fit values.

        Parameters
        ----------
        x,y : arrays
            Specify the abscissa and ordinate values of the data points.
        yerr : array, optional
            Error of data values.
        X0 : list, optional
            The initial guess. If not provided, it will be assumed
            that the current parameter values already contains
            the initial guess.
        minAlgo : callable or string, optional
            The minimization algorithm. If not specified, scipy's 'fmin'
            implementation will be used. If a callable is given, it
            must adhere to funcFit's internal optimizer model.
            Valid strings are:
              - 'spfmin' : scipy.optimize.fmin
              - 'fufnm' : funcFit's implementation of the Nelder-Mead simplex algorithm.
        mAA : dictionary, optional
            Keyword arguments handed to the constructor of minAlgo, i.e.,
            minAlgo Arguments (mAA). Valid keywords depend on the choice
            of the fit algorithm.
        fminArgs : dict 
            Keywords passed to the minimization method
            (e.g., `xtol` or `ftol` for scipy.optimize.fmin).  
        fminPars :
            Non-keyword arguments passed to the  minimization method
            (e.g., fprime in scipy.optimize.fmin_ncg).
        miniFunc : None, string {"chisqr", "cash79"}, or function
            Function to be minimized. If None or "chisqr" is given,
            chi-square will be minimized. If "cash79" is given, the
            Cash statistics (Cash 1979, ApJ 228, 939, Eq. 5) will
            be minimized. If a function is specified, that, potentially
            user defined, function will be used to calculated the
            statistics, which will be minimized.
        printTime: boolean, optional
            If True, the number of seconds needed to carry out the fit
            is printed. Default is False. At any rate, the time in seconds
            is stored in the "requiredTime" attribute.
        """
        # Assign attributes and check x, y, and yerr.
        if (x is not None) and (y is not None):
            # This is a not-so-beautiful way to allow
            # calls from within without redefining the
            # data object. This way, the API can remain
            # unchanged.
            self._fufDS = FufDS(x, y, yerr)
        # Choose minimization algorithm
        self.mAA = copy.copy(mAA)
        self.minAlgo = self._resolveMinAlgo(minAlgo, default="spfmin", mAA=mAA)

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
        self.fminArgs, self.fminPars = fminArgs, fminPars
        # Carry out fit
        fitStartTime = timestamp()
        # For "historical" reasons, the fitResult must hold a list of the
        # best-fit parameter values as its first item and the final value of
        # the objective function as its second item.
        self.fitResult = self.minAlgo(self.miniFunc, self.pars.getFreeParams(), *self.fminPars,
                                      **self.fminArgs)
        # Set parameters and model to best-fit values
        self.pars.setFreeParams(self.fitResult[0])
        self.updateModel()
        self._stepparEnabled = True

        self.requiredTime = timestamp() - fitStartTime
        if printTime:
            print("The fit took " + str(self.requiredTime) + " seconds.")

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
            raise(PE.PyAOrderError("Before you can use steppar, you must call a function, which enables its use (e.g., `fit`).",
                                   solution="Call the `fit` method first and then try again."))
        if isinstance(pars, six.string_types):
            # Make it a list
            pars = [pars]
        # Check parameter consistency
        for p in pars:
            # Check existence
            tmp = self[p]
            if not p in ranges:
                raise(PE.PyAValError("There is no range for parameter: " + p,
                                     solution="Specify a range; e.g., {'xyz':[0.5,1.9,20,'lin']}"))
        # Function to extract function value from the fit result
        if extractFctVal is None:
            self._extractFctVal = self.__extractFunctionValue
        else:
            if not hasattr(extractFctVal, "__call__"):
                raise(PE.PyAValError("`extractFctVal` needs to be callable!",
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
                    raise(PE.PyAValError("If the range has 4 entries, the fourth must be a string specifying the mode.",
                                         solution="Use either 'lin' or 'log' as the fourth entry."))
                mode = r[3]
            if mode == 'lin':
                rs.append(np.linspace(r[0], r[1], r[2]))
            elif mode == 'log':
                # Calculate factor
                s = np.power((r[1] / r[0]), 1.0 / r[2])
                rs.append(r[0] * np.power(s, np.arange(r[2])))
            else:
                raise(PE.PyAValError("Unknown mode: " + str(mode),
                                     solution="Use either 'lin' or 'log'."))
        # Save state of object
        saveObj = self.saveState()
        saveFitResult = self.fitResult
        saveModel = self.model.copy()
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
            self.fit(None, None, yerr=None, minAlgo=self.minAlgo, miniFunc=self.miniFunc,
                     *self.fminPars, **self.fminArgs)
            # Build up result
            ppr = []
            for par in pars:
                ppr.append(self[par])
            try:
                ppr.append(self._extractFctVal(self.fitResult))
            except Exception as e:
                PE.warn(PE.PyAValError("The call to the `extractFctVal` function failed. Using full output." +
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
        self.model = saveModel
        return result

    def errorConfInterval(self, par, dstat=2.706, statTol=1e-2, hardLimit=None, maxiter=100, scale=None):
        """
        Calculate confidence interval for a parameter.

        This function uses linear extrapolation (similar
        to Newton's method) to find the points, where the
        objective function reaches an offset of `dstat`
        from the minimum value. 

        .. note:: You need to have carried out a fit before you can
                  use `errorConfInterval`.

        Parameters
        ----------
        par : string
            The parameter name for which to determine the error.
        dstat : float, optional
            The offset from the minimum to be reached in the objective
            function's value. For chi-square, a value of 1.0 corresponds
            to a one sigma (68%) confidence interval whereas the default
            of 2.706 refers to a 90% confidence interval.
        statTol : float, optional
            The acceptable (absolute) tolerance in achieving the
            target value in the objective function. The default is 0.01.
        hardLimit : list of two floats, optional
            Hard limits for the confidence interval. You can specify hard
            limits independently of parameter restrictions defined otherwise.
            To specify limits, specify a list of two floats defining
            the lower and upper limit. Use None where no limit is needed:
            e.g., use hardLimits=[0.0, None].
        maxiter : int, optional
            The maximum number of acceptable iterations. If exceeded,
            the calculation of the confidence interval is stopped and
            None is returned. The default is 100.
        scale : float, optional
            This number should reflect the typical scale or extent of
            the problem. If nothing is specified, a tenth of the current
            parameter value is used.

        Returns
        -------
        Confidence interval : dictionary
            The dictionary contains the following keys:
              - "limits" : A list of two floats specifying the lower and upper
                           limits of the confidence interval.
              - "OFVals" : A list of two floats containing the values of the
                           objective function at the lower and upper limit.
              - "OFMin" : The minimum value of the objective function.
              - "parMinVal" : The parameter value at the minimum value of the
                              objective function (best-fit value).
              - "iters" : The number of iterations needed to reach the result.
                          Note that the function returns None if the maximum
                          number of iterations is exceeded.
        """
        if not self._stepparEnabled:
            raise(PE.PyAOrderError("No data available. Did you carry out a fit already?",
                                   solutions="Fit the data first using `fit` method."))
        # Set up hard limit
        if hardLimit is None:
            hardLimit = [None, None]
        # Check whether hard limits are already violated
        if hardLimit[0] is not None:
            if self[par] < hardLimit[0]:
                raise(PE.PyAValError("The parameter value (" + str(self[par]) + ") is below the lower hard limit" +
                                     " of " + str(hardLimit[0]),
                                     solution="Check parameter values and hard limits."))
        if hardLimit[1] is not None:
            if self[par] > hardLimit[1]:
                raise(PE.PyAValError("The parameter value (" + str(self[par]) + ") is above the upper hard limit" +
                                     " of " + str(hardLimit[1]),
                                     solution="Check parameter values and hard limits."))
        # Set up scale parameter
        if scale is None:
            scale = abs(self[par] * 0.1)
        if scale == 0.0:
            raise(PE.PyAValError("The `scale` is 0.0.",
                                 solution="Provide a scale via the corresponding keyword."))
        # Save state of object
        saveObj = self.saveState()
        saveFitResult = self.fitResult
        saveModel = self.model.copy()
        # Get the minimum value of objective function
        # from fit.
        cvalmin = self.__extractFunctionValue(self.fitResult)
        parmin = self[par]
        iters = 2
        self.freeze(par)
        # Find lower limit
        # Calculate first step
        x0 = parmin
        y0 = cvalmin
        self[par] = parmin - scale
        # Note that mAA is dispensable, because self.minAlgo will be a callable.
        self.fit(None, None, yerr=None, minAlgo=self.minAlgo,
                 miniFunc=self.miniFunc)
        x1 = self[par]
        y1 = self.__extractFunctionValue(self.fitResult)
        while True:
            if y1 < cvalmin:
                PE.warn(PE.PyAParameterConflict(
                    "Found a new minimum for the objective function during error search."))
            iters += 1
            if iters > maxiter:
                PE.warn(PE.PyAValError("Maximum number if iterations (" + str(maxiter) + ") exceeded. Returning None.",
                                       solution=["Increase maximum number of iterations via `maxiter` flag.",
                                                 "Adapt scale parameter via `scale` keyword.",
                                                 "Use `steppar` to explore behavior of objective function."]))
                return None
            # Calculate parameters of straight line connecting points
            m = (y1 - y0) / (x1 - x0)
            b = y1 - m * x1
            # Extrapolate new guess from straight line
            newguess = ((dstat + cvalmin) - b) / m
            # Enforce hard limits
            atLimit = False
            if hardLimit[0] is not None:
                if newguess < hardLimit[0]:
                    newguess = hardLimit[0]
                    atLimit = True
            if hardLimit[1] is not None:
                if newguess > hardLimit[1]:
                    newguess = hardLimit[1]
                    atLimit = True
            self[par] = newguess
            self.fit(None, None, yerr=None, minAlgo=self.minAlgo, miniFunc=self.miniFunc,
                     *self.fminPars, **self.fminArgs)
            # Update guess
            x0 = x1
            y0 = y1
            x1 = self[par]
            y1 = self.__extractFunctionValue(self.fitResult)
            if abs(y1 - (dstat + cvalmin)) < statTol:
                break
            if atLimit and (abs(y1 - cvalmin) < dstat):
                break
        # Save lower limit and value of objective function
        lowerLimit = self[par]
        lowerOF = y1

        # Find upper limit
        x0 = parmin
        y0 = cvalmin
        self[par] = parmin + scale
        self.fit(None, None, yerr=None, minAlgo=self.minAlgo, miniFunc=self.miniFunc,
                 *self.fminPars, **self.fminArgs)
        x1 = self[par]
        y1 = self.__extractFunctionValue(self.fitResult)
        while True:
            iters += 1
            if iters > maxiter:
                PE.warn(PE.PyAValError("Maximum number if iterations (" + str(maxiter) + ") exceeded. Returning None.",
                                       solution=["Increase maximum number of iterations via `maxiter` flag.",
                                                 "Adapt scale parameter via `scale` keyword.",
                                                 "Use `steppar` to explore behavior of objective function."]))
                return None
            m = (y1 - y0) / (x1 - x0)
            b = y1 - m * x1
            newguess = ((dstat + cvalmin) - b) / m
            # Enforce hard limits
            atLimit = False
            if hardLimit[0] is not None:
                if newguess < hardLimit[0]:
                    newguess = hardLimit[0]
                    atLimit = True
            if hardLimit[1] is not None:
                if newguess > hardLimit[1]:
                    newguess = hardLimit[1]
                    atLimit = True
            self[par] = newguess
            self.fit(None, None, yerr=None, minAlgo=self.minAlgo, miniFunc=self.miniFunc,
                     *self.fminPars, **self.fminArgs)
            x0 = x1
            y0 = y1
            x1 = self[par]
            y1 = self.__extractFunctionValue(self.fitResult)
            if abs(y1 - (dstat + cvalmin)) < statTol:
                break
            # Check hard limit
            if atLimit and (abs(y1 - cvalmin) < dstat):
                break
        # Save upper limit and value of objective function
        upperLimit = self[par]
        upperOF = y1
        # Restore old state of object
        self.restoreState(saveObj)
        self.fitResult = saveFitResult
        self.model = saveModel
        return {"limits": [lowerLimit, upperLimit], "OFVals": [lowerOF, upperOF], "OFMin": cvalmin,
                "parMinVal": parmin, "iters": iters}


#  def annealFit(self, x, y, scheduler, yerr=None, miniFunc=None):
#    """
#      Carry out a fit using scipy's anneal function.
#    """
#    raise(PE.PyANotImplemented("This is under construction..."))
#    # Assign internal data properties.
#    self.x = x.copy()
#    self.y = y.copy()
#    if yerr is not None:
#      self.yerr = yerr.copy()
#    else:
#      self.yerr = None
#    # Determine function to be minimized
#    if miniFunc is None:
#      self.miniFunc = self.__chiSqr()
#    else:
#      self.miniFunc = miniFunc
#
#    pyaan = PyAAnneal(self)
#    pyaan.anneal(self.miniFunc, scheduler)
#
#    self.assignValue(pyaan.states["best"]["pars"])
#    self.evaluate(x)
#
#    print pyaan.whichStop

#    print "Results of annealing: "
#    print "  ", self.fitResult
#    # Set parameters to best fit
#    self.pars.setFreeParams(self.fitResult[0])
#    # Calculate the best-fit model
#    self.updateModel()


def sampleEMCEE(fpns, fv0, lnp, largs=None, nwalker=None, scales=None, sampleArgs=None, dbfile="chain.emcee", ps=None, emcp=None):
    """
    MCMC sampling from specific density using the emcee package.

    This function may be used to use emcee to sample from any user-specified
    density, which does not have to be normalized. The resulting Markov Chains
    can be analyzed using the trace analysis package.

    Parameters
    ----------
    fpns : list of strings
        Names of parameters for which Markov Chains are constructed. 
    fv0 : dictionary
        A dictionary mapping parameter name to starting value. This
        dictionary may contain any number of additional key-value pairs
    lnp : callable
        A function (or equivalent callable) which returns the (natural)
        logarithm of the (generally unnormalized) density. The
        first (and only mandatory) argument to the function is a dictionary
        holding the current parameter values for which the posterior
        density is to be evaluated. The function may take any number of
        additional keyword arguments, which can be specifies by the `largs`
        parameter.
    largs : dictionary, optional
        A set of additional arguments passed to the `lnp` callable.
    nwalker : int, optional
        The number of walker to be used. By default, four times the
        number of free parameters is used.
    scales : dictionary, optional
        The scales argument can be used to control the initial distribution
        of the walkers. By default, all walkers are distributed around the
        location given by the current state of the object, i.e., the current
        parameter values. In each direction, the walkers are randomly distributed
        with a Gaussian distribution, whose default standard deviation is one.
        The scales argument can be used to control the width of the Gaussians used
        to distribute the walkers.
    sampleArgs : dictionary, optional
        Number controlling the sampling process. Use 'burn' (int) to specify
        the number of burn-in iterations (default is 0). Via 'iters' (int)
        the numbers of iterations after the burn-in can be specified (default 1000).
        The 'process' (int) key can be used to control the number of iterations after
        which the progress bar is updated (default is iters/100). Note that the
        'progressbar' package must be installed to get a progress bar. Otherwise
        more mundane print statements will be used.  
    dbfile : string, optional
        The result of the sampling, i.e., the chain(s), the corresponding
        values of the posterior, and the names of the free parameters are
        saved to the specified file (by default 'chain.emcee' is used).
        The traces stored there can be analyzed using the 'TraceAnalysis'
        class. Set this parameter to 'None' to avoid saving the results.
    ps : tuple, optional
        A tuple holding the current position and state of the sampler. This
        tuple is returned by this method. The `ps` argument can be used
        to continue sampling from the last state. Note that no burn-in will
        ne carried out and the other arguments should be given as previously
        to continue sampling successfully. 
    emcp : dictionary, optional
        Extra arguments handed to `EnsembleSampler` object.

    Returns
    -------
    pos, state : state of emcee sample
        These information may be used to continue the sampling from previous position.
    """

    if not ic.check["emcee"]:
        raise(PE.PyARequiredImport("Could not import the 'emcee' package.",
                                   solution="Please install 'emcee'."))

    # Number of dimensions
    ndims = len(fpns)

    if ndims == 0:
        raise(PE.PyAValError("At least one free parameter is required for sampling.",
                             where="sampleEMCEE",
                             solution="TBD."))

    if not dbfile is None:
        if re.match(".*\.emcee$", dbfile) is None:
            PE.warn(PE.PyAValError("The db filename (" + str(dbfile) + ") does not end in .emcee. TraceAnalysis will not recognize it as an emcee trace file.",
                                   solution="Use a filename of the form *.emcee"))

    # Number of walkers
    if nwalker is None:
        nwalker = ndims * 4

    if nwalker < ndims * 2:
        raise(PE.PyAValError("The number of walkers must be at least twice the number of free parameters.",
                             where="sampleEMCEE",
                             solution="Increase the number of walkers."))
    if nwalker % 2 == 1:
        raise(PE.PyAValError("The number of walkers must be even.",
                             where="sampleEMCEE",
                             solution="Use an even number of walkers."))

    # Set default values for sampleArgs
    if sampleArgs is None:
        sampleArgs = {}
    if not "burn" in sampleArgs:
        sampleArgs["burn"] = 0
    if not "iters" in sampleArgs:
        sampleArgs["iters"] = 1000
    if not "progress" in sampleArgs:
        sampleArgs["progress"] = sampleArgs["iters"] / 100

    # Dictionary used to call the provided log posterior function
    kvs = copy.copy(fv0)

    # Manage None is argument to largs
    if largs is None:
        _largs = {}
    else:
        _largs = largs

    # Generate log posterior function required by emcee
    def logpost(x):
        """
        Calls user-specified log posterior function with
        properly updated parameter dictionary.
        """
        for i, n in enumerate(fpns):
            # Assign parameter values, which are variable
            kvs[n] = x[i]
        return lnp(kvs, **_largs)

    if emcp is None:
        emcp = {}

    if ps is None:

        # Generate the sampler
        emceeSampler = emcee.EnsembleSampler(nwalker, ndims, logpost, **emcp)

        if scales is None:
            scales = {}

        # Generate starting values
        pos = []
        for _ in smo.range(nwalker):
            pos.append(np.zeros(ndims))
            for i, n in enumerate(fpns):
                if not n in scales:
                    s = 1.0
                else:
                    s = scales[n]
                pos[-1][i] = np.random.normal(fv0[n], s)

        # Default value for state
        state = None

        if sampleArgs["burn"] > 0:
            # Run burn-in
            pos, prob, state = emceeSampler.run_mcmc(pos, sampleArgs["burn"])
            # Reset the chain to remove the burn-in samples.
            emceeSampler.reset()

    else:
        # Generate the sampler
        emceeSampler = emcee.EnsembleSampler(nwalker, ndims, logpost, **emcp)
        # Assign position and state from previous run
        pos, state = ps

    if (not sampleArgs["progress"] is None) and ic.check["progressbar"]:
        widgets = ['EMCEE progress: ', progressbar.Percentage(), ' ', progressbar.Bar(marker=progressbar.RotatingMarker()),
                   ' ', progressbar.ETA()]
        pbar = progressbar.ProgressBar(
            widgets=widgets, maxval=sampleArgs["iters"]).start()

    n = 0
    for pos, prob, state in emceeSampler.sample(pos, rstate0=state, iterations=sampleArgs["iters"], thin=1, storechain=True):
        n += 1
        if (not sampleArgs["progress"] is None) and (n % sampleArgs["progress"] == 0):
            if ic.check["progressbar"]:
                pbar.update(n)
            else:
                print("EMCEE: Reached iteration ",
                      n, " of ", sampleArgs["iters"])

    # Save the chain to a file
    if not dbfile is None:
        np.savez_compressed(open(dbfile, 'wb'), chain=emceeSampler.chain, lnp=emceeSampler.lnprobability,
                            pnames=np.array(fpns, dtype=np.unicode_))

    return pos, state
