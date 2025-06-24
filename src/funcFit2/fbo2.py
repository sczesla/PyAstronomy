# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
import re
import copy
import types
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import pyaC
import six
import six.moves as smo
import bidict
import collections
import sys
import scipy.optimize as sco
import inspect
import sys
import fnmatch
from PyAstronomy import modelSuite as _ms

_fbo2module = sys.modules[__name__]

from PyAstronomy.funcFit import _pymcImport, _scoImport, ic
from PyAstronomy import funcFit as _fuf

if ic.check["progressbar"]:
    import progressbar
if ic.check["emcee"]:
    import emcee


class PyAPa(object):
    
    def setVal(self, v):
        self._value = v
    
    def getVal(self):
        return self._value
    
    def setFree(self, v):
        if (v == True) and (not self.thawable):
            raise(PE.PyAValError("Cannot thaw requested parameter.", \
                                 solution="Some parameters cannot be thawed. See the documentation for details."))
        self._free = v
    
    def getFree(self):
        return self._free
    
    value = property(getVal, setVal)
    free = property(getFree, setFree)
    
    def __init__(self, value=0.0):
        self._value = value
        self._free = False
        self.affects = set()
        self.relation = None
        self.dependsOn = None
        self.thawable = True
        self.restriction = (None, None)
    
    def restrict(self, lower=None, upper=None):
        lin = not lower is None
        uin = not upper is None
        if lin+uin == 0:
            # No limits defined (or removed)
            pass
        elif lin+uin == 2:
            if lower >= upper:
                raise(PE.PyAValError("The lower limit has to be smaller than the upper bound in restriction."))
        self.restriction = (lower, upper)
    
    def isRestricted(self):
        return any([(not r is None) for r in self.restriction])
    
    def isRelated(self):
        return (not self.relation is None)
    
    def updateRelation(self):
        """
        Update value according to relation
        """
        if self.isRelated():
            self._value = self.relation.update()

      
class PyARelation(object):
    
    def __init__(self, ivs, func):
        """
        Parameters
        ----------
        ivs : list of Parameter instances
            The independent parameters
        func : Callable
            The functional relation
        """
        self.ivs = ivs
        self.func = func
    
    def update(self):
        return self.func(*[p.value for p in self.ivs])
    
    
class PyAPrior(object):
    
    def __init__(self, ivs, func, descr=None):
        """
        Parameters
        ----------
        ivs : list of Parameter instances
            The independent parameters
        func : Callable
            Function returning log(prior) depending on the value of the
            independent parameters
        descr : string, optional
            Description of the prior
        """
        if not hasattr(ivs, "__iter__"):
            ivs = [ivs]
        self.ivs = ivs
        self.func = func
        if descr is None:
            self.descr = "Pya Prior"
        else:
            self.descr = descr
    
    def evaluate(self):
        return self.func(*[p.value for p in self.ivs])
    
    def __str__(self):
        return str(self.descr)


class PyAUniformPrior(PyAPrior):
    
    def __init__(self, ivs, lower=None, upper=None, descr=""):

        if (lower is None) and (upper is None):
            raise(PE.PyAValError("At least one of lower and upper must be specified"))

        if (lower is None) and (upper is None):
            f = lambda x:0.0
        elif (lower is None) and (not upper is None):
            f = lambda x:0.0 if x <= upper else -np.inf
        elif (not lower is None) and (upper is None):
            f = lambda x:0.0 if x >= lower else -np.inf
        else:
            if upper <= lower:
                raise(PE.PyAValError("'lower' must be smaller than 'upper'", \
                                     where="addUniformPrior"))
            r = upper - lower
            def f(x):
                y = np.log(r) if (x >= lower) and (x <= upper) else -np.inf
                return y
        
        PyAPrior.__init__(self, ivs, f, descr=descr)


class PyANormalPrior(PyAPrior):
    
    def __init__(self, ivs, mean=None, std=None, descr=""):

        if (mean is None) and (std is None):
            raise(PE.PyAValError("Both mean and std must be specified"))

        const = -0.5*np.log(2*np.pi*std**2)
        durch = 2*std**2

        f = lambda x: const - (x-mean)**2/durch
        
        PyAPrior.__init__(self, ivs, f, descr=descr)

        
class PyASmoothUniformPrior(PyAPrior):
      
    def __init__(self, ivs, lower=None, upper=None, scale=1e-9, descr=""):
        
        if (lower is None) and (upper is None):
            raise(PE.PyAValError("At least one of lower and upper must be specified"))      
      
        pih = np.pi/2.0
        flmax = sys.float_info.max/1e3
        def f(x):
            if not lower is None:
                if x < lower:
                    r = -np.arctan((lower-x)/scale)/pih*flmax
                    return r
            if not upper is None:
                if x > upper:
                    r = -np.arctan((x-upper)/scale)/pih*flmax
                    return r 
            return 0.0
        
        PyAPrior.__init__(self, ivs, f, descr=descr)
      
      
class PyAScalePrior(PyAPrior):
    
    def __init__(self, ivs, m=-1., descr=""):

        def f(x):
            if x <= 0.0:
                return -np.inf
            return m*np.log(x)
        
        PyAPrior.__init__(self, ivs, f, descr=descr)


class PyACauchyPrior(PyAPrior):
    
    def __init__(self, ivs, scale, x0=0., side="both", descr=""):

        def f(x):
            
            pdf = np.log(scale) - np.log(np.pi) - np.log(scale**2 + (x-x0)**2)
            if side == "both":
                return pdf
            elif side == "positive":
                if (x < x0):
                    return -np.inf
                # An extra factor accounts for the 'missing' other side
                return np.log(2.) + pdf
            elif side == "negative":
                if (x > x0):
                    return -np.inf
                # An extra factor accounts for the 'missing' other side
                return np.log(2.) + pdf
            else:
                raise(PE.PyAValError("Unknown side parameter in PyACauchyPrior (" + str(side) + ")", \
                                     solution="Use either of 'both', 'positive', 'negative'"))
        
        PyAPrior.__init__(self, ivs, f, descr=descr)
 
      
class PyABPS(object):
    """
    Basic Parameter Set
    """
    
    def addParam(self, name, value=0.0):
        """
        Add parameter to the list
        
        Parameters
        ----------
        name : string
            Name of the parameter
        value : optional
            Default value of the parameter
        """
        if name in self.pmap:
            raise(PE.PyAValError("Name '" + str(name) + "' already exists!", \
                                 where="addPyAPS::Param"))
        self.pmap[name] = PyAPa(value)
    
    def _checkParam(self, n):
        """ Throws and exception if parameter 'n' does not exist """
        if not n in self.pmap:
            raise(PE.PyAValError("No such parameter: " + str(n),
                                 solution="Use one of: " + ', '.join(list(self.pmap))))
    
    def __getitem__(self, n, ref=False):
        """ Get reference to PyaPa """
        self._checkParam(n)
        if ref:
            return self.pmap[n]
        else:
            return self.pmap[n].value
    
    def __setitem__(self, n, v):
        """ Set value and update related if necessary """
        self._checkParam(n)
        self.pmap[n].value = v
        if len(self.pmap[n].affects) > 0:
            for a in self.pmap[n].affects:
                a.updateRelation()

    def copy(self):
        """ Return a shallow copy of the object """
        c = PyABPS([], self.rootName, number=self.number)
        c.pmap = copy.copy(self.pmap)
        return c

    def __init__(self, pns, rootName, number=0):
        self.rootName = rootName
        self.number = number
        self.pmap = bidict.OrderedBidict()
        for n in pns:
            self.addParam(n, 0.0)
        
      
      
class PyABaSoS(object):
    """
    Basic Set of Basic Parameter Sets
    """
    
    def _allRoots(self):
        """ Get root names of all basic sets """
        return [s.rootName for s in self.bss]
    
    def _rootNumbers(self, rn):
        """ All number associated with a certain root """
        return set([b.number if b.rootName == rn else 0 for b in self.bss])
    
    def _updatepmap(self):
        """ Update parameter names according to root name and number. """
        ar = self._allRoots()
        if len(ar) == 1:
            # Only one subset. Apply no name updating
            self.pmap = bidict.OrderedBidict()
            for n in list(self.bss[0].pmap):
                self.pmap[n] = self.bss[0].pmap[n]
        else:
            # More than one subset
            # Root Count
            rc = {r:ar.count(r) for r in set(ar)}
            self.pmap = bidict.OrderedBidict()
            for b in self.bss:
                useNumber = int(rc[b.rootName] > 1)
                for n in list(b.pmap):
                    useName = self._composePN(n, b.rootName, b.number*useNumber)
                    self.pmap[useName] = b.pmap[n]
    
    def addBPS(self, s):
        """
        Add Basic Parameter Set to set of sets
        
        Parameters
        ----------
        s : PyABPS
            Instance of PyABPS
        """
        snew = s.copy()
        rns = self._rootNumbers(snew.rootName)
        if len(rns) == 0:
            snew.number = 1
        else:
            snew.number = max(rns)+1
        self.bss.append(snew)
        self._updatepmap()
    
    def _checkParam(self, n):
        """ Throws and exception if parameter 'n' does not exist """
        if not n in self.pmap:
            raise(PE.PyAValError("No such parameter: " + str(n),
                                 solution="Use one of: " + ', '.join(list(self.pmap))))
    
    def _lazyNameMatch(self, n):
        """
        Search by 'part of name'
        
        Searches the parameter map for occurrences of the specified
        sub-string, tests whether it is unique, and returns the full
        valid name if so. 
        
        Parameters
        ----------
        n : string
            (Unique) Part of the parameter name
            
        Returns
        -------
        True name : string
        """
        try:
            self._checkParam(n)
            return n
        except PE.PyAValError as e:
            matches = [pn  for pn in list(self.pmap) if pn.find(n) != -1]
            if len(matches) == 1:
                return matches[0]
            elif len(matches) > 1:
                raise(PE.PyAValError("More than one match in lazy name matching for parameter '" + str(n) + "'. Matches: " + ', '.join(matches)))
            else:
                raise(PE.PyAValError("No match in lazy name matching for parameter '" + str(n) + "'."))
            
    
    def __getitem__(self, n):
        """ Get value of parameter """
        n = self._lazyNameMatch(n)
        return self.pmap[n].value
    
    def getPRef(self, n):
        """ Get reference to PyaPa """
        self._checkParam(n)
        return self.pmap[n]
    
    def __setitem__(self, n, v):
        """ Set value and update related if necessary """
        n = self._lazyNameMatch(n)
        self.pmap[n].value = v
        if len(self.pmap[n].affects) > 0:
            for a in self.pmap[n].affects:
                a.updateRelation()

    def rename(self, old, new):
        """
        Rename parameter
        
        Parameters
        ----------
        old, new : strings
            Old and new name
        """
        self._checkParam(old)
        if new in self.pmap:
            raise(PE.PyAValError("Parameter " + str(new) + " already exists.", \
                                 solution="Use a unique name."))
        tmp = self.pmap[old]
        del self.pmap[old]
        self.pmap[new] = tmp
        
    def relate(self, dv, idv, func=None):
        """
        Establish functional relations
        
        Parameters
        ----------
        dv : string
            Dependent variable
        idv : string or list of strings
            Independent variables
        func : callable
            Callable taking the independent variables as arguments
            and returns the value of the dependent variable.
        """
        if isinstance(idv, six.string_types):
            # Convert single string into list
            idv = [idv]
        # Check all parameters
        self._checkParam(dv)
        [self._checkParam(n) for n in idv]
        
        # By default, use equal
        if func is None:
            func = lambda x:x
        
        # Manage dependent variable
        self.pmap[dv].dependsOn = [self.pmap[n] for n in idv]
        self.pmap[dv].relation = PyARelation(self.pmap[dv].dependsOn, func)
        # Manage affected variables
        for n in idv:
            self.pmap[n].affects.update([self.pmap[dv]])
        # Trigger update of the value of the related parameter
        self.pmap[dv].updateRelation()
        # Freeze the dependent variable
        self.freeze(dv)
    
    def copy(self):
        """ Return a shallow copy of the object """
        c = PyABaSoS()
        c.pmap = copy.copy(self.pmap)
        c.bss = copy.copy(self.bss)
        return c
    
    @classmethod
    def combine(cls, left, right):
        """
        Combine parameter sets
        
        Parameters
        ----------
        left, right : PyABaSoS
            Instances of PyABaSoS holding the individual sets
        
        Returns
        -------
        Combined object : PyABaSoS
        """
        c = cls()
        c.bss = []
        for b in left.bss:
            c.addBPS(b)
        for b in right.bss:
            c.addBPS(b)
        return c
    
    def _composePN(self, b, r, n):
        """
        Compose parameter name
        
        Parameters
        ----------
        b, r : strings
            Base and root
        n : int
            Number
        
        Returns
        -------
        Name : string
            Complete parameter name
        """
        result = b
        if (not r is None) and (r != ""):
            result += "_" + r
        if (not n is None) and (n != 0):
            result += ("(%d)" % n)
        return result
    
    def _decomposePN(self, n):
        """
        Decompose parameter name into Base, Root, Number
        """
        r = re.match("([^_]+)((_([^\(]+)?)?(\(([0-9]+)\))?)?", n)
        if r is None:
            return None, None, None
        else:
            return r.group(0), r.group(3), r.group(5)
   
    def parameters(self):
        """
        Get all parameter names and values
        
        Returns
        -------
        Parameters : OrderedDict (collections)
            (Ordered) Dictionary mapping parameter name to value
        """
        ps = collections.OrderedDict()
        for k, v in six.iteritems(self.pmap):
            ps[k] = v.value
        return ps
    
    def freeParamNames(self):
        """ Get list of names of free parameters """
        return [p for p in list(self.pmap) if self.pmap[p].free]
    
    def freeParameters(self):
        """
        Get a dictionary mapping parameter name to value
        """
        fpn = self.freeParamNames()
        fpv = self.freeParamVals()
        return dict(zip(fpn, fpv))
    
    def getRestrictions(self):
        """ Get dictionary mapping parameter names to restrictions """
        return {p:self.pmap[p].restriction for p in list(self.pmap) if self.pmap[p].isRestricted()}
    
    def setRestriction(self, restrictions):
        """
        Add restrictions for parameters
        
        Parameters
        ----------
        restrictions : dictionary
            Maps parameter to restriction. A restriction is a two-tuple
            with lower and upper limit. If None is specified, no restriction
            applies in this direction.
        """
        for k, v in six.iteritems(restrictions):
            self._checkParam(k)
            self.getPRef(k).restrict(lower=v[0], upper=v[1])
    
    def freeParamVals(self):
        """ Get list of names of free parameters """
        return [self.pmap[p].value for p in list(self.pmap) if self.pmap[p].free]
    
    def setFreeParamVals(self, vals):
        """
        Assign values to free parameters.
        
        The order is fixed (use, e.g., freeParamNames to get it).
        
        Parameters
        ----------
        vals : list or array
            The values to be assigned
        """
        try:
            # Handle 'singleton' arrays
            if len(vals.shape) == 0:
                vals = np.array([vals])
        except:
            pass
        for i, p in enumerate(self.freeParamNames()):
            self[p] = vals[i]
    
    def _thawfreeze(self, pns, free):
        """
        Thaw or freeze parameters
        
        Parameters
        ----------
        pns : string or list of strings
            Relevant parameters
        free : boolean
            True to thaw and False to freeze
        """
        if isinstance(pns, six.string_types):
            pns = [pns]
        for p in pns:
            self._checkParam(p)
            self.pmap[p].free = free
    
    def thaw(self, pns):
        """
        Thaw parameters
        
        Parameters
        ----------
        pns : string or list of strings
            Relevant parameters
        """
        self._thawfreeze(pns, True)
    
    def freeze(self, pns):
        """
        Freeze parameters
        
        Parameters
        ----------
        pns : string or list of strings
            Relevant parameters
        """
        self._thawfreeze(pns, False)
    
    def __init__(self, *args):
        self.bss = []
        for s in args:
            self.addBPS(s)
        self._updatepmap()
      


def gaussLogL(x, y, yerr, m):
    """
    Gaussian likelihood
    
    Parameters
    ----------
    x, y : arrays
        The x and y coordinates of the data points.
    yerr : array or float
        Uncertainty of data points
    m : array
        The model
    
    Returns
    -------
    lnl : float
        The natural logarithm of the likelihood.
    """
    if hasattr(yerr, "__iter__"):
        # yerr is an array
        lnl = -len(x)/2.0*np.log(2.*np.pi) - np.sum(np.log(yerr)) - 0.5 * np.sum((m-y)**2/(yerr**2))
        return lnl
    else:
        # yerr is a float
        lnl = -len(x)/2.0*np.log(2.*np.pi) - len(x)*np.log(yerr) - 0.5 * np.sum((m-y)**2/(yerr**2))
        return lnl
    


def _gaussLogL(self, *args, **kwargs):
    """
    Default implementation of Gaussian likelihood.
    
    Parameters
    ----------
    x, y : arrays
        The x and y coordinates of the data points.
    yerr : array or float, optional
        If not specified, a value of 1.0 will be assumed. Otherwise interpreted as
        the error of the data points.
    
    Returns
    -------
    lnl : float
        The natural logarithm of the likelihood based on a model f(x), data points (x,y), and
        normal errors, yerr.
    """
    if len(args) == 2:
        x, y = args[0], args[1]
        yerr = 1.0
    elif len(args) >= 3:
        x, y, yerr = args[0], args[1], args[2]
    else:
        raise(PE.PyAValError("Invalid call to _gaussLogL. Received " + str(len(args)) + " arguments but takes 2 or 3 (x, y, [yerr])."))

    if "_currentModel" in kwargs:
        m = kwargs["_currentModel"]
    else:
        m = self.evaluate(*args, **kwargs)
    
    return gaussLogL(x, y, yerr, m)

        
def getchisqrobjf(sumfct=np.sum):
    
    def chisqrobjf(self, pars, *args, **kwargs):
        """
        Default implementation of chi-square objective function
        
        Parameters
        ----------
        x, y : arrays
            The x and y coordinates of the data points.
        yerr : array or float, optional
            If not specified, a value of 1.0 will be assumed. Otherwise interpreted as
            the error of the data points.
        
        Returns
        -------
        chi square : float
            Returns chi square if uncertainty is specified and sum of squared residuals
            otherwise.
        """
        if len(args) == 2:
            x, y = args[0], args[1]
            yerr = 1.0
        elif len(args) >= 3:
            x, y, yerr = args[0], args[1], args[2]
        else:
            raise(PE.PyAValError("Invalid call to _chisqr. Received " + str(len(args)) + " arguments but takes 2 or 3 (x, y, [yerr])."))
    
        if "_currentModel" in kwargs:
            m = kwargs["_currentModel"]
        else:
            m = self.evaluate(*args, **kwargs)
        
        return sumfct( (m-y)**2/yerr**2 )
    return chisqrobjf
    

class MBO(object):
    """
    Model Base Object
    
    Concepts
    --------
    
    Parameter names: Internal vs. external
        
    
    
    """
    
    def __init__(self, pars=None, rootName="", **kwargs):
        
        if pars is None:
            raise(PE.PyAValError("You need to specify the parameter names via 'pars'.", \
                                 where="MBO2", \
                                 solution="Specify something along the lines of 'pars = ['pn1', 'pn2', ...]'."))
        
        self.pars = (PyABaSoS(PyABPS(pars, rootName)))
        self.priors = []
        
        # Penalty factor used in objective function
        self._pf = 1e20
        
        self.rootName = rootName
        self.leftCompo = None
        self.rightCompo = None
     
    def _combineMBOs(self, right):
        """
        Combine two MBO2s into a new one

        Parameters
        ----------
        right : MBO2 instance
            MBO object on the right hand side of the operation
        
        Returns
        -------
        combination : MBO2
        """
        r = MBO([], rootName="combined")
        r.pars = PyABaSoS.combine(self.pars, right.pars)
        r.leftCompo = self
        r.rightCompo = right
        return r
     
    def __add__(self, right):
        result = self._combineMBOs(right)
        result.evaluate = types.MethodType(lambda self, *args, **kwargs: \
                                           self.leftCompo.evaluate(*args, **kwargs) + self.rightCompo.evaluate(*args, **kwargs), result)
        return result
    
    def __sub__(self, right):
        result = self._combineMBOs(right)
        result.evaluate = types.MethodType(lambda self, *args, **kwargs: \
                                           self.leftCompo.evaluate(*args, **kwargs) - self.rightCompo.evaluate(*args, **kwargs), result)
        return result
    
    def __mul__(self, right):
        result = self._combineMBOs(right)
        result.evaluate = types.MethodType(lambda self, *args, **kwargs: \
                                           self.leftCompo.evaluate(*args, **kwargs) * self.rightCompo.evaluate(*args, **kwargs), result)
        return result
    
    def __div__(self, right):
        result = self._combineMBOs(right)
        result.evaluate = types.MethodType(lambda self, *args, **kwargs: \
                                           self.leftCompo.evaluate(*args, **kwargs) / self.rightCompo.evaluate(*args, **kwargs), result)
        return result
    
    def __truediv__(self, right):
        result = self._combineMBOs(right)
        result.evaluate = types.MethodType(lambda self, *args, **kwargs: \
                                           self.leftCompo.evaluate(*args, **kwargs) / self.rightCompo.evaluate(*args, **kwargs), result)
        return result

    def __pow__(self, right):
        result = self._combineMBOs(right)
        result.evaluate = types.MethodType(lambda self, *args, **kwargs: \
                                           self.leftCompo.evaluate(*args, **kwargs) ** self.rightCompo.evaluate(*args, **kwargs), result)
        return result
             
    def parameterSummary(self, toScreen=True, prefix=""):
        """
        Parameter summary
        
        Parameters
        ----------
        toScreen : boolean, optional
            If True (default), the parameter summary is printed to screen.
        prefix : string, optional
            If specified, the prefix is added at the front of all lines.
        
        Returns
        -------
        Summary : list of lines
            The summary as a list of lines
        """
        # Max length of parameter name
        mpl = max([len(n) for n in list(self.parameters())])
        lines = []
        tf = {True:"T", False:"F"}
        for k, v in self.pars.pmap.items():
            l = prefix + "    " + ("%" + str(mpl) + "s") % k + " = % 12g" % v.value
            l += ", free: " + tf[v.free] + ", restricted: " + tf[v.isRestricted()]
            l += ", related: " + tf[v.isRelated()]
            lines.append(l)
#         lines.append(prefix + "Priors: ")
#         if len(self.priors) > 0:
#             for i, p in enumerate(self.priors, 1):
#                 lines.append(prefix + " "*4 + ("%2d) " % i) + str(p))
#         else:
#             lines.append("    All uniform (=1)")
        
        if toScreen:
            mll = max([len(l) for l in lines])
            print("-" * (mll//3) + " Parameter summary " + "-"*(mll-mll//3-19) )
            for l in lines:
                print(l)
            print("-" * mll)
        return lines
        
    def relate(self, dv, idv, func=None):
        self.pars.relate(dv, idv, func)
    
    def __getitem__(self, n):
        return self.pars[n]
    
    def __setitem__(self, n, v):
        nns = self._pex(n)
        for nn in nns: 
            self.pars[nn] = v
    
    def getPRef(self, n):
        return self.pars.getPRef(n)
    
    def freeParamVals(self):
        return self.pars.freeParamVals()
    
    def freeParamNames(self):
        return self.pars.freeParamNames()
    
    def freeParameters(self):
        return self.pars.freeParameters()
    
    def setFreeParamVals(self, vals):
        self.pars.setFreeParamVals(vals)
        
    def assignValues(self, vals):
        """
        Assign values to parameters
        
        Parameters
        ----------
        vals : dictionary
            Maps parameter name to new value
        """
        for k, v in six.iteritems(vals):
            self.pars[k] = v
    
    def _pex(self, pns):
        """
        Apply unix filename-like pattern matching and return updated list of parameters
        
        Returns
        -------
        updated list : list of strings
            Possibly expanded list of parameter names. For parameter names without expansion (and
            possibly without match in the list of model parameter names), the input string is returned.
        """
        # Ensure it is a list of strings
        if isinstance(pns, six.string_types):
            pns = [pns]
        r = []
        for p in pns:
            e = fnmatch.filter(list(self.parameters()), p) 
            r.extend(e if len(e)>0 else [p])
        return r
       
    def thaw(self, pns):
        pns = self._pex(pns)
        self.pars.thaw(pns)
        
    def freeze(self, pns):
        pns = self._pex(pns)
        self.pars.freeze(pns)
    
    def setRestriction(self, restricts, ap=True, usesmooth=False, sscale=1e-9):
        """
        Restrict parameter ranges
        
        Sets lower and upper bounds for parameters where desired. The bounds
        are saved to be used in constrained optimization, e.g., via
        fmin_cobyla. Additionally, a uniform prior is added to take effect
        on the posterior probability distribution. If desired, the edges
        of the priors are 'smoothed' to alleviate convergence problems with
        unconstrained optimization algorithms working on the posterior.
        
        Parameters
        ----------
        restricts : dictionary
            Maps parameter name to resctriction, i.e., a two-tuple giving
            lower and upper bound. Use None where no bound is desired.
        ap : boolean, optional
            Apply priors? If True (default), priors will be defined according
            to the restrictions.
        usesmooth : boolean, optional
            If True, a uniform prior with 'smoothed' edges is used instead
            if a sharply-edged uniform prior. The scale if the edge is given
            by `sscale`. Default is False.
        sscale : float, optional
            Scale of the edge used for smooth uniform prior (if desired).
        """
        self.pars.setRestriction(restricts)
        if ap:
            # Add prior(s)
            for k, v in six.iteritems(restricts):
                self.pars._checkParam(k)
                if any([not b is None for b in v]):
                    # Discard cases were both limits are none
                    if not usesmooth:
                        self.addUniformPrior(k, lower=v[0], upper=v[1])
                    else:
                        self.addSmoothUniformPrior(k, lower=v[0], upper=v[1], scale=sscale)
    
    def getPenaltyFactor(self):
        return self._pf
    
    def setPenaltyFactor(self, pf):
        if pf < 0:
            raise(PE.PyAValError("Penalty factor must not be negative."))
        self._pf = pf
            
    def getRestrictions(self):
        return self.pars.getRestrictions()
    
    def parameters(self):
        return self.pars.parameters()
    
    def availableParameters(self):
        """
        Provides a list of existing parameter names.

        Returns
        -------
        Parameters : list of strings 
            A list with the names of available parameters.
        """
        return list(self.parameters())

    def frozenParameters(self):
        """
        Get names and values of frozen parameters.

        Returns
        -------
        Frozen parameters: dict
             Dictionary containing the names and values
             of all frozen parameters ({"parName":value, ...}).
        """
        fps = {}
        fpns = self.freeParamNames()
        for k, v in six.iteritems(self.parameters()):
            if not k in fpns:
                fps[k] = v
        return fps
    
    def setlogL(self, logl):
        """ Assign logL method """
        
        if isinstance(logl, six.string_types): 
            # Specification by string
            if logl.lower() == "1dgauss":
                self.logL = types.MethodType(_gaussLogL, self)
            else:
                raise(PE.PyAValError("Unknown logl specification: " + str(logl)))
        elif hasattr(logl, "__call__"):
            # It is a callable
            self.logL = types.MethodType(logl, self)
        elif logl is None:
            # Do nothing here
            pass
        else:
            raise(PE.PyAValError("logl is neither a string nor a callable nor None.", \
                                 where="MBO2"))
        
    def logL(self, *args, **kwargs):
        raise(PE.PyANotImplemented("""
The log(likelihood) method 'logL' needs to be implemented.

The likelihood function calculates the (natural) logarithm of the likelihood for the
current set of parameter values. It takes no obligatory arguments. Usually, however,
it requires a data set as input.
        """))

    def evaluate(self, *args, **kwargs):
        pass

    def logPrior(self, *args, **kwargs):
        """ Returns natural logarithm of prior for current parameter values """
        result = 0.0
        for p in self.priors:
            result += p.evaluate()
            if np.isinf(result):
                return result
        return result
    
    def _too(self, x, d):
        """ Returns d if x is None and x otherwise """
        if x is None:
            return d
        return x
    
    def addUniformPrior(self, pn, lower=None, upper=None):
        """
        Add a uniform prior
        
        Parameters
        ----------
        pn : string
            Name of the parameter
        lower, upper : float, optional
            Lower and upper bounds
        """
        self.pars._checkParam(pn)
        descr = "Uniform prior on '" + str(pn) + "' (lower = %g, upper = %g)" % (self._too(lower, -np.inf), self._too(upper, np.inf))
        self.addPrior(PyAUniformPrior(self.getPRef(pn), lower=lower, upper=upper, descr=descr))
      
    def addSmoothUniformPrior(self, pn, lower=None, upper=None, scale=1e-9):
        """
        Add a uniform prior with 'smoothed' (arctan) edges
        
        Parameters
        ----------
        pn : string
            Name of the parameter
        lower, upper : float, optional
            Lower and upper bounds
        """
        self.pars._checkParam(pn)
        descr = "Smooth uniform prior on '" + str(pn) + "' (lower = %g, upper = %g, scale = %g)" % \
                (self._too(lower, -np.inf), self._too(upper, np.inf), scale)
        self.addPrior(PyASmoothUniformPrior(self.getPRef(pn), lower=lower, upper=upper, descr=descr))
        
    def addScalePrior(self, pn, m):
        """
        Add scale prior of the form x**m
        
        Parameters
        ----------
        pn : string
            Name of the parameter
        m : float
            Exponent (e.g., -1 for x**-1 prior distribution)
        """
        self.pars._checkParam(pn)
        descr = "Scale prior on '" + str(pn) + "' (exponent = %g)" % m
        self.addPrior(PyAScalePrior(self.getPRef(pn), m, descr=descr))
        
    def addCauchyPrior(self, pn, x0, scale, side):
        """
        Add Cauchy prior of the form scale / (pi*(scale**2 + (x-x0)**2)).
        
        The density is multiplied by two of a one-sided version is desired.
        
        Parameters
        ----------
        pn : string
            Name of the parameter
        x0 : float
            Location parameter of the distribution
        scale : float
            Scale parameter of the distribution
        side : string, {both, positive, negative}
            Use 'positive' to rule out values lower than x0 and negative for
            the reverse. Use 'both' to allow both sides of the center.
        """
        self.pars._checkParam(pn)
        descr = "Cauchy prior on '" + str(pn) + "' (scale = %g, location = %g, side = %s)" % (scale, x0, side)
        self.addPrior(PyACauchyPrior(self.getPRef(pn), scale, x0, side, descr=descr))
        
    def addNormalPrior(self, pn, mean, std):
        """
        Add normal (Gaussian) prior
        
        Parameters
        ----------
        pn : string
            Name of the parameter
        mean : float
            Location parameter of the distribution
        std : float
            Standard deviation of the distribution
        """
        self.pars._checkParam(pn)
        descr = "Normal prior on '" + str(pn) + "' (mean = %g, std = %g)" % (mean, std)
        self.addPrior(PyANormalPrior(self.getPRef(pn), mean, std, descr=descr))
        
    def addPrior(self, pyaprior):
        """
        Add a prior
        
        Parameters
        ----------
        pyaprior : PyAPrior
            Instance of PyAPrior
        """
        self.priors.append(pyaprior)
    
    def mlD(self, *args, **kwargs):
        """
        Get marginal likelihood of the data (if known)
        
        Returns
        -------
        mlD : float
            Natural logarithm of the likelihood of the data if known
            and None otherwise.
        """
        return None
        
    def logPost(self, *args, **kwargs):
        """
        Get natural logarithm of the posterior
        
        Parameters
        ----------
        alllogs : boolean, optional
            If given and set True, also a tuple is returned, holding the logarithm
            of the priors and the likelihood separately.
        
        Returns
        -------
        lnp : float
            The (natural) logarithm of the posterior. The posterior is normalized if
            the marginal likelihood of the data is known and specified by the `mlD`
            method.
        pandl : tuple of floats, optional
            Tuple holding the natural logarithms of the priors and the likelihood.
            Only returned if alllogs is given as a keyword argument and True.
        """
        
        lp = self.logPrior(*args, **kwargs)
        if np.isinf(lp):
            if ("alllogs" in kwargs) and (kwargs["alllogs"] == True):
                return -np.inf, (-np.inf, -np.inf)
            else:
                return -np.inf
        
        ll = self.logL(*args, **kwargs)
        md = self.mlD(*args, **kwargs)
        po = lp + ll  - (md or 0.0)
        
        if ("alllogs" in kwargs) and (kwargs["alllogs"] == True):
            return po, (lp, ll)
        
        return po

    def emceeLogPost(self, *args, **kwargs):
        """
        Get the logarithm of the posterior assuming that the first parameter is an array of free parameter values
        
        Sets the values of the free parameters to the specified values
        and evaluates the posterior probability. An implementation which
        can easily combined with sampling by emcee.
        
        Parameters
        ----------
        P : array
            An array holding values for the free parameters
        
        Returns
        -------
        lp : float
            Natural logarithm of the posterior
        """
        self.setFreeParamVals(args[0])
        kwargs["alllogs"] = True
        r = self.logPost(*args[1:], **kwargs)
        return r
  
    def addSPLikeObjf(self, f, name):
        """
        Add a SciPy-Like (SPL) objective function to the instance.
        
        Parameters
        ----------
        f : callable or string
            The scipy-like objective function to be assigned. A scipy-like objective function is one
            which takes as its first argument an array holding the values of the free parameters. The
            function may take any number of additional arguments and keywords.
        """
        
        if isinstance(f, six.string_types):
            if f == "-logl":
                def nln(self, *args, **kwargs):
                    """ -natural_log(Likelihood) """
                    return -self.logL(*args[1:], **kwargs)
                f = nln
            elif f == "-logpost":
                def nln(self, *args, **kwargs):
                    """ -natural_log(Posterior) """
                    return -self.logPost(*args[1:], **kwargs)
                f = nln
            elif f == "chisqr":
                f = getchisqrobjf()
            elif f == "nanchisqr":
                f = getchisqrobjf(sumfct=np.nansum)
            else:
                raise(PE.PyAValError("Unknown objective function string: '"+f+"'"))
        
        def objf(self, *args, **kwargs):
            # Make sure the current parameters are assigned
            self.setFreeParamVals(args[0])
            v = f(self, *args, **kwargs)
            if not np.isfinite(v):
                PE.PyAValError("Infinite value encountered in objective function for parameters: " + str(args[0]) + ", free parameters : " + ",".join(self.freeParamNames()))
            v += self.getRPenalty()
            return v
        
        objf.__doc__ = f.__doc__
        setattr(self, name, types.MethodType(objf, self))

    def grad(self, *args, **kwargs):
        raise(PE.PyANotImplemented("To use derivatives, implement the 'grad' function."))
      
    def getRPenalty(self):
        """
        Get penalty for violating restrictions
        
        Parameters
        ----------
        pf : float
            If a boundary, b, specified in a restriction is violated, the
            penalty is calculated as abs(b-v) x pf, where v is the current
            parameter value.
        
        Returns
        -------
        Penalty : float
            Penalty for current restrictions and parameter values
        """
        if self._pf == 0:
            return 0.0
        x = 0.0
        for p, r in six.iteritems(self.getRestrictions()):
            if (not r[0] is None) and (self[p] < r[0]):
                x += self._pf*abs(self[p]-r[0])
                continue
            if (not r[1] is None) and (self[p] > r[1]):
                x += self._pf*abs(self[p]-r[1])
                continue
        return x
    
#     def objfPenalize(self, pf=1e20):
#         """
#         Add restriction penalties to objective function
#         """
#         self._nonpenobjf = self.getSPLikeObjf()
#         def pobj(self, *args, **kwargs):
#             x = self._nonpenobjf(*args, **kwargs)
#             x += self.getRPenalty()
#             return x
#         pobj.__doc__ = self._nonpenobjf.__doc__ + " (penalized)"
#         self.objf = pobj
        

class MBOEv(MBO):
    """
    Model Base Object with an evaluate method and a predefined chi square objective function
    """
    
    def __init__(self, pars=None, rootName="", **kwargs):
        
        if pars is None:
            raise(PE.PyAValError("You need to specify the parameter names via 'pars'.", \
                                 where="MBOEv", \
                                 solution="Specify something along the lines of 'pars = ['pn1', 'pn2', ...]'."))
        
        MBO.__init__(self, pars, rootName, **kwargs)
        
        # Use likelihood based on Gaussian errors with std yerr
        self.setlogL("1dgauss")
        self.addSPLikeObjf("chisqr", "chisqr")
        self.addSPLikeObjf("nanchisqr", "nanchisqr")
    
    def evaluate(self, *args, **kwargs):
        """
        Evaluate model
        
        Parameters
        ----------
        
        """
        pass
    
    def _combineMBOs(self, right):
        """
        Combine two MBOs into a new one

        Parameters
        ----------
        right : MBO2 instance
            MBO object on the right hand side of the operation
        
        Returns
        -------
        combination : MBO2
        """
        r = MBOEv([], rootName="combined")
        r.pars = PyABaSoS.combine(self.pars, right.pars)
        r.leftCompo = self
        r.rightCompo = right
        return r
     



def _introdefarg(f, **kwargs):
    """
    Analyze and update keyword arguments of callable.
    
    Parameters
    ----------
    f : callable
        A function to be inspected (e.g., sco.fmin)
    kwargs : dictionary
        Dictionary of arguments to be updated.
    
    Returns
    -------
    kwargs, original kwargs : dictionaries
        Updated kwargs (including specifications from kwargs) and original
        kwargs as provided by the callable.
    """
    # Get dictionary of default arguments and default values
    fi = inspect.getargspec(f)
    odefargs = dict(zip(fi.args[-len(fi.defaults):],fi.defaults))
    
    defargs = odefargs.copy()
    # Update entries for which kwargs holds values
    for k in list(defargs):
        if k in kwargs:
            defargs[k] = kwargs[k]
    return defargs, odefargs


def fitfmin(m, objf, x, y, yerr=None, **kwargs):
    """
    Use scipy's fmin to fit model.
    """
    # Get keywords and default arguments
    defargs, _ = _introdefarg(sco.fmin, **kwargs)
    
    if not yerr is None:
        defargs["args"] = (x, y, yerr)
    else:
        defargs["args"] = (x, y)
    
    defargs["full_output"] = True
    
    fr = sco.fmin(objf, m.freeParamVals(), **defargs)
    m.setFreeParamVals(fr[0])
    return fr


def fitfmin_cobyla(m, objf, x, y, yerr=None, cons=None, **kwargs):
    """
    Use scipy's fitfmin_cobyla to fit model.
    """
    # Get keywords and default arguments
    defargs, _ = _introdefarg(sco.fmin_cobyla, **kwargs)
    
    if not yerr is None:
        defargs["args"] = (x, y, yerr)
    else:
        defargs["args"] = (x, y)
    
    if cons is None:
        cons = []
    
    rs = m.getRestrictions()
    # Loop over freeParamNames (get order right)
    for i, p in enumerate(m.freeParamNames()):
        if p in rs:
            # There is a restriction for this parameter
            lower, upper = rs[p]
            # A function factory to close over 'j' (i.e., 'i' in the loop)
            # Constraints have to be positive if valid
            def getf(j):
                if (not lower is None) and (not upper is None):
                    # Upper and lower limit
                    f = lambda x:min((x[j]-lower), (upper-x[j]))
                elif (not lower is None):
                    # Only lower limit
                    f = lambda x:x[j]-lower
                elif (not upper is None):
                    # Only upper limit
                    f = lambda x:upper-x[j]
                    #f = lambda *x:
                return f
            cons.append(getf(i))

    defargs["cons"] = cons
    # Otherwise, the constraint function will be called with the same
    # 'args' as the objective function
    defargs["consargs"] = ()
    
    fr = sco.fmin_cobyla(objf, m.freeParamVals(), **defargs)
    m.setFreeParamVals(fr)
    return fr


def fitfmin_powell(m, objf, x, y, yerr=None, **kwargs):
    """
    Use scipy's fmin_powell to fit model.
    """
    # Get keywords and default arguments
    defargs, _ = _introdefarg(sco.fmin_powell, **kwargs)
    
    if not yerr is None:
        defargs["args"] = (x, y, yerr)
    else:
        defargs["args"] = (x, y)
    
    defargs["full_output"] = True
    
    fr = sco.fmin_powell(objf, m.freeParamVals(), **defargs)
    m.setFreeParamVals(fr[0])
    return fr


def fitfmin_l_bfgs_b(m, objf, x, y, yerr=None, **kwargs):
    """
    Use scipy's fmin_l_bfgs_b to fit model.
    
    Parameters
    ----------
    userBounds : dictionary, optional
        Maps parameter name to a tuple of the form (lower, upper) bound.
        If given, the bounds given there will overrule those specified
        in the form of restrictions (should there by any).
    """
    # Get keywords and default arguments
    defargs, _ = _introdefarg(sco.fmin_l_bfgs_b, **kwargs)
    
    if not yerr is None:
        defargs["args"] = (x, y, yerr)
    else:
        defargs["args"] = (x, y)
    
    userBounds = kwargs.get("userBounds") or {}
    
    rs = m.getRestrictions()
    bounds = []
    # Loop over freeParamNames (get order right)
    for i, p in enumerate(m.freeParamNames()):
        if p in userBounds:
            bounds.append(userBounds[p])
        elif p in rs:
            # There is a restriction for this parameter
            bounds.append(rs[p])
        else:
            bounds.append((None, None))
            
    defargs["bounds"] = bounds
    defargs["approx_grad"] = True
    
    fr = sco.fmin_l_bfgs_b(objf, m.freeParamVals(), **defargs)
    m.setFreeParamVals(fr[0])
    return fr


def sampleEMCEE2(m, pargs=(), walkerdimfac=4, scales=None,
              sampleArgs=None, dbfile="chain.emcee", ps=None, emcp=None, toMAP=True):
    """
    MCMC sampling using emcee package.
    
    Sample from the posterior probability distribution using the emcee
    package. 
    
    The emcee sampler can be accessed via the `emceeSampler` attribute,
    which may be used to continue or manipulate sampling.
    
    Parameters
    ----------
    walkerdimfac : int, optional
        Determines the number of walkers. By default, four times the
        number of free parameters is used.
    pargs : tuple
        Content of 'args' argument handed to EnsembleSampler
    scales : dictionary, optional
        The scales argument can be used to control the initial distribution
        of the walkers. By default, all walkers are distributed around the
        location given by the current state of the object, i.e., the current
        parameter values. In each direction, the walker are randomly distributed
        with a Gaussian distribution, whose default standard deviation is 1e-8.
        The scales argument can be used to control the width of Gaussians used
        to distribute the walkers.
    sampleArgs : dictionary, optional
        Keyword controlling the sampling process. Use 'burn' (int) to specify
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
        be carried out and the other arguments should be given as previously
        to continue sampling successfully. 
    emcp : dictionary, optional
        Extra arguments handed to `EnsembleSampler` object.
    toMAP : boolean, optional
        If True (default), the object is set to the maximum posterior probability point sampled.
        Otherwise, it remains in a random state.
    
    Returns
    -------
    (pos, state) : tuple
        Position and state of the sampler. Can be used to continue sampling.
    sampler : emceeSampler
        The emcee object used for sampling.
    """
    
    if not ic.check["emcee"]:
        raise(PE.PyARequiredImport("Could not import the 'emcee' package.",
                                   solution="Please install 'emcee'."))
    
    # Names of free parameters
    fpns = m.freeParamNames()
    # Values of of free parameters
    fpvs = m.freeParamVals()
    # Names and values of free parameter (in no order)
    fps = dict(zip(fpns, fpvs))
    # Number of dimensions
    ndims = len(fpns)
    
    if ndims == 0:
        raise(PE.PyAValError("At least one free parameter is required for sampling.",
                             where="sampleEMCEE2",
                             solution="Use 'thaw' to free same parameters."))
    
    if not dbfile is None:
        if re.match(".*\.emcee$", dbfile) is None:
            PE.warn(PE.PyAValError("The db filename (" + str(dbfile) + ") does not end in .emcee. TraceAnalysis will not recognize it as an emcee trace file.",
                                   solution="Use a filename of the form *.emcee"))
    
    # Number of walkers
    nwalker = ndims * walkerdimfac
    
    if nwalker < ndims * 2:
        raise(PE.PyAValError("The number of walkers must be at least twice the number of free parameters.",
                             where="sampleEMCEE2",
                             solution="Increase the number of walkers."))
    if nwalker % 2 == 1:
        raise(PE.PyAValError("The number of walkers must be even.",
                             where="sampleEMCEE2",
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
    
    if ps is None:
    
        if emcp is None:
            emcp = {}
    
        # Generate the sampler
        emceeSampler = emcee.EnsembleSampler(nwalker, ndims, m.emceeLogPost, args=pargs, **emcp)
    
        if scales is None:
            scales = {}
    
        # Generate starting values
        pos = []
        for _ in smo.range(nwalker):
            pos.append(np.zeros(ndims))
            for i, n in enumerate(fpns):
                # Use 1e-8 unless scale is defined
                s = (scales.get(n) or 1e-8)
                pos[-1][i] = np.random.normal(fps[n], s)
                
        # Default value for state
        state = None
    
        if sampleArgs["burn"] > 0:
            # Run burn-in
            pos, prob, state, blobs = emceeSampler.run_mcmc(pos, sampleArgs["burn"])
            # Reset the chain to remove the burn-in samples.
            emceeSampler.reset()
    
    else:
        # Assign position and state from previous run
        pos, state = ps
    
    if (not sampleArgs["progress"] is None) and ic.check["progressbar"]:
        widgets = ['EMCEE progress: ', progressbar.Percentage(), ' ', progressbar.Bar(marker=progressbar.RotatingMarker()),
                   ' ', progressbar.ETA()]
        pbar = progressbar.ProgressBar(
            widgets=widgets, maxval=sampleArgs["iters"]).start()
    
    n = 0
    for pos, prob, state, blobs in emceeSampler.sample(pos, rstate0=state, iterations=sampleArgs["iters"], thin=1):
        n += 1
        if (not sampleArgs["progress"] is None) and (n % sampleArgs["progress"] == 0):
            if ic.check["progressbar"]:
                pbar.update(n)
            else:
                print("EMCEE: Reached iteration ",
                      n, " of ", sampleArgs["iters"])

    # Save the chain to a file
    if not dbfile is None:
        
        def deblob(b):
            """ Convert blob format (list of lists) into numpy arrays """
            # Shape is (no. of walkers, length of chain)
            lnprior = np.zeros( (len(b[0]), len(b)) )
            lnl = np.zeros( (len(b[0]), len(b)) )
            for i, blob in enumerate(b):
                lnprior[::, i] = np.array([x[0] for x in blob])
                lnl[::, i] = np.array([x[1] for x in blob])
            return lnprior, lnl
        
        lnprior, lnl = deblob(emceeSampler.blobs)
        np.savez_compressed(open(dbfile, 'wb'), chain=emceeSampler.chain, lnpost=emceeSampler.lnprobability, lnprior=lnprior, lnl=lnl,
                            pnames=np.array(fpns, dtype=np.str_))
    
    if toMAP:
        # Set to Maximum-A-Posteriori solution
        indimin = np.argmax(emceeSampler.lnprobability)
        for i, p in enumerate(fpns):
            m[p] = emceeSampler.flatchain[indimin, i]
    
    return (pos, state), emceeSampler
    
    
def _fuf2ModelFactory(odc, rn):
    
    class F2M(MBOEv):
    
        def __init__(self, *args, **kwargs):
            self._odm = odc(*args, **kwargs)
            MBOEv.__init__(self, pars=list(self._odm.parameters()), rootName=rn)
            self.setlogL("1dgauss")
 
        def evaluate(self, *args, **kwargs):
            self._odm.assignValues(self.parameters())
            return self._odm.evaluate(args[0])       
        
    return F2M
       

# Transform OneDFit objects into MBOs
# List gives tuples with fuf2-name, fuf-class, and root name to be assigned
_ml = [("GaussFit1d", _fuf.GaussFit1d, "GF"), \
       ("MultiGaussFit1d", _fuf.MultiGauss1d, "MGF"),
       ("CauchyLorentz1d", _fuf.CauchyLorentz1d, "CauchyLorentz"), \
       ("Voigt1d", _fuf.voigt1d, "Voigt"), \
       ("MultiVoigt1d", _fuf.MultiVoigt1d, "MultiVoigt"), \
       ("SinusFit1d", _fuf.SinusFit1d, "SinusFit"), \
       ("ExpDecayFit1d", _fuf.ExpDecayFit1d, "ExpDecayFit"), \
       ("PolyFit1d", _fuf.PolyFit1d, "PolyFit"), \
       ("ConstantFit1d", _fuf.ConstantFit1d, "ConstantFit"), \
       ("GaussFit2d", _fuf.GaussFit2d, "GaussFit2d"), \
       ("MultiGauss2d", _fuf.MultiGauss2d, "MultiGauss2d"), \
       ("Circle2d", _fuf.Circle2d, "Circle2d")]

# Provide versions/names without the 1d ending
_eml = []
for m in _ml:
    if m[0].endswith("1d"):
        _eml.append( (m[0][0:-2], m[1], m[2]) )
_ml.extend(_eml)
       
# Add models from modelSuite
_ml.extend([("PalLC", _ms.palTrans.PalLC, "PalLC"), \
            ("PalLCKep", _ms.palTrans.PalLCKep, "PalLCKep"), \
            ("MandelAgolLC", _ms.forTrans.MandelAgolLC, "MandelAgolLC"), \
            ("MandelAgolNLLC", _ms.forTrans.MandelAgolNLLC, "MandelAgolNLLC"), \
            ("LimBrightTrans", _ms.LimBrightTrans, "LimBrightTrans"), \
            ("RmcL", _ms.RmcL, "RMcL"), \
            ("RmcLell", _ms.RmcLell, "RmcLell"), \
            ("SinRadVel", _ms.SinRadVel, "SinRadVel"), \
            ("KeplerEllipseModel", _ms.keplerEllipseModel, "KeplerEllipseModel"), \
            ("KeplerRVModel", _ms.KeplerRVModel, "KeplerRVModel"), \
            ("LLGauss", _ms.LLGauss, "LLGauss"), \
            ("VoigtAstroP", _ms.VoigtAstroP, "VoigtAstroP"), \
            ("LyaTransmission", _ms.LyaTransmission, "LyaTransmission"), \
            ("RotBroadProfile", _ms.RotBroadProfile, "RotBroadProfile") 
     ])

for m in _ml:
    setattr(_fbo2module, m[0], _fuf2ModelFactory(m[1], m[2]) )


class Poly2d(MBOEv):
    """
    Bivariate (2d) polynomial model
    
    Parameters
    ----------
    nx, ny : int
        Number of coefficients in first and second variable
    """

    def __init__(self, nx, ny):
        if (nx < 0) or (ny < 0):
            raise(PE.PyAValError("nx and ny must be non-negative integers."))
        self.nx, self.ny = nx, ny
        pars = []
        for x in range(nx+1):
            for y in range(ny+1):
                pars.append("c%d-%d" % (x,y))
        # 'pars' specifies parameter names in the model
        MBOEv.__init__(self, pars=pars, rootName="Poly2d")

    def evaluate(self, x, *args):
        """
        Evaluate model at points (x[::,0], x[::,1])
        """
        c = np.zeros( (self.nx+1, self.ny+1) )
        for i in range(self.nx+1):
            for j in range(self.ny+1):
                c[i,j] = self["c%d-%d" % (i,j)]
        return np.polynomial.polynomial.polyval2d(x[::,0], x[::,1], c)

# 
# # class GF1d(MBO2):
# #     
# #     def __init__(self):
# #         self._gf = GaussFit1d()
# #         MBO2.__init__(self, pars=list(self._gf.parameters()), rootName=GF1d)
# #         self.setLogL("1dgauss")
# #     
# #     def evaluate(self, t):
# #         self._gf.assignValues(self.parameters())
# #         return self._gf.evaluate(t)
# #     
# # 
# #     
# class Poly2(MBO2):
# 
#     def __init__(self):
#         MBO2.__init__(self, pars=["c0", "c1", "c2"], rootName="Poly2")
#         self.setlogL("1dgauss")
#     
#     def evaluate(self, x, **kwargs):
#         s = self._imap
#         return s["c0"] + s["c1"]*x + s["c2"]*x**2
#     
# 
# class CeleriteModel(MBO2):
#     
#     def _nc(self, s):
#         """
#         Adjust parameter names
#         
#         Parameters
#         ----------
#         s : string
#             Parameter name from celerite
#         
#         Returns
#         -------
#         name : string
#             Adjusted name (here, 'kernel: ' removed)
#         """
#         s = s.replace("kernel:", "")
#         return s
#     
#     def __init__(self, gp, nc=None):
#         
#         # Name conversion
#         if nc is None:
#             # Use default naming convention
#             nc = self._nc
#         
#         # Save reference to GP
#         self.gp = gp
#         # Save parameter names (and the order)
#         self._pns = [nc(n) for n in gp.get_parameter_names(include_frozen=True)]
#         
#         MBO2.__init__(self, self._pns, rootName="celmo", logl=None)
#         
#         vals = gp.get_parameter_vector(include_frozen=True)
#         for i, n in enumerate(self._pns):
#             self[n] = vals[i]
#         
#         for i, b in enumerate(self.gp.get_parameter_bounds(include_frozen=True)):
#             self.setRestriction({self._pns[i]:b})
#         
#     def _updateGPParams(self):
#         """ Assign current model parameters """
#         pv = np.zeros(len(self._pns))
#         for i, _ in enumerate(self.gp.get_parameter_names()):
#             pv[i] = self[self._pns[i]]
#         self.gp.set_parameter_vector(pv, include_frozen=True)
#     
#     def logL(self, *args, **kwargs):
#         """ Likelihood function for celerite GP """
#         self._updateGPParams()
#         ll = self.gp.log_likelihood(args[1])
#         return ll
#     
#     def evaluate(self, x):
#         return None
#     
#     def predict(self, px, y):
#         """
#         Get predicted mean and variance of GP
#         
#         Parameters
#         ----------
#         px : array
#             Position at which to evaluate the prediction
#         y : array
#             The data
#         
#         Returns
#         -------
#         mean, variance, std : arrays
#             Prediction
#         """
#         self._updateGPParams()
#         pred_mean, pred_var = self.gp.predict(y, px, return_var=True)
#         return pred_mean, pred_var, np.sqrt(pred_var)
# 
#     
