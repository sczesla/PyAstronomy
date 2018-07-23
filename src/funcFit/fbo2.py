# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
import re
import copy
import types
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import pyaC
from .fufDS import FufDS
from .extFitter import NelderMead
import six
import six.moves as smo
import bidict

from PyAstronomy.funcFit import _pymcImport, _scoImport, ic

if _pymcImport:
    import pymc
if _scoImport:
    import scipy.optimize as sco
if ic.check["emcee"]:
    import emcee
if ic.check["progressbar"]:
    import progressbar


class PyAPa(object):
    
    def setVal(self, v):
        self._value = v
    
    def getVal(self):
        return self._value
    
    value = property(getVal, setVal)
    
    def __init__(self, value=0.0):
        self._value = value
        self.free = False
        self.affects = set()
        self.affectsPriors = set()
        self.relation = None
        self.dependsOn = None
        self.freeable = True
        
    def updateRelation(self):
        """
        Update value according to relation
        """
        if not self.relation is None:
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
    
    def __init__(self, ivs, func):
        """
        Parameters
        ----------
        ivs : list of Parameter instances
            The independent parameters
        func : Callable
            Function returning log(prior) depending on the value of the
            independent parameters
        """
        self.ivs = ivs
        self.func = func
    
    def update(self):
        return self.func(*[p.value for p in self.ivs])
    
      
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
        """
        """
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
        if not r is None:
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
    
    def addPrior(self, pyap):
        for p in pyap.ivs:
            p.affectsPriors.update(pyap)
    
    def __init__(self, *args):
        self.bss = []
        for s in args:
            self.addBPS(s)
        self._updatepmap()
      
      
    

class MBO2(object):
    """
    Model Base Object
    
    Concepts
    --------
    
    Parameter names: Internal vs. external
        
    
    
    """
    
    def __init__(self, pars, rootName=""):
        self.pars = (PyABaSoS(PyABPS(pars, rootName)))
        self._imap = self.pars.copy()
        
        self.rootName = rootName
        self.leftCompo = None
        self.rightCompo = None
     
    def _combineMBOs(self, right):
        """
        Combine two MBO2s into a new one

        Parameters
        ----------
        right : MBO2 instance
            Right side of the operation
        """
        r = MBO2([], rootName="combined")
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

    def parameterSummary(self):
        for k, v in self.pars.pmap.items():
            print(k, v.value)
        
    def evaluate(self, *args, **kwargs):
        pass
    
    def logL(self, *args, **kwargs):
        pass
    
    def relate(self, dv, idv, func=None):
        self.pars.relate(dv, idv, func)
    
    def _defLogL(self, *args, **kwargs):
        if len(args) == 2:
            x, y = args[0], args[1]
            m = self.evaluate(x)
            return np.sum((m - y)**2)
        elif len(args) == 3:
            x, y, yerr = args[0], args[1], args[2]
            m = self.evaluate(x)
            return np.sum((m-y)**2/yerr**2)
        else:
            raise(PE.PyAValError("Invalid call to _defLogL. No. of arguments: " + str(len(args))))
    
    def __getitem__(self, n):
        return self.pars[n]
    
    def __setitem__(self, n, v):
        self.pars[n] = v
        
    
    
class Poly2(MBO2):
    
    def __init__(self):
        MBO2.__init__(self, ["c0", "c1", "c2"], rootName="Poly2")
    
    def evaluate(self, *args, **kwargs):
        s = self._imap
        x = args[0]
        return s["c0"] + s["c1"]*x + s["c2"]*x**2
    
    def logL(self, *args, **kwargs):
        return super()._defLogL(*args)
    



