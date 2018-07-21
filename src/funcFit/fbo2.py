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
import collections
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
        self.relation = None
        self.dependsOn = None
        self.freeable = True
        
    def updateRelation(self):
        """
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
    
        

class PyAPS(object):
    
    def addParam(self, name, value=0.0):
        """
        """
        if name in self.pmap:
            raise(PE.PyAValError("Name '" + str(name) + "' already exists!", \
                                 where="addPyAPS::Param"))
        self.pmap[name] = PyAPa(value)
    
    def _checkParam(self, n):
        if not n in self.pmap:
            raise(PE.PyAValError("No such parameter: " + str(n),
                                 solution="Use one of: " + ', '.join(list(self.pmap))))
    
    def __getitem__(self, n):
        self._checkParam(n)
        return self.pmap[n]
    
    def __setitem__(self, n, v):
        self._checkParam(n)
        self.pmap[n].value = v
        if len(self.pmap[n].affects) > 0:
            for a in self.pmap[n].affects:
                a.updateRelation()

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
        # Manage independent variables
        for n in idv:
            self.pmap[n].affects.update([self.pmap[dv]])
        # Trigger update of the value of the related parameter
        self.pmap[dv].updateRelation()
    
    def __init__(self, *args, **kwargs):
        self.pmap = bidict.OrderedBidict()
        for n in args:
            self.addParam(n, 0.0)
        for n, v in kwargs.items():
            self.addParam(n, v)
            

class FBO2(object):
    
    def __init__(self, pars, rootName=None):
        self.pars = PyAPS(*pars)
        self.rootName = rootName
        self.leftCompo = self
        self.rightCompo = None
     
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
        
    
    
class Poly2(FBO2):
    
    def __init__(self):
        super().__init__(["c0", "c1", "c2"], rootName="Poly2")
    
    def evaluate(self, *args, **kwargs):
        x = args[0]
        return self["c0"] + self["c1"]*x + self["c2"]*x**2
    
    def logL(self, *args, **kwargs):
        return super()._defLogL(*args)
    



