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
    
    def __init__(self, value=0.0):
        self.value = value
        self.free = False
        self.affects = []
        self.givenBy = None


class PyAPS(object):
    
    
    def addParam(self, name, value=0.0):
        """
        """
        if name in self._pmap:
            raise(PE.PyAValError("Name '" + str(name) + "' already exists!", \
                                 where="addPyAPS::Param"))
        self._pmap[name] = PyAPa(value)
    
    def __init__(self, *args, **kwargs):
        self._pmap = collections.OrderedDict()
        for n in args:
            self.addParam(n, 0.0)
        for n, v in kwargs.items():
            self.addParam(n, v)
            

class FBO2(object):
    
    pass



