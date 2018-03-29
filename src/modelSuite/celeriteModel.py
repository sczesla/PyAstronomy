from __future__ import print_function, division
from PyAstronomy import funcFit as fuf
import new
from PyAstronomy.pyaC import pyaErrors as PE
import numpy as np

class CeleriteModel(fuf.OneDFit):
    
    def _nc(self, s):
        """
        Adjust parameter names
        
        Parameters
        ----------
        s : string
            Parameter name from celerite
        
        Returns
        -------
        name : string
            Adjusted name (here, 'kernel: ' removed)
        """
        s = s.replace("kernel:", "")
        return s
    
    def __init__(self, gp, nc=None):
        
        if nc is None:
            # Use default naming convention
            nc = self._nc
        
        # Save reference to GP
        self.gp = gp
        # Save parameter names (and the order)
        self._pns = [nc(n) for n in gp.get_parameter_names()]
        
        fuf.OneDFit.__init__(self, self._pns, rootName="celmo")
        
        @fuf.MiniFunc(self)
        def mini(self, P):
            # Negative likelihood
            return -self.objFctLL()
        
        @fuf.MiniFunc(self)
        def miniPlus(self, P):
            # Likelihood
            return self.objFctLL()
        
        # Preserve the initial fit function
        self._initalfit = self.fit
        
        # Monkey patch to replace the default miniFunc by likelihood
        def newfit(_, *args, **kwargs):
            if (not "miniFunc" in kwargs) or (kwargs["miniFunc"] is None):
                kwargs["miniFunc"] = mini
            return self._initalfit(*args, **kwargs)
        self.fit = new.instancemethod(newfit, self, CeleriteModel)
        
        # Preserve initial fitEMCEE
        self._initialfitEMCEE = self.fitEMCEE
        
        # Monkey patch to replace the default likelihood in fitEMCEE
        def newfitEMCEE(_, *args, **kwargs):
            if (not "llfunc" in kwargs) or (kwargs["llfunc"] is None):
                kwargs["llfunc"] = miniPlus
            return self._initialfitEMCEE(*args, **kwargs)
        self.fitEMCEE = new.instancemethod(newfitEMCEE, self, CeleriteModel)
        
    def objFctLL(self):
        ll = self.gp.log_likelihood(self._fufDS.y)
        return ll
    
    def evaluate(self, x):
        pvs = [self[p] for p in self._pns]
        self.gp.set_parameter_vector(pvs)
        
        if (not hasattr(self, "_gpcomx")) or (not self._gpcomx is self._fufDS.x):
            # Call to compute only if necessary
            self.gp.compute(self._fufDS.x)
            self._gpcomx = self._fufDS.x
    
    