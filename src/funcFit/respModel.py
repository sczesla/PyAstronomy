#import PyAstronomy.funcFit as fuf
import matplotlib.pyplot as plt
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
import copy
from .onedfit import OneDFit

class ConvolutionModel(OneDFit):
    """
      Provides functionality to convolve a funcFit with a response
      while retaining the full fitting capabilities of funcFit models.
      
      Example:
      --------
      import PyAstronomy.funcFit as fuf
      import matplotlib.pyplot as plt
      import numpy as np

      # Create an response based on a Gaussian profile      
      class Broad(fuf.GaussFit1d):  
        def convolve(self, x, y):
            tmp = self.evaluate(x)
            tmp/=np.sum(tmp) * self["A"]
            #out = np.convolve(tmp, y)[len(x)/2:len(x)+len(x)/2]
            return np.convolve(tmp, y, 'same')

        resp = Broad()
        
        # The full OneDFit functionality remains intact:
        pm = fuf.GaussFit1d() + fuf.PolyFit1d(1)
        m = ConvolutionModel(resp, pm)
        
        # Set some parameters
        m["A_Gaussian"] = 1.0
        m["sig_Gaussian"] = 1.0
        m["sig_Response"] = 1.0
        m["A_Response"] = 1.0
        m.parameterSummary()
        
        # Plot the results
        x = np.linspace(-5, 5, 101)
        yy = m.evaluate(x)
        plt.plot(x, yy, label="Convolved")
        pm["A_Gaussian"] = 1.0
        pm["sig_Gaussian"] = 1.0
        plt.plot(x, pm.evaluate(x), label="input")
        plt.legend()
        plt.show()
      
    """
    def __init__(self, resp, pmodel):
        """
        
          Parameters
          ----------
          resp - class instance; must provide the method 'convolve', which
                 must be of the form convolve(x, y) with len(convolve(x, y))==len(x).
                 Represents the instrument response.
                 Class can be derived from OneDFit so that its parameters
                 can be varied during the fit simultaneously with the 
                 parameters of the physical model
          pmodel - OneDFit instance
                   Physical model
        """
        
        # check if response model has convolve method
        if not hasattr(resp, "convolve"):
            raise(PE.pyaErrors.pyaOtherErrors.PyAAlgorithmFailure("Response model has no method 'convolve'.", \
                              where="ConvolutionModel::__init__"))
        if not isinstance(pmodel, OneDFit):
            raise(PE.pyaErrors.pyaOtherErrors.PyAAlgorithmFailure("Physical model must be derived from OneDFit", \
                              where="ConvolutionModel::__init__"))
        # Sollte man hier ein copy.deepcopy verwenden?
        self.left = copy.deepcopy(resp)
        self.right = copy.deepcopy(pmodel)
        if isinstance(self.left, OneDFit):
            self.left.setRootName("Response")
            self._combineRemapping(self.left, self.right)
            parameters = self.left.parameters().keys() + self.right.parameters().keys()
        else:
            parameters = self.right.parameters().keys()
        OneDFit.__init__(self, parameters)
        self.convolve = self.left.convolve
        self.setRootName("ConvolutionModel")
    def evaluate(self, x):
        if isinstance(self.left, OneDFit):
          for p in self.left.parameters().keys():
            self.left[p] = self[p]
        for p in self.right.parameters().keys():
          self.right[p] = self[p]          
        return self.convolve(x, self.right.evaluate(x))
