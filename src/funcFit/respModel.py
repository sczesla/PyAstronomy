from __future__ import print_function, division
from PyAstronomy.pyaC import pyaErrors as PE
import copy
from .onedfit import OneDFit
import six

class ConvolutionModel(OneDFit):
    """
      Generate a convolved funcFit model with full fitting functionality.
      
      In many cases, it can be useful to modify a model by a convolution with
      some function. Examples if this comprise a spectral response function or
      instrumental broadening in optical high-resolution spectra.  
      
      This model 'wraps around' the original model, uses its output, and convolves
      it with some kernel. 
      
      Parameters
      ----------
      pmodel - OneDFit instance
          FuncFit model representing the model to be folded with the response.
      resp : class instance
          A class instance, which must provide a 'convolve' method, which
          must be of the form convolve(x, y) and the result must satisfy
          the condition len(convolve(x, y)) == len(x). The convolve method
          is called whenever the model is evaluated. In particular, the first
          parameter is the x-values it which to evaluate and the second parameter
          is the result of the evaluation of 'pmodel'. If the class instance is
          also of type 'OneDFit', i.e., a funcFit object, it can have free
          parameters adapted in the fitting process.
        
      
      Example
      -------
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
        # check if response model has convolve method
        if not hasattr(resp, "convolve"):
            raise(PE.PyAAlgorithmFailure("Response model has no method 'convolve'.", \
                              where="ConvolutionModel::__init__"))
            
        if not isinstance(pmodel, OneDFit):
            raise(PE.PyAAlgorithmFailure("Physical model must be derived from OneDFit", \
                              where="ConvolutionModel::__init__"))
            
        # Copy the input models
        self.response = copy.deepcopy(resp)
        self.pmodel = copy.deepcopy(pmodel)
        
        if isinstance(self.response, OneDFit):
            # The response is also a funcFit model and, therefore, could
            # already contribute fitting parameters
            # 
            # Change root name to "Response"
            self.response.setRootName("Response")
            # Ensure unique parameter names    
            self._combineRemapping(self.response, self.pmodel)
            # List of parameters for the this (the convolved) model
            parameters = list(self.response.parameters()) + list(self.pmodel.parameters())
        else:
            # List of parameters for the this (the convolved) model
            parameters = list(self.pmodel.parameters())
        
        # Set up funcFit model and assign
        OneDFit.__init__(self, parameters)
        
        # Define root name
        if self.pmodel.naming._root == "":
          # Use default 
          rn = "model"
        else:
          rn = self.pmodel.naming._root
        
        self.setRootName("Conv-" + rn)
        
    def evaluate(self, x):
        
        if isinstance(self.response, OneDFit):
            # Assign parameters if OneDFit instance
            for p in six.iterkeys(self.response.parameters()):
                self.response[p] = self[p]
        for p in six.iterkeys(self.pmodel.parameters()):
          self.pmodel[p] = self[p]
               
        return self.response.convolve(x, self.pmodel.evaluate(x))