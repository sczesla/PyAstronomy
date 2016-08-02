import PyAstronomy.funcFit as fuf
import matplotlib.pyplot as plt
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
import copy

class ConvolutionModel(fuf.OneDFit):
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
        if not isinstance(pmodel, fuf.OneDFit):
            raise(PE.pyaErrors.pyaOtherErrors.PyAAlgorithmFailure("Physical model must be derived from OneDFit", \
                              where="ConvolutionModel::__init__"))
        # Sollte man hier ein copy.deepcopy verwenden?
        self.left = copy.deepcopy(resp)
        self.right = copy.deepcopy(pmodel)
        if isinstance(self.left, fuf.OneDFit):
            self.left.setRootName("Response")
            self.combineRemapping(self.left, self.right)
            parameters = self.left.parameters().keys() + self.right.parameters().keys()
        else:
            parameters = self.right.parameters().keys()
        fuf.OneDFit.__init__(self, parameters)
        self.convolve = self.left.convolve
        self.setRootName("ConvolutionModel")
    def evaluate(self, x):
        if isinstance(self.left, fuf.OneDFit):
          for p in self.left.parameters().keys():
            self.left[p] = self[p]
        for p in self.right.parameters().keys():
          self.right[p] = self[p]          
        return self.convolve(x, self.right.evaluate(x))

## ===============================================
## Example 1
#class Response1():
    #def __init__(self):
        #pass
    #def convolve(self, x, y):
        #return y*1.2
    
#resp = Response1()
#r = fuf.GaussFit1d()
#r["sig"] = 1.0
#r["A"] = 1.0

#m = ConvolutionModel(resp, r)

#m["sig"] = 1.0
#m["A"] = 1.0

#x = np.linspace(-5, 5, 100)
#yy = m.evaluate(x)
#plt.plot(x, yy)
#plt.plot(x, r.evaluate(x), label="right")
#plt.legend()
#plt.show()

## ===============================================
## Example 2
    
#class Response2(fuf.OneDFit):
    #def __init__(self):
        #fuf.OneDFit.__init__(self,["scaling"])
    #def convolve(self, x, y):
        #tmp = self.evaluate(x)
        #return np.convolve(tmp, y, mode='same')
    #def evaluate(self, x):
        #return self["scaling"]

#resp = Response2()
#pm = fuf.GaussFit1d()

#m = ConvolutionModel(resp, pm)

#m["sig_Gaussian"] = 1.0
#m["A_Gaussian"] = 1.0
#m["scaling_Response"] = 1.2
#m.parameterSummary()

#x = np.linspace(-5, 5, 100)
#yy = m.evaluate(x)
#plt.plot(x, yy)
#pm["A"] = 1.0
#pm["sig"] = 1.0
#plt.plot(x, pm.evaluate(x), label="right")
#plt.legend()
#plt.show()

      
      
## ===============================================
## Example 3
    
#class Broad(fuf.GaussFit1d):  
  #def convolve(self, x, y):
    #tmp = self.evaluate(x)
    #tmp/=np.sum(tmp) * self["A"]
    #return np.convolve(tmp, y, mode='same')

#resp = Broad()

#pm = fuf.GaussFit1d()
#pm["A"] = 1.0
#pm["sig"] = 1.0

#m = ConvolutionModel(resp, pm)

#m["sig_Gaussian"] = 2.0
#m["A_Gaussian"] = 1.0
#m["sig_Response"] = 1.0
#m["A_Response"] = 1.0
#m.parameterSummary()

#x = np.linspace(-5, 5, 100)
#yy = m.evaluate(x)
#plt.plot(x, yy)
#plt.plot(x, pm.evaluate(x), label="right")
#plt.legend()
#plt.show()


      
## ===============================================
## Example 3
    
#class Broad(fuf.GaussFit1d):  
  #def convolve(self, x, y):
    #tmp = self.evaluate(x)
    #tmp/=np.sum(tmp) * self["A"]
    ##out = np.convolve(tmp, y)[len(x)/2:len(x)+len(x)/2]
    #return np.convolve(tmp, y, 'same')

#resp = Broad()

#pm = fuf.GaussFit1d() + fuf.PolyFit1d(1)

#m = ConvolutionModel(resp, pm)
#m["A_Gaussian"] = 1.0
#m["sig_Gaussian"] = 1.0
#m["sig_Response"] = 1.0
#m["A_Response"] = 1.0

#m.parameterSummary()
#x = np.linspace(-5, 5, 101)
#yy = m.evaluate(x)
#plt.plot(x, yy)
#pm["A_Gaussian"] = 1.0
#pm["sig_Gaussian"] = 1.0
#plt.plot(x, pm.evaluate(x), label="input")
#pm["A_Gaussian"] = 1.0
#pm["sig_Gaussian"] = np.sqrt(2)*1.0
#plt.plot(x, pm.evaluate(x), label="output")
#plt.legend()
#plt.show()
