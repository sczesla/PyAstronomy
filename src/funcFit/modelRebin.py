from __future__ import print_function, division
import numpy
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo

def turnIntoRebin(CO):
  """
    Turn a "normal" fitting object into rebinning fitting object.
  
    This function accepts a class object representing a model and
    returns another class object extended by the rebinning functionality.
  
    Parameters
    ----------
    CO : A class object
        The class object describing the model to use
        rebinning.
        
    Returns
    -------
    Rebinned model : Fitting object
        Another class object extended by the rebinning functionality.
  """
  
  # Check for evaluate
  if not hasattr(CO, "evaluate"):
    raise(PE.PyANotImplemented("The class object has no evaluate method. Is it really a model?", \
                               where="turnIntoRebin"))
  # Check for other properties and methods
  for p in ["setRebinArray_Ndt", "rebinTimes", "rebinIdent"]:
    if hasattr(CO, p):
      PE.warn(PE.PyANameClash("Class object already has attribute: `"+p+"`. Potentially harmful..", \
              solution="Check class object, maybe rename attribute", where="turnIntoRebin" ))
  
  class _ModelRebin(CO):
    """
      Base class providing rebinning functionality.
      
      The model is evaluated at more points than actually needed.
      Several points are than averaged to obtain a "binned" model,
      which can, for example, account for finite integration times
      in observations.
      
      Attributes
      ----------
      rebinTimes : array
          Defined the abscissa values at which
          to evaluate the model to be rebinned afterwards.
      rebinIdent : dict
          A dictionary associating bin number (in the unbinned
          model) with a list holding the bins in `rebinTimes`,
          which are to be averaged to obtain the binned model.
    """
  
    def __init__(self, *args, **kwargs):
      CO.__init__(self, *args, **kwargs)
      self.rebinTimes = None
      self.rebinIdent = None

    def setRebinArray_Ndt(self, time, N, dt):
      """  
        Defines the overbinning parameters (`rebinTimes`, `rebinIdent`).
        
        It is assumed that the time points given in the `time`
        array refer to the center of the time bins and every bin has
        length `dt`. The bins are then subdivided into `N` subintervals;
        the center of each such subinterval becomes
        a point in the overbinned time axis (`rebinTimes`).
  
        Parameters
        ----------
        time : array
            The time axis of the "observed" (not overbinned)
            transit light-curve.
        N : int
            The number of point into which to subdivide each time
            bin of length `dt`.
        dt : float
            The length of each time bin (on the original not
            oversampled time axis).
      """
      
      self.rebinTimes = numpy.zeros(time.size * N)
      self.rebinIdent = {}
      for i in smo.range(time.size):
        self.rebinTimes[i*N:(i+1)*N] = \
              (time[i] - dt/2.0) + (numpy.arange(N)*dt)/float(N) + dt/float(N)/2.0
        self.rebinIdent[i] = list(range(i*N,(i+1)*N))
  

    def evaluate(self, x):
      """
        Calculate the model.

        Parameters
        ----------
        x : array
            The abscissa values.
        
        Returns
        -------
        model : array,
            The binned model.

        Notes
        -----
        This function calculates the model at those time points
        specified by the `rebinTimes` property and saves the result in the
        class property `unbinnedModel`. Then it bins
        according to the definitions in `rebinIdent` and save the resulting model
        in the `binnedModel` property.
      """
      if (self.rebinTimes is None) or (self.rebinIdent is None):
        raise(PE.PyAValError("Rebinning parameters (properties rebinTimes and/or rebinIdent) not appropriately defined.",
                             solution=["Use setRebinArray_Ndt method.", "Define the properties by accessing them directly."]))
      # Calculate the model using current parameters at the time points
      # defined in the `rebinTimes` array
      self.unbinnedModel = CO.evaluate(self, self.rebinTimes)
      # Build up rebinned model
      self.binnedModel = numpy.zeros(x.size)
      for i in smo.range(x.size):
        self.binnedModel[i] = numpy.mean(self.unbinnedModel[self.rebinIdent[i]])
      # Return the resulting model
      return self.binnedModel

  return _ModelRebin


# The code below is only relevant for Sphinx documentation

class _DummyModel:
  def evaluate(self):
    pass
  
_ModelRebinDocu = turnIntoRebin(_DummyModel)
_ModelRebinDocu.__name__="_ModelRebinDocu"