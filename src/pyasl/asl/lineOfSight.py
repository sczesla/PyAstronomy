from __future__ import print_function, division
import re
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
import six

class LineOfSight:
  """
    Define the "line of sight".
    
    Constructs a line of sight, i.e., a unit vector in the
    given direction.
    
    Parameters
    ----------
    los : string, {'+x','-x','+y','-y','+z','-z'}, or array
        The line of sight given in string format.
    
    Attributes
    ----------
    
  """
  
  def __init__(self, los='-z'):

    if isinstance(los, six.string_types):
      # The input is a string. Check the content.
      r = re.match("^([-+])([xyz])$", los)
      if r is None:
        raise(PE.PyAValError("A line of sight encoded by " + los + " cannot be constructed."))
      if r.group(2) == "x":
        self.los = np.array([1.0,0.0,0.0])
      elif r.group(2) == "y":
        self.los = np.array([0.0,1.0,0.0])
      elif r.group(2) == "z":
        self.los = np.array([0.0,0.0,1.0])
      
      if r.group(1) == "-":
        self.los *= -1.0
      
      self.losName = los
      return
    
    # Try to convert input to numpy array
    convlos = np.array(los)
    # Check length
    if len(convlos) != 3:
      raise(PE.PyAValError("A line of sight argument does not have the right length of 3."))
    # Check that it is a unit vector
    if np.abs(np.sum(convlos**2) - 1.0) > 1e-8:
      raise(PE.PyAValError("A line of sight argument is not a unit vector."))
    
    self.los = convlos
    self.losName = str(convlos)
      
      
  
  def __str__(self):
    """
      String output of the class.
    """
    return "Line of sight (" + self.losName + "): " + str(self.los)