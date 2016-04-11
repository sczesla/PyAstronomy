from __future__ import print_function, division
from PyAstronomy.pyaC import pyaErrors as PE
import os
import six

class AtomicNo:
  
  def __init__(self):
    """
    """
    
    path = os.path.dirname(__file__)
    
    self._data = {}
    self._elNames = {}
    for l in open(os.path.join(path,"atomicNo.dat")):
      l = l.rstrip("\n")
      s = l.split()
      self._data[int(s[2])] = s[1]
      self._elNames[int(s[2])] = s[0]
    # Add reverse mapping to _data dict, i.e., element symbol->atomic number
    self._data.update( dict((self._data[k], k) for k in self._data) )
    
  def getElSymbol(self, atn):
    """
      Convert atomic number into elemental symbol.
      
      Parameters
      ----------
      atn : int
          Atomic number
          
      Returns
      -------
      Symbol : string
          Elemental symbol
    """
    try:
      atn = int(atn)
    except ValueError:
      raise(PE.PyAValError("`atn` needs to be an integer. You gave: " + str(atn)))
    
    if not atn in self._data:
      raise(PE.PyAValError("No such atomic number in the data base: " + str(atn)))
    
    return self._data[atn]
  
  def getAtomicNo(self, sym):
    """
      Convert elemental symbol into atomic number.
      
      Parameters
      ----------
      sym : string
          Elemental symbol
          
      Returns
      -------
      Atomic number : int
          Atomic number
    """
    if not isinstance(sym, six.string_types):
      raise(PE.PyAValError("`sym` must be a string, e.g., 'Fe'."))
    
    if not sym in self._data:
      raise(PE.PyAValError("No such symbol in the data base: " + str(sym)))
    
    return self._data[sym]
  
  def showAll(self):
    """
      Print atomic number, elemental symbol, and element name.
    """
    nos = [k for k in self._data.keys() if type(k) == int]
    for k in sorted(nos):
      if not type(k) == int:
        continue
      print("%3d | %3s | " % (k, self._data[k]) + self._elNames[k])
  
  def getElementName(self, atn):
    """
      Convert atomic number into elemental name.
      
      Parameters
      ----------
      atn : int
          Atomic number
          
      Returns
      -------
      Name : string
          Name of element
    """
    try:
      atn = int(atn)
    except ValueError:
      raise(PE.PyAValError("`atn` needs to be an integer. You gave: " + str(atn)))
    
    if not atn in self._data:
      raise(PE.PyAValError("No such atomic number in the data base: " + str(atn)))
    
    return self._elNames[atn]