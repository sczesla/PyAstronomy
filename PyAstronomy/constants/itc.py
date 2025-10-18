from __future__ import print_function, division
from PyAstronomy.pyaC import ImportCheck as _ImportCheck
from PyAstronomy.pyaC import pyaErrors as _PE
import PyAstronomy as _PyA
from PyAstronomy.constants import _systems
import os
from six.moves import configparser as _ConfigParser
import numpy
import six

_ic = _ImportCheck(["quantities"], ["quantities"])

import quantities as _pq


def _displaySummary(scope):
  """
    Display unit system and constants.
    
    Parameters
    ----------
    scope : module or class
        The scope from which to get the data to
        be displayed.
  """
  print("Current unit system: ", getattr(scope, "getSystem")())
  print()
  
  cmod = scope
  inventory = getattr(scope, "inventory")
    
  maxlenSym = 0
  maxlenDescr = 0
  maxlenVal = 0
  for k in six.iterkeys(inventory):
    maxlenSym = max(maxlenSym, len(k))
    maxlenDescr = max(maxlenDescr, len(inventory[k]["descr"]))
    maxlenVal = max(maxlenVal, len("%10.5e %s" % (getattr(cmod, "f_"+k).magnitude, getattr(cmod, "f_"+k).dimensionality)))
  maxlenSym = max(maxlenSym, len("Symbol"))
  maxlenDescr = max(maxlenDescr, len("Description"))
  maxlenVal = max(maxlenVal, len("Value"))
  
  totlen = maxlenSym + maxlenDescr + maxlenVal + 6
  
  print(("%" + str(maxlenSym) + "s | %" + str(maxlenDescr) + "s | %" + \
         str(maxlenVal) + "s") % ("Symbol", "Description", "Value"))
  print("-" * maxlenSym + "---" + '-' * maxlenDescr + "---" + '-' * maxlenVal)

  for k in sorted(inventory.keys()):
    print((("%" + str(maxlenSym) + "s | ") % k) + \
          ("%" + str(maxlenDescr) + "s") % inventory[k]["descr"] + " | " + \
          "%10.5e %s" % (getattr(cmod, "f_"+k).magnitude, getattr(cmod, "f_"+k).dimensionality))
  print("-" * totlen)


class PyAConstants:
  """
    Class scope for constants.
    
    Parameters
    ----------
    fn : string, optional
        If given, constants will be read from the specified file.
    unitSystem : string, optional, {"SI", "cgs"}
        The system of units to be used. Default is "cgs".
  """
  
  def _inventoryToScope(self):
    """
      Translate information in `inventory` into class attributes.
    """
    for k, u in six.iteritems(self.inventory):
      val = _pq.Quantity(numpy.float64(u["valueSI"]), u["units"]["SI"])
      value = val.rescale(u["units"][self._unitSystem])
      setattr(self, u["symbol"], value.magnitude)
      setattr(self, "f_" + u["symbol"], value)
      err = _pq.Quantity(numpy.float64(u["errSI"]), u["units"]["SI"])
      err = err.rescale(u["units"][self._unitSystem])
      setattr(self, "f_" + u["symbol"] + "_err", err)

  def summary(self):
    """
      Print a summary of available constants to screen.
    """
    _displaySummary(self)

  def setSystem(self, system):
    """
      Set the unit system.
      
      Parameters
      ----------
      system : string, {SI, cgs}
          The unit system in use.
    """
    if not system in _systems:
      raise(_PE.PyAValError("No such system: '" + system + "'.", \
                            solution="Choose one of: "+str(_systems)))
    self._unitSystem = system
    self._inventoryToScope()

  def inUnitsOf(self, const, units):
    """
      Get a constant and apply unit conversion.
      
      Parameters
      ----------
      const : string
          Name of constant.
      units : string
          The units to be converted into.
      
      Returns
      -------
      Constant : Quantity object (`quantities` package)
          The value and unit of the constant.
    """
    if not hasattr(self, "f_" + const):
      raise(_PE.PyAValError("No such constant: '" + str(const) + "'"))
    val = getattr(self, "f_" + const)
    return val.rescale(units)

  def cleanUp(self):
    """
      Clean up loaded constants.
      
      Empties the `inventory` and deletes all
      related module attributes.
    """
    for k, u in six.iteritems(self.inventory):
      delattr(self, u["symbol"])
      delattr(self, "f_" + u["symbol"])
      delattr(self, "f_" + u["symbol"] + "_err")
    self.inventory = {}

  def constantDetails(self, const):
    """
      Print all information for a specific constant.
      
      Parameters
      ----------
      const : string
          The constant.
      
      Returns
      -------
      info : dictionary
          A dictionary containing all information about the
          specified constant.
    """
    if not const in self.inventory:
      raise(_PE.PyAValError("No such constant: '" + str(const) + "'."))
    print("Constant              : '" + const + "'")
    print("Description           : " + self.inventory[const]["descr"])
    print("Current value and unit: " + str(getattr(self, "f_"+const)))
    print("Uncertainty           : " + str(getattr(self, "f_"+const+"_err")))
    print("Source                : " + self.inventory[const]["source"])
    return self.inventory[const]

  def getSystem(self):
    """
      Get current unit system.
      
      Returns
      -------
      System : string, {SI, erg}
          The unit system currently in use.
    """
    return self._unitSystem

  def load(self, fn, constClobber=False):
    """
      Load constants file.
      
      Parameters
      ----------
      fn : string
          The file name.
      constClobber : boolean, optional
          If True, existing constants will be overwritten.
    """
    config = _ConfigParser.RawConfigParser()
    config.read(fn)
    for s in config.sections():
      nc = None
      try:
        nc = {}
        nc["symbol"]  = config.get(s, "symbol")
        nc["descr"]   = config.get(s, "descr")
        nc["valueSI"] = config.get(s, "valueSI")
        nc["errSI"]   = config.get(s, "errSI")
        nc["units"] = {}
        for u in _systems:
          nc["units"][u] = config.get(s, "unit"+u)
        nc["source"] = config.get(s, "source")
      except _ConfigParser.NoOptionError as e:
        _PE.warn(_PE.PyAValError("The definition of the constant defined in section '" + s + "' of file '" + \
                                 fn + "' is incomplete. Ignoring this constant. Caught the error:\n" + \
                                 str(e)))
        nc = None
      if nc is not None:
        if (nc["symbol"] in self.inventory) and (not constClobber):
          raise(_PE.PyANameClash("A constant with the symbol '" + nc["symbol"] + "' is already defined.", \
                                 solution=["Choose another symbol.", "Set the `constClobber` flag to True to overwrite existing definitions of constants."])) 
        self.inventory[nc["symbol"]] = nc
      self._inventoryToScope()
    
  
  def __init__(self, fn=None, unitSystem="cgs"):
    if fn is None:
      # Use default file.
      fn = os.path.join(os.path.dirname(_PyA.__file__), "constants/cdat.dat")
    self.inventory = {}
    self._unitSystem = unitSystem
    self.load(fn)
    
