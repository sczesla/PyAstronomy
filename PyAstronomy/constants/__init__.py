from __future__ import print_function, division
from PyAstronomy.pyaC import ImportCheck as _ImportCheck
from PyAstronomy.pyaC import pyaErrors as _PE
import PyAstronomy as _PyA
import os, sys
from six.moves import configparser as _ConfigParser
import numpy
import six

_ic = _ImportCheck(["quantities"], ["quantities"])

import quantities as _pq

_systems = ["SI", "cgs"]
_unitSystem = "cgs"
inventory = {}

from .itc import PyAConstants, _displaySummary



def summary():
  """
    Print a summary of available constants to screen.
  """
  _displaySummary(sys.modules[__name__])

    
def constantDetails(const):
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
  if not const in inventory:
    raise(_PE.PyAValError("No such constant: '" + str(const) + "'."))
  cmod = sys.modules[__name__]
  print("Constant              : '" + const + "'")
  print("Description           : " + inventory[const]["descr"])
  print("Current value and unit: " + str(getattr(cmod, "f_"+const)))
  print("Uncertainty           : " + str(getattr(cmod, "f_"+const+"_err")))
  print("Source                : " + inventory[const]["source"])
  return inventory[const]

def inUnitsOf(const, units):
  """
    Get a constant and apply unit conversion.
    
    Returns
    -------
    Constant : Quantity object (`quantities` package)
        The value and unit of the constant.
  """
  cmod = sys.modules[__name__]
  if not hasattr(cmod, "f_" + const):
    raise(_PE.PyAValError("No such constant: '" + str(const) + "'"))
  val = getattr(cmod, "f_" + const)
  return val.rescale(units)

def _inventoryToModule():
  """
    Translate information in `inventory` into module attributes.
  """
  cmod = sys.modules[__name__]
  for k, u in six.iteritems(inventory):
    val = _pq.Quantity(numpy.float64(u["valueSI"]), u["units"]["SI"])
    value = val.rescale(u["units"][_unitSystem])
    setattr(cmod, u["symbol"], value.magnitude)
    setattr(cmod, "f_" + u["symbol"], value)
    err = _pq.Quantity(numpy.float64(u["errSI"]), u["units"]["SI"])
    err = err.rescale(u["units"][_unitSystem])
    setattr(cmod, "f_" + u["symbol"] + "_err", err)

def cleanUp():
  """
    Clean up loaded constants.
    
    Empties the `inventory` and deletes all
    related module attributes.
  """
  global inventory
  cmod = sys.modules[__name__]
  for k, u in six.iteritems(inventory):
    delattr(cmod, u["symbol"])
    delattr(cmod, "f_" + u["symbol"])
    delattr(cmod, "f_" + u["symbol"] + "_err")
  inventory = {}

def setSystem(system):
  """
    Set the unit system.
    
    Parameters
    ----------
    system : string, {SI, cgs}
        The unit system in use.
  """
  global _unitSystem
  if not system in _systems:
    raise(_PE.PyAValError("No such system: '" + system + "'.", \
                          solution="Choose one of: "+str(_systems)))
  _unitSystem = system
  _inventoryToModule()

def getSystem():
  """
    Get current unit system.
    
    Returns
    -------
    System : string, {SI, erg}
        The unit system currently in use.
  """
  global _unitSystem
  return _unitSystem

def load(fn, constClobber=False):
  """
    Load information about constants from file.
    
    A file containing information about constants is a
    simple Python configuration file. Each constant
    must be described in its own (arbitrarily named)
    section.
    
    An example for a valid entry looks like this: 
    
    ::
    
      [NewtonGravitation]
      descr   = Newton's constant of gravitation
      symbol  = G
      valueSI = 6.67428e-11 
      errSI   = 6.7e-15
      unitSI  = m^3 * kg^-1 * s^-2
      unitcgs = cm^3 * g^-1 * s^-2
      source  = http://maia.usno.navy.mil/NSFA/IAU2009_consts.html
    
    Note that all fields are mandatory.
    
    The units (and conversion) are managed by the `quantities`
    package.
    
    Parameters
    ----------
    fn : string
        Name of the file.
    constClobber : boolean, optional
        If this is set to true, previously loaded constants will be
        overwritten by newly loaded definitions. Otherwise, an error
        will be issued.
  """
  global inventory
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
      if (nc["symbol"] in inventory) and (not constClobber):
        raise(_PE.PyANameClash("A constant with the symbol '" + nc["symbol"] + "' is already defined.", \
                               solution=["Choose another symbol.", "Set the `constClobber` flag to True to overwrite existing definitions of constants."])) 
      inventory[nc["symbol"]] = nc
    _inventoryToModule()



load(os.path.join(os.path.dirname(_PyA.__file__), "constants/cdat.dat"))