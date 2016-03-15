from . import pyaErrors as PE
import importlib

def pyaimport(mn, pak, globs, cra=True, rf=True):
  """
  """
  try:
    n = importlib.import_module("." + mn, pak)
    try:
      # If __all__ is defined, use it!
      to_import = n.__all__
    except AttributeError:
      # Otherwise, import all names not starting with an underscore
      to_import = [name for name in n.__dict__ if not name.startswith('_')]
      if cra:
        # Check reassignment
        for name in to_import:
          if name in globs:
            print("Reassigning: ", name, " from ", n)
      
    globs.update({name: n.__dict__[name] for name in to_import})
  except:
    if rf:
      # Report failure
      print("FAILED to import module '" + str(mn) + "' from package '" + str(pak) + "'.")

class ImportCheck:
  
  def __init__(self, modules, required=None):
    """
      Checks whether individual modules can be imported.
    
      Parameters
      ----------
      modules : list of strings
          The names of the modules to be checked.
      required : list of strings, optional
          List of modules which are required. If not found,
          an import error is raised.             
    """
    if required is None:
      required = []
    # List of required modules that could not be imported,
    # list of all modules that could not be imported
    self.requiredFail = []
    self.fail = []
    
    self.check = {}
    for module in modules:
      self.check[module] = True
      try:
        __import__(module)
      except ImportError:
        self.check[module] = False
        self.fail.append(module)
      except PE.PyARequiredImport:
        self.check[module] = False
        self.fail.append(module)
      
      if (not self.check[module]) and (module in required):
        self.requiredFail.append(module)
      
    if len(self.requiredFail) > 0:
      ms = ', '.join(self.requiredFail)
      raise(PE.PyARequiredImport("Could not import required module(s): "+ms, solution="Please install "+ms))
        
  