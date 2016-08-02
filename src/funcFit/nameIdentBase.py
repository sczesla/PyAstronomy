from __future__ import print_function, division
from PyAstronomy.pyaC import pyaErrors as PE
import re
import six

class ModelNameIdentBase:
  """
    Managing the naming of model components.
    
    This class handles the names of models or model components.
    Individual names or identifiers are composed of a "root name"
    and a "component counter". The root name is supposed to be a
    concise string summarizing the type of model, while the component
    counter is used to distinguish between components
    with the same root name in composed models.
    
    Parameters
    ----------
    rootName : string, optional
        A concise name for the model (default="").
    
    Notes
    -----
    The term "specifier" is used to indicate either a string
    containing a variable name or a tuple composed of
    (property, root name, component counter), which
    specifies a variable. In some cases, parts of the specifier
    tuple may be left out.
  """
  
  def __init__(self, rootName=""):
    self._root = rootName
    self.componentCounter = 0
  
  def setComponentCounter(self, c):
    self.componentCounter = c
  
  def getComponentCounter(self):
    return self.componentCounter
    
  def setRootName(self, root):
    """
      Define the root used to build the identifier.
      
      Parameters
      ----------
      root : string
          Concise description of the model.
    """
    if re.match("^[^_^\(^\)]*$", root) is None:
      raise(PE.PyAValError("The proposed root '"+root+"' contains forbidden characters: '_()'."))
    self._root = root
  
  def identifier(self):
    """
      Return the appropriate identifier.
      
      Returns
      -------
          identifier : string
      
      Notes
      -----
      The "identifier" reads "root-name(c)",
      where the later part, "(c)", will only be present if the `componentCounter` property
      is larger than zero; in this case, 'c' is that number.
    """
    result = self._root
    if self.componentCounter != 0:
      # Append component counter in parentheses if necessary
      result += "(" + str(self.componentCounter) + ")"
    return result
  
  def getRoot(self):
    """
      Returns the root name.
    """
    return self._root
    
  def composeVariableName(self, property, rootName=None, counter=None):
    """
      Combine property, root name, and counter into a variable name.
      
      Parameters
      ----------
      property : string
          The property name.
      rootName : string, optional,
          A root name. If None, the instance's own root
          name will be used.
      counter : int, optional,
          The component counter. If None, the instance's
          counter will be used.
      
      Returns
      -------
      Variable name : string
          The resulting variable name.

      Notes
      -----
      The standard way used to compose variable names is:
      "property_rootname(counter)"
    """
    if rootName is None:
      rootName = self._root
    if counter is None:
      counter = self.componentCounter
    if rootName == "" and counter == 0:
      return property
    result = property + "_" + rootName
    if counter == 0:
      return result
    result += "("+str(counter)+")"
    return result

  def convertSpecifier(self, specifier):
    """
      Decompose specifier and compose variable name.
      
      Parameters
      ----------
      specifier : string or tuple,
          Either a string giving a variable name or a tuple
          specifying property, root name, and counter
          (the later two are optional).
      
      Returns
      -------
      Decomposed specifier : tuple
          A tuple of (variable name, property, root name,
          component counter).
    """
    component = ""
    counter = 0
    prop = ""
    if isinstance(specifier, tuple):
      prop = specifier[0]
      if len(specifier) >= 2:
        # A component is given
        component = specifier[1]
      if len(specifier) == 3:
        # A counter is also given
        counter = specifier[2]
    else:
      prop = specifier
    return (self.composeVariableName(prop, component, counter), prop, component, counter)

  def specifierToName(self, input):
    """
      Convert specifier(s) to variable names.
      
      Parameter
      ---------
      input : Specifier, list, or dict,
          The input is a single specifier, a list of specifiers, or
          a dictionary in which specifiers are used as keys.
      
      Returns
      -------
      out : Same type as input,
          Converts all specifiers to variable names. The output is
          the same type as the input.
    """
    if isinstance(input, list):
      # Turn all element into variable names
      result = []
      for el in input:
        result.append(self.convertSpecifier(el)[0])
      return result
    if isinstance(input, dict):
      # Convert all keys to variable names
      result = {}
      for k, v in six.iteritems(input):
        result[self.convertSpecifier(k)[0]] = v
      return result
    # Assume it is a single specifier
    return self.convertSpecifier(input)[0]

  def decomposeVariableName(self, name):
    """
      Decomposes variable name into constituents.
    
      Expects a name of the from "property_root(c)" and returns
      the individual building blocks of the name, i.e., property,
      root name, component counter in that order. If one or more
      of these building blocks is not present, None will
      be returned instead.
      
      Parameters
      ----------
      name : string
          The variable name.
    """
    r = re.match("([^_]+)(_([^\(]*)(\(([0-9]+)\))?)?", name)
    return r.group(1), r.group(3), r.group(5)