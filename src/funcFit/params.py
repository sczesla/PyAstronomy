# -*- coding: utf-8 -*-
from __future__ import print_function, division
import re
from PyAstronomy.pyaC import pyaErrors as PE
import pickle
import os
import uuid
import six
import six.moves as smo


def equal(dependsOn):
  return dependsOn


class Params:
  """
    Manage a set of parameters.
    
    This class provides a framework for managing parameter values
    and several associated aspects such as restrictions and relation.
    Parameter management is implemented in view of later use by a
    a fit routine.
    
    Parameters
    ----------
    paramNames : list of strings
        Specifies the names of the parameters to be managed.
    
    
    Attributes
    ----------
    __params : dict
        A dictionary containing entries of the kind [variable-name:value].
        The `__params` dictionary may only be access by the
        assignValue function to prevent causing inconsistencies,
        which may occur especially if relations exist.                     .
    paramNum : dict
        A dictionary assigning every parameter name to a number
        (e.g., paramNum[2] = "XYZ"). Such a numbering is important to
        guarantee the correct order of variables in
        the fitting process. The ordering is the same as the order in
        which the constructor obtains the parameter names.
    isFree : dict
        A dictionary associating each parameter with a boolean saying
        whether it shall be considered a free parameter during
        optimization.
    isRestricted : dict
        A dictionary associating each parameter with a boolean
        saying whether there is a restriction to the allowed
        range (see `restrictions`).
    restrictions : dict
        A dictionary containing entries of the form
        [parName:[lower, upper]. Here 'lower' and 'upper' are floats
        defining lower and upper limits for the variable's value.
    relations : dict
        Parameters may be interrelated, e.g., A = B + C.
        This dictionary holds the definition of such relations.
        In particular, the entries have the form
        {"b":[ ["a", func, ["b", "c"]], [...] ]...}.
        This means that a = func(b, c). func is
        a function pointer. Whenever b is assigned a new value,
        'a' has to be updated
        as well. Note that it is important that the independent
        variables "know" about the relation, because if their value
        changes, the value of the dependent
        variable (in this case `a`) has to be updated.
    conditionalRestrictions : dict
        A dictionary holding the `conditional restrictions', i.e.,
        complex restrictions, which may, e.g., depend in the values
        of other parameters. The dictionary key is a unique ID
        generated, when a conditional restriction is added. For each
        key, the dictionary holds a tuple of which the first entry
        is a list holding the names of the parameters on which the
        conditional restrictions depends and the second is a
        callable, which is called with the values of those parameters
        specified in the first entry. The callable must return a
        float that specifies the penalty (or reward) depending on
        the given parameter values. Because conditional restrictions
        are referred to using a unique ID, their name (i.e., ID) does
        not change if models are combined.
    
    Notes
    -----
    
    Different models make different demands on the variables.
    For example, only certain ranges may be valid, some are constant
    and others not, or there may be a functional dependence between
    different variables. This class provides a framework to
    manage a parameter set appropriately.
    
    Depending on what kind of model is currently adapted, the number,
    names, allowed ranges, and interdependencies of variables can differ.
    On creation, this class is given a list with parameter names
    to manage. Those can than be assigned values. Parameters can be
    "thawed", i.e., regarded free during the fit process, or
    frozen. The allowed ranges can be restricted either on one or
    both sides, and interdependencies can be defined via the `relate`
    method.
  """

  def __add__(self, right):
    """
      This member allows to combine the properties of two Param class \
      instances. All parameter properties will be copied into the new \
      Parameter class instance, which is returned by this function.
      
      .. caution:: All parameter names must be unique!
    """
    paLeft = self.parameters()
    paRight = right.parameters()
    # Check whether parameters have unique names
    for p in paLeft:
      if p in paRight:
        raise(PE.PyANameClash("At least two parameters share the same name.", where="Params::__add__"))
    # List of new names
    nlist = list(paLeft.keys()); nlist.extend(paRight.keys())
    # Create new parameter objects
    result = Params(nlist)
    # Combine properties
    result.paramNum = self.paramNum.copy()
    pn = len(result.paramNum)
    for i in smo.range(len(right.paramNum)):
      result.paramNum[pn + i] = right.paramNum[i]
    result.isFree = self.isFree.copy(); result.isFree.update(right.isFree)
    result.isRestricted = self.isRestricted.copy(); result.isRestricted.update(right.isRestricted)
    result.restrictions = self.restrictions.copy(); result.restrictions.update(right.restrictions)
    result.relations = self.relations.copy(); result.relations.update(right.relations)
    result.conditionalRestrictions = self.conditionalRestrictions.copy()
    result.conditionalRestrictions.update(right.conditionalRestrictions)
    result.assignValue(self.parameters())
    result.assignValue(right.parameters())
    return result
  
  def renameParameter(self, old, new):
    """
      Rename an existing parameter.
      
      Parameters
      ----------
      old : string
          The existing (old) name.
      new : string
          The new name.
    """
    self.__checkForParam(old)
    if new in self.__params:
      raise(PE.PyAValError("Parameter already exists: "+new, where="Params::renameParameter"))
    self.__params[new] = self.__params[old]; del self.__params[old]
    for i in smo.range(len(self.paramNum)):
      if self.paramNum[i] == old:
        self.paramNum[i] = new
        break
    self.isFree[new] = self.isFree[old]; del self.isFree[old]
    self.isRestricted[new] = self.isRestricted[old]; del self.isRestricted[old]
    self.restrictions[new] = self.restrictions[old]; del self.restrictions[old]
    self.relations[new] = self.relations[old]; del self.relations[old]
    # Loop through relations, search, and replace occurrences of the old name.
    for p in six.iterkeys(self.__params):
      relations = self.relations[p]
      if relations == []: continue
      for k in smo.range(len(relations)):
        relat = relations[k]
        if relat[0] == old:
          relat[0] = new
        for i in smo.range(len(relat[2])):
          if relat[2][i] == old:
            relat[2][i] = new
        self.relations[p][k] = relat
    # Loop over conditional restrictions and replace occurrences of the old names
    for name, v in six.iteritems(self.conditionalRestrictions):
      for i, p in enumerate(v[0]):
        # Loop over input-parameter names and replace if necessary
        if p == old:
          # Needs to be replaced
          self.conditionalRestrictions[name][0][i] = new
  
  def __init__(self, paramNames):
    i = 0
    # Define instance properties
    self.__params = {}
    self.paramNum = {}
    self.isFree = {}
    self.isRestricted = {}
    self.restrictions = {}
    self.relations = {}
    self.conditionalRestrictions = {}
    for n in paramNames:
      self.__params[n] = 0.0
      self.isFree[n] = False
      self.paramNum[i] = n
      i += 1
      # Firstly, the lower bound, secondly, the upper bound
      self.isRestricted[n] = [False, False]
      self.restrictions[n] = [None, None]
      self.relations[n] = []
  
  def addConditionalRestriction(self, pars, func):
    """
      Define a conditional restriction.
      
      Conditional restrictions can be used to modify the
      behavior in a more complex manner. For instance,
      penalties can be added to the objective function
      depending on the relation of one or more parameters.
      
      The given function is evaluated in each iteration
      and its return value (a float) is added to the
      objective function (e.g., chi square).
      
      Parameters
      ----------
      pars : list of strings
          The names of the parameters the given function
          needs as input.
      func : callable object
          This callable object must take the specified
          parameters (in that exact order) as input. It
          must return a float, which is added to the
          value of the objective function.
      
      Returns
      -------
      identifier : string
          A unique ID used to refer to the conditional
          restriction.
    """
    for p in pars:
      self.__checkForParam(p)
    name = uuid.uuid4()
    self.conditionalRestrictions[name] = (pars[:], func)
    return name
  
  def removeConditionalRestriction(self, id):
    """
      Remove an existing conditional constraint.
      
      Parameters
      ----------
      id : string
          The identifier used to refer to the conditional
          constraint (returned by `addConditionalRestriction`).
    """
    if not id in self.conditionalRestrictions:
      raise(PE.PyAValError("No conditional restriction with ID '" + str(id)))
    del self.conditionalRestrictions[id]
  
  def showConditionalRestrictions(self, toScreen=True):
    """
      Show conditional restrictions.
      
      Parameters
      ----------
      toScreen : boolean, optional
          If True (default), the output is written to stdout.
      
      Returns
      -------
      Output : list of strings
          The output as a list of strings.
    """
    ll = []
    nc = 80
    ll.append("-" * nc)
    ll.append("    Conditional restrictions")
    ll.append("-" * nc)
    for name, v in six.iteritems(self.conditionalRestrictions):
      s = "ID: " + str(name) + ", parameters: "
      s += ', '.join(v[0])
      ll.append(s)
    ll.append("-" * nc)
    if toScreen:
      for l in ll:
        print(l)
    return ll
      
  
  def applyConditionalRestrictions(self, fullout=False):
    """
      Apply all conditional restrictions.
      
      Parameters
      ----------
      fullout : boolean, optional
          If True, a dictionary holding the values
          of the individual restrictions is returned.
          The IDs are used as dictionary keys. The
          default is False.
      
      Returns
      -------
      Modification : float
          The summed value of the existing conditional restrictions.
      Individual values : dictionary, optional
          The contributions of the individual conditional
          restrictions.
    """
    result = 0.0
    if fullout:
      fo = {}
    for name, v in six.iteritems(self.conditionalRestrictions):
      vs = ()
      for p in v[0]:
        vs += (self[p],)
      val = v[1](*vs)
      if val is None:
        raise(PE.PyAValError("Conditional restriction with ID '" + str(name) + "' returned None.", \
                             solution="Did you forget to return 0.0, in case the condition is not fulfilled?"))
      if fullout:
        fo[name] = val
      result += val
    if fullout:
      return result, fo
    else:
      return result
  
  def __toList(self, x):
    """
      Return a list of x is a string and do nothing of x is already a list.
    """
    if isinstance(x, six.string_types):
      return [x]
    return x
  
  def __checkForParam(self, name):
    """
      Checks whether parameter exists and throws an exception if this is not the case.
      No return value.
    """
    if not self.hasParam(name):
      raise(PE.PyAValError("No such parameter: \""+name+"\".\n  Available parameters: "+', '.join(list(self.parameters().keys())) , where="Params::__checkForParam"))
  
  def hasParam(self, name):
    """
      Check whether parameter exists.
      
      Parameters
      ----------
      name : string
          The parameter name.
      
      Returns
      -------
      flag : boolean
          True if parameter exists, False otherwise.
    """
    if not name in self.__params:
      return False
    return True

  def availableParameters(self):
    """
      Provides a list of existing parameters.
    
      Returns
      -------
      Parameters : list of strings 
          A list with the names of available parameters.
    """
    return list(self.__params.keys())

  def parameters(self):
    """
      Obtain parameter names and values.
      
      Returns
      -------
      Name-value : dict
          A dictionary with the names and values of all parameters
          ({"parName":value, ...}).
    """
    return self.__params.copy()

  def __getitem__(self, name):
    """
      Access parameter values.
      
      Parameters
      ----------
      name : string
          Parameter name.
    
      Returns
      -------
      Value : float
          The value of the parameter specified by `name`.
    """
    self.__checkForParam(name)
    return self.__params[name]

  def assignValue(self, namval):
    """
      Define new parameter values.
      
      Parameters
      ----------
      namval : dict
          A dictionary containing ['name':'value'] pairs.
    """
    for n, v in six.iteritems(namval):
      self.__checkForParam(n)
      self.__params[n] = v
      # Check for relations
      if len(self.relations[n]) > 0:
        # Yes, there are relations
        for tie in self.relations[n]:
          arg = []
          if len(tie[2]) == 1:
            # Function takes a single argument
            # Get value of dependent variable
            arg = [self[tie[2][0]]]
          else:
            # Function takes a list of arguments
            for name in tie[2]:
              arg.append(self[name])
          # Assign value to dependent variable
          self.assignValue({tie[0]:tie[1](*arg)})
  
  def thaw(self, name):
    """
      Thaw (regard as free) a parameter.
      
      Parameters
      ----------
      name : string or list of strings
          The name(s) of the parameter(s) to be thawed.
    """
    names = self.__toList(name)
    for n in names:
      # Check whether the parameter is a dependent variable in
      # a relation
      for pn, rels in six.iteritems(self.relations):
        for rel in rels:
          if rel[0] == n:
            raise(PE.PyAParameterConflict("You tried to free the parameter '" + n + \
                                          "', which is a dependent variable in a relation.", \
                                          solution="Use 'untie' first to remove the relation."))
      self.__checkForParam(n)
      self.isFree[n] = True

  def freeze(self, name):
    """
      Freeze parameter(s) (contrary of thaw).

      Parameters
      ----------
      name : string or list of strings
          The name(s) of the parameter(s) to be frozen.
    """
    names = self.__toList(name)
    for n in names:
      self.__checkForParam(n)
      self.isFree[n] = False

  def freeParameters(self):
    """
      Get names and values of free parameters.
      
      Returns
      -------
      Free parameters: dict
           Dictionary containing the names and values
           of all free parameters ({"parName":value, ...}).
    """
    tuplist = []
    for i in smo.range(len(self.paramNum)):
      name = self.paramNum[i]
      if self.isFree[name]:
        tuplist.append((name,self.__params[name]))
    result=dict(tuplist)
    return result
  
  def frozenParameters(self):
    """
      Get names and values of frozen parameters.
      
      Returns
      -------
      Frozen parameters: dict
           Dictionary containing the names and values
           of all frozen parameters ({"parName":value, ...}).
    """
    tuplist = []
    for i in smo.range(len(self.paramNum)):
      name = self.paramNum[i]
      if not self.isFree[name]:
        tuplist.append((name,self.__params[name]))
    result=dict(tuplist)
    return result

  def numberOfFreeParams(self):
    """
      Get number of free parameters.
    
      Returns
      -------
      n : int
          The number of free parameters (determined by `isFree`).
    """
    c = 0
    for v in six.itervalues(self.isFree):
      if v: c += 1
    return c

  def freeParamNames(self):
    """
      Get the names of the free parameters.
      
      Returns
      -------
      Free parameters : list of strings
          The names of the free parameters.
          The order is determined by the
          `paramNum` attribute.
    """
    result = []
    for i in smo.range(len(self.paramNum)):
      name = self.paramNum[i]
      if self.isFree[name]: result.append(name)
    return result

  def setFreeParams(self, X):
    """
      Change the values of the free parameters.
    
      Parameters
      ----------
      X : list of floats
          Contains the values for all free parameters. Note
          that the numbering is according to the `paramNum` attribute.

      Notes
      -----
      This method is primarily implemented to be used by fit routines.
    """
    if len(X) != self.numberOfFreeParams():
      raise(PE.PyAValError("Number of supplied parameters does not match number of free parameters.", where="Params::setFreeParams"))
    c = 0
    for i in smo.range(len(self.__params)):
      name = self.paramNum[i]
      if self.isFree[name]:
        self.assignValue({name:X[c]})
        c += 1

  def getFreeParams(self):
    """
      Get values of free parameters.
      
      Returns
      -------
      Values : list of floats
          The values of the free parameters. Note that the order is
          determined by the `paramNum` attribute.
    """
    result = []
    for i in smo.range(len(self.__params)):
      name = self.paramNum[i]
      if self.isFree[name]: result.append(self.__params[name])
    return result
  
  def setRestriction(self, restricts):
    """
      Apply restrictions to parameter ranges.
      
      Parameters
      ----------
      restricts : dict
          A dictionary associating name and [lower-bound, upper-bound].
          If no boundary shall exist on one side, use 'None'.
    """
    for name,v in six.iteritems(restricts):
      self.__checkForParam(name)
      if (v[0] is not None) and (v[1] is not None):
        if v[0] >= v[1]:
          raise(PE.PyAValError("Lower bound < upper bound must be fulfilled!", where="Params::setRestriction"))
      self.isRestricted[name] = [(v[0] is not None),(v[1] is not None)]
      self.restrictions[name] = v

  def delRestriction(self, parName):
    """
      Delete restriction
      
      Parameters
      ----------
      parName : string
          Name of restricted parameter
    """
    self.__checkForParam(parName)
    self.isRestricted[parName] = [False, False]
    self.restrictions[parName] = [None, None]
    

  def getRestrictions(self):
    """
      Get all restrictions.
      
      Returns
      -------
      Restrictions : dict
          Dictionary associating parameter name and restriction
          (see `restrictions` attribute).
    """
    return self.restrictions

  def untie(self, parName, forceFree=False):
    """
      Remove all relations of parameter "parName".
      
      After this operation, the parameter no longer depends 
      on other parameters. Unless foreFree is True, the parameter
      referred to by "parName" is "frozen".
      
      Parameters
      ----------
      parName : string
          The name of the dependent variable, which shall be untied.
      forceFree : boolean
          Set parName to "free" instead of "frozen" if set to True.
    """      
    for p in six.iterkeys(self.relations):
        self.relations[p] = [r for r in self.relations[p] if r[0] != parName]           
    self.isFree[parName] = forceFree

  def relate(self, parName1, pars, func = None, force=False):
    """
      Apply functional relation between parameters.
    
      Relates parameters, i.e., par1 = func(pars). The two values are related by 'func'. \
      In order to be tied, all involved parameters have to be free.
      
      Parameters
      ----------
      parName1 : string
          The name of the dependent variable.
      pars : list of strings
          Names of independent variables.
      func : callable
          The function that determines the form of the dependence.
          The default is "equal".
      force : boolean, optional
          Set to True in order to suppress error when the dependent
          variable is frozen. 
    """
    self.__checkForParam(parName1)
    if not self.isFree[parName1] and not force:
      raise(PE.PyAValError(parName1+" is not free.", where="Params::relate", solution="Use 'thaw' to free parameter."))
    if isinstance(pars, six.string_types): pars = [pars]
    for p in pars:
      self.__checkForParam(p)
    # Check whether the dependent variable is in the list of independent ones
    if parName1 in pars:
      raise(PE.PyAValError("The dependent variable ("+parName1+") cannot be in the list of independent variables."))
    # Use equal function if None is specified.
    if func is None: func = equal
    # Append to 'relations' list
    for p in pars:
      # If parameter 'p' is set, redefine 'parName1' using function 'func' with parameters 'pars'
      self.relations[p].append( [parName1, func, pars] )
    # Update parameter values (on re-assigned the values of the independent
    # parameters, the value of the dependent parameter will be updated)
    for p in pars:
      self.assignValue({p:self.__params[p]})
    self.isFree[parName1] = False

  def getPenalty(self, penaltyFact=1e20):
    """
      Get applied penalty for current parameter set.
    
      Parameters
      ----------
      penaltyFact : float, optional
          The higher the number the higher the penalty for small deviations
          (default is 10**20).
    
      Returns
      -------
      Penalty : float
          The applied penalty for current parameter set.
      Penalties : dict
          A dictionary with a key for every parameter values
          for which a bound is violated. The value is the
          amount by which the bound (upper or lower) is
          violated.
    """
    result = {}
    totval = 0.0
    for name, r in six.iteritems(self.isRestricted):
      if r[0]:
        # There is a lower bound
        if self.__params[name] < self.restrictions[name][0]:
          result[name] = abs(self.__params[name] - self.restrictions[name][0])
          totval += result[name]
          # It cannot be smaller and greater...continue
          continue
      if r[1]:
        if self.__params[name] > self.restrictions[name][1]:
          result[name] = abs(self.__params[name] - self.restrictions[name][1])
          totval += result[name]
          continue
    totval *= penaltyFact
    return (totval, result)

  def getRelationsOf(self, parName):
    """
      Obtain relations for a parameter.
    
      Parameters
      ----------
      parName : string
          The name of the parameter of which relations shall be searched.
    
      Returns
      -------
      Relations : list of relations
          Those relations in which `parName` is the dependent variable, i.e.,
          parName = f(x,y).
    """
    self.__checkForParam(parName)
    result = []
    for k in six.iterkeys(self.__params):
      for r in self.relations[k]:
        if r[0] == parName:
          if not r in result:
            result.append(r)
    return result

  def saveState(self, fn=None, clobber=False):
    """
      Save the state of the fitting object.
      
      This method collects the parameter values, the applied restrictions,
      and the information whether parameters are free or frozen and saves
      them to the specified file (if given) using pickle.
      The saved state can be restored using
      the `restoreState` method.
      
      .. note:: Unfortunately, "relations" cannot be saved.
      
      Parameters
      ----------
      fn : string, optional
          The filename to which the state shall be written.
          If None, the output will not be written to a file.
      clobber : boolean, optional
          If True, existing files will be overwritten (default is False).
      
      Returns
      -------
      Saved data : dict
          The dictionary containing the data saved to the file.
    """
    dat = {}
    dat["parameters"] = self.parameters()
    dat["restrictions"] = self.getRestrictions()
    dat["frozenParameters"] = self.frozenParameters()
    if isinstance(fn, six.string_types):
      if os.path.isfile(fn) and (not clobber):
        raise(PE.PyANameClash("The file '"+fn+"' exists.", where="saveState", \
                              solution="Change filename or set clobber to True."))
      pickle.dump(dat, open(fn, 'wb'))
    return dat
  
  def restoreState(self, resource):
    """
      Restores parameter values from file or dictionary.
      
      Parameters
      ----------
      resource : string or dictionary
          If string, it is interpreted as filename of a pickle file holding the
          data dictionary. If dictionary, it uses the data saved in it; note that
          a valid data dictionary is returned by `saveState`.
    """
    if isinstance(resource, six.string_types):
      pd = pickle.load(open(resource, 'rb'))
    else:
      pd = resource
    # Free ALL values
    for k in six.iterkeys(self.parameters()):
      rels = self.getRelationsOf(k)
      if len(rels) == 0:
        self.thaw([k])
    # Assign parameter values
    self.assignValue(pd["parameters"])
    # Assign restrictions
    self.setRestriction(pd["restrictions"])
    # Freeze the frozen parameters
    self.freeze(pd["frozenParameters"])

  def parameterSummary(self, lines=None, toScreen=True, prefix="", onlyParams=False):
    """
      Writes a summary of the parameters in text form.
      
      Parameters
      ----------
      lines : list of strings, optional
          If given, the output will be attached to this list.
      toScreen : boolean, optional
          If False, screen output will be suppressed
          (default is True).
      prefix : string, optional
          A prefix applied to every output line (e.g., '#')
      
      @FIXME - This is antique...
          
      Attaches the text lines
      to the lines list if given. If lines is not given, 'toScreen' will be set to True
      and the result is written to the screen.
      Returns either a new list containing the text, or (if given) the 'lines' list
      with result appended.
      A 'prefix' may be specified (e.g., '#'), which will preceded every line.
    """
    if lines is None:
      lines = []
      toScreen = True
    
    # Find maximal length of parameter name
    plen = 0
    for k,v in six.iteritems(self.__params):
      if len(k) > plen: plen = len(k)
    
    # Comppars - dictionary mapping Component number -> list of variables
    compPars = {}
    for name in six.iterkeys(self.__params):
      r = re.match(".*_([0-9]+)", name)
      if r is not None:
        cn = int(r.group(1))
      else:
        cn = -1
      if cn in compPars:
        compPars[cn].append(name)
      else:
        compPars[cn] = [name]

    for k, v in six.iteritems(compPars):
      compPars[k] = sorted(v)

    components = sorted(compPars.keys())
    for com in components:
      if com != -1 and not onlyParams:
        lines.append(prefix + "Parameters for component no. "+str(com))
        lines.append(prefix + "-------------------------------------")
      for i in smo.range(len(compPars[com])):
  #      name = self.paramNum[i]
        name = compPars[com][i]
        
        restricted = (self.isRestricted[name][0] or self.isRestricted[name][1])
        line = (prefix + "Parameter: %"+str(plen)+"s, value:  %13.8g, free: %5s, restricted: %5s") % (name, self.__params[name], str(self.isFree[name]), str(restricted))
        if restricted:
          line += ", lower bound: "
          if self.restrictions[name][0] is not None:
            line += "%12.6f" % self.restrictions[name][0]
          else:
            line += "None"
          line += ", "+"upper bound: "
          if self.restrictions[name][1] is not None:
            line += "%12.6f" % self.restrictions[name][1]
          else:
            line += "None"
  
        relas = self.getRelationsOf(name)
        if len(relas) > 0:
          for r in relas:
            line += ",  Relation: "+r[0]+" = f("
            for k in range(len(r[2])-1):
              line += r[2][k] + ", "
            line += r[2][len(r[2])-1] + ")"
        lines.append(line)
    
    if toScreen:
      for l in lines:
        print(l)
    
    return lines