from __future__ import print_function
import copy
import six

class PyaErrTemplate(Exception):
  
  def __init__(self, what, errType, where=None, why=None, solution=None, addInfo=None):
    """
      The PyA error class template.
      
      Parameters:
       - `what` - What has happened?
       - `errType` - Description of the type of error.
       - `where` - Where did it happen?
       - `why` - Why did it happen?
       - `solution` - string or list of strings, How can the problem be solved?
       - `addInfo` - Whatever additional information may be available.
      
      All parameters are strings of not stated otherwise.
      Of the parameters, only the first (what) and second (errType) are mandatory; the latter
      should be provided by every more specialized exception class.
      
      .. note:: The template should never be raised directly, but specialized derived classes
                should be used.
    """
    self.errorType = errType
    self.what = what
    self.where = where
    self.why = why
    # Check whether solution is iterable
    self.solution = solution
    if solution is not None:
      if not isinstance(solution, six.string_types):
        self.solution = copy.deepcopy(solution)
      else:
        self.solution = [solution]    
    self.addInfo = addInfo
  
  def __str__(self, head=True):
    """
      Return error message as string.
      
      Parameters:
        - `head` - boolean, If True, a header will be printed says "PyA error". Setting False is, e.g.,
                   used in printing warnings.
    """
    e = ""
    if head:
      e = "\n"
      e += "---------------------\n"
      e += "A PyA error occurred:\n"
      e += "---------------------\n"
      e += "Type of error: " + self.errorType + "\n"
    e += "What happened?\n"
    e += "    " + self.what + "\n"
    if not self.where is None:
      e += "Where did it happen?\n"
      e += "    " + self.where + "\n"
    if not self.why is None:
      e += "Why did it happen?\n"
      e += "    " + self.why + "\n"
    if not self.solution is None:
      e += "What are possible solutions?\n"
      for s in self.solution:
        e += "  - " + s + "\n"
    if not self.addInfo is None:
      e += "Additional information:\n"
      e += "    " + self.addInfo + "\n"
    return e