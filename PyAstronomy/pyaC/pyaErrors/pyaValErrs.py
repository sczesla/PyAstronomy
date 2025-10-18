from __future__ import absolute_import
from .pyaErrTemplate import PyaErrTemplate

class PyAValError(PyaErrTemplate):
  
  def __init__(self, what, **keys):
    """
      When to be raised?
      
      Whenever an unexpected value is encountered.
    """
    PyaErrTemplate.__init__(self, what, "PyA Value Error", **keys)


class PyANameClash(PyaErrTemplate):

  def __init__(self, what, **keys):
    """
      When to be raised?
      
      Whenever an unexpected doubling of names (e.g., key values in a dictionary) is encountered. 
    """
    PyaErrTemplate.__init__(self, what, "PyA Name Clash", **keys)


class PyAFloatingPointError(PyaErrTemplate):

  def __init__(self, what, **keys):
    """
      When to be raised?
      
      Whenever a floating point error is caught.
    """
    PyaErrTemplate.__init__(self, what, "PyA Floating Point Error", **keys)