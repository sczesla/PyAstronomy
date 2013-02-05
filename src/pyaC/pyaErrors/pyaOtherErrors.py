from pyaErrTemplate import PyaErrTemplate

class PyAOrderError(PyaErrTemplate):
  
  def __init__(self, what, **keys):
    """
      When to be raised?
      
      Whenever operations seem to be carried out in the
      wrong order.
    """
    PyaErrTemplate.__init__(self, what, "PyA Order Error", **keys)


class PyAParameterConflict(PyaErrTemplate):
  
  def __init__(self, what, **keys):
    """
      When to be raised?
      
      This exception should be raised when conflicting/mutually exclusive
      parameters are received.
    """
    PyaErrTemplate.__init__(self, what, "PyA Parameter Conflict", **keys)


class PyANetworkError(PyaErrTemplate):

  def __init__(self, what, **keys):
    """
      When to be raised?
      
      This exception should be raised when an action through the network
      has failed (e.g., a download).
    """
    PyaErrTemplate.__init__(self, what, "PyA Network Error", **keys)
  

class PyANotImplemented(PyaErrTemplate):
  
  def __init__(self, what, **keys):
    """
      When to be raised?
      
      This exception should be raised when the function/method called has not been implemented. \
      Such a situation often occurs when a member is to be implemented in a derived class \
      (abstract base class concept in c++).
    """
    PyaErrTemplate.__init__(self, what, "PyA Not Implemented", **keys)
    

class PyARequiredImport(PyaErrTemplate, ImportError):
  
  def __init__(self, what, **keys):
    """
      When to be raised?
      
      If a definitely needed package (e.g., numpy) cannot be imported.
    """
    PyaErrTemplate.__init__(self, what, "PyA import error", **keys)


class PyAUnclassifiedError(PyaErrTemplate):
  
  def __init__(self, what, **keys):
    """
      When to be raised?
      
      If an error occurred that cannot or shall not be specified further. 
    """
    PyaErrTemplate.__init__(self, what, "PyA unclassified error", **keys)


class PyADownloadError(PyaErrTemplate):
  
  def __init__(self, what, **keys):
    """
      When to be raised?
      
      A download could not successfully be carried out.
    """
    PyaErrTemplate.__init__(self, what, "PyA download error", **keys)