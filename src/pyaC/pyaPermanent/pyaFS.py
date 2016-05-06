from __future__ import print_function, division
import os
from PyAstronomy.pyaC import pyaErrors as PE
from .pyaConfig import *
import six.moves.urllib as urllib

class PyAFS:
  """
    Manage access to PyA's data directory.
    
    This class provides a convenient interface to create,
    read, and modify files in PyA's data directory.
  """

  def _analyzeFilename(self, fn, pedantic=True):
    """
      Split and check file or directory name
      
      Parameters
      ----------
      fn : string
          The file or directory name. Either absolute or relative
          to PyA's data path.
      pedantic : boolean, optional
          If True, the method will check whether the file is
          located within PyA's data path and raise an exception
          otherwise. Default is True.
      
      Returns
      -------
      Components : dictionary
          A dictionary containing "basename", "relPath", "absPath",
          and "fullname", where "relPath" is the path relative to
          PyA's data directory. "fullname" contains the total
          filename with absolute path. 
    """
    path, basename = os.path.split(fn)
    path = os.path.normpath(path)
    isAbsPath = os.path.isabs(path)
    if not isAbsPath:
      fullpath = os.path.join(self.dpath, path)
    else:
      fullpath = path
    fullpath = os.path.normpath(fullpath)
    if not fullpath.startswith(self.dpath) and pedantic:
      raise(PE.PyAValError("PyAFS has been asked to access the file '" + fn + "', which is not located in " + \
                           "a subdirectory of PyA's data path. Forbidden access.", \
                           solution="Use another method to access the file or check whether name is correct."))
    relativePath = os.path.relpath(fullpath, self.dpath)
    fullname = os.path.join(fullpath, basename)
    return {"basename":basename, "relPath":relativePath, "absPath":fullpath, "fullname":fullname}
  
  def createSubfolder(self, folder):
    """
      Create subfolders (relative to data root directory).
      
      Ignores already existing folders. This method can handle
      recursive directory creation.
    """
    # Assume it is a path. Append path separator
    cfold = self._analyzeFilename(folder + os.sep)
    if os.path.isdir(cfold["absPath"]):
      return
    os.makedirs(cfold["absPath"])

  def composeFilename(self, rfn):
    """
      Returns full filename.
      
      Parameters
      ----------
      rfn : string
          Relative file name (in PyA data path).
      
      Returns
      -------
      filename : string
          Absolute filename.
    """
    return os.path.join(self.dpath, rfn)

  def removeFile(self, fn):
    """
      Remove a file from PyA's data directory.
      
      If the file does not exist, the request is ignored.
      
      Parameters
      ----------
      fn : string
          The filename either absolute or relative.
    """
    ana = self._analyzeFilename(fn)
    if os.path.isfile(ana["fullname"]):
      os.remove(ana["fullname"])
  
  def fileExists(self, fn):
    """
      Check whether file exists
      
      Parameters
      ----------
      fn : string
          The filename either relative to PyA's data path or
          absolute.
      
      Returns
      -------
      Flag : boolean
          True if file exists and False otherwise.
    """
    ana = self._analyzeFilename(fn, False)
    return os.path.isfile(ana["fullname"])
  
  def _checkOnline(self, url="http://www.google.com", raiseNOC=True):
    """
      Check whether network can be reached.
      
      Parameters
      ----------
      url : string, optional
          The reference URL tried to be reached to check
          network availability.
      raiseNOC : boolean, optional
          If True (default), raise an exception if network
          cannot be reached.
      
      Returns
      -------
      Online flag : boolean
          True if network could be reached.
    """
    online = True
    try:
      _ = urllib.request.urlopen(url)
    except Exception as e:
      # No connection seems to be available
      online = False
      if raiseNOC:
        raise(PE.PyANetworkError("You seem to be offline. Could not reach URL: '" + str(url) + "'.", \
                                 solution="Get online", \
                                 tbfe=e))
    return online
      
  def _checkContext(self, url="http://www.google.com", context=None, raiseOther=True):
    """
      Check whether `context` parameter is available in urllib.
      
      Parameters
      ----------
      url : string, optional
          The reference URL tried to be reached to check
          network availability.
      context : context, optional
          Context used in attempt. Default is None.
      raiseOther : boolean, optional
          If True (default), an exception is raised if an error other
          than the anticipated TypeError is raised.
    """
    contextWorks = True
    try:
      _ = urllib.request.urlopen(url, context=context)
    except TypeError as e:
      # Context appears not to work
      contextWorks = False
    except Exception as e:
      # Another error occurred 
      if raiseOther:
        raise(PE.PyANetworkError("Unknown network error occurred.", \
                                 tbfe=e))
    return contextWorks
  
  def downloadToFile(self, url, fn, clobber=False, verbose=True, openMethod=open, context=None):
    """
      Download content from URL.
      
      Parameters
      ----------
      url : string
          The location of the content.
      fn : string
          The relative or absolute name of the file
          to which the download shall be saved.
      clobber : boolean, optional
          If True, an existing file will be overwritten.
      verbose : boolean, optional
          If True, information on the download will be
          printed to the screen.
      openMethod : callable
          The method used to open the file to write to (default is
          open, other choices may be gzip.open or io.ipen)
      context : ssl context
          SSL context parameter handed to urlopen.
    """
    if self.fileExists(fn) and (clobber == False):
      return
    
    def download(url, context, nocontext=False):
      if not nocontext:
        # Use context
        response = urllib.request.urlopen(url, context=context)
      else:
        # Disregard context
        response = urllib.request.urlopen(url)
      data = response.read()     # a `bytes` object
      self.requestFile(fn, 'wb', openMethod).write(data)
    
    ana = self._analyzeFilename(fn, True)
    self.touchFile(ana["fullname"])
    try:
      if verbose:
        print("PyA download info:")
        print("  - Downloading from URL: " + str(url))
      download(url, context)
    except (KeyboardInterrupt, SystemExit):
      self.removeFile(ana["fullname"])
      raise
    except TypeError as e:
      # Possibly, context is not supported
      cs = self._checkContext()
      if not cs:
        # Network is all right, but context parameter must
        # not be specified.
        if verbose:
          print("PyA download info:")
          print("  - Downloading from URL: " + str(url) + ", (no context)")
        download(url, context, nocontext=True)
    except Exception as e:
      self.removeFile(ana["fullname"])
      sols = ["Check whether URL exists and is spelled correctly."]
      # Check whether network can be reached.
      netreach = self._checkOnline(url, raiseNOC=False)
      if not netreach:
        sols.append("Network could not be reached. Check your network status.")
      raise(PE.PyADownloadError("Could not download data from URL: " + str(url) + ".\n", \
            solution=sols, \
            tbfe=e, \
            addInfo="Could network be reached (online)? " + {True:"yes", False:"No"}[netreach]))
    if verbose:
      print("  - Downloaded " + str(os.path.getsize(ana["fullname"]) / 1000.0) + " kb")
      print("    to file: " + ana["fullname"])

  def requestFile(self, relName, mode='r', openMethod=open, *args, **kwargs):
    """
      Obtain a file object within the PyA data path.
      
      This method opens and creates files within PyA's data path.
      If the directory in which a file shall be created does not
      yet exist, it will be created.
      
      Parameters
      ----------
      relName : string
          The filename. Usually, it will be given relative to the PyA
          data path, but it can also be an absolute filename.
      mode : string, optional
          The opening mode (e.g., "r" or "w"). This flag is given to
          the `openMethod`. The default is "r" for "read".
      openMethod : method, optional
          The method used to create the file object. The default
          is Python's built-in `open` method. Another example
          could be `gzip.open`.
      
      Returns
      -------
      The file object.
    """
    writingMode = (('w' in mode) or ('a' in mode))
    ana = self._analyzeFilename(relName, pedantic=writingMode)
    fullName = ana["fullname"]
    if writingMode:
      self.createSubfolder(ana["absPath"])
    try:
      fileObject = openMethod(fullName, mode=mode, *args, **kwargs)
    except IOError as ioe:
      raise(PE.PyAValError("Could not access file (relative name): "+relName+", full name: "+fullName + "\n" + \
                           "Caught exception: "+str(ioe)))
    return fileObject

  def touchFile(self, relName):
    """
      Create an empty file in PyA's data directory.
      
      Parameters
      ----------
      relName : string
          The filename. Usually, it will be given relative to the PyA
          data path, but it can also be an absolute filename.
    """
    fullName = self._analyzeFilename(relName)["fullname"]
    if not os.path.isfile(fullName):
      f = self.requestFile(fullName, 'w')
      f.close()
    return
    
  def __init__(self):
    self.conf = PyAConfig()
    # The next call will throw exception, if there is
    # no valid date directory.
    self.dpath = self.conf.getDataRoot()
