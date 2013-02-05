import os
from PyAstronomy.pyaC import pyaErrors as PE
from pyaConfig import *
import urllib

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
          The filename wither relative to PyA's data path or
          absolute.
      
      Returns
      -------
      Flag : boolean
          True if file exists and False otherwise.
    """
    ana = self._analyzeFilename(fn, False)
    return os.path.isfile(ana["fullname"])
    
  def downloadToFile(self, url, fn, clobber=False, verbose=True):
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
    """
    if self.fileExists(fn) and (clobber == False):
      return
    
    ana = self._analyzeFilename(fn, True)
    self.touchFile(ana["fullname"])
    try:
      if verbose:
        print "PyA download info:"
        print "  - Downloading from URL: " + str(url)
      filename, header = urllib.urlretrieve(url, ana["fullname"])
    except (KeyboardInterrupt, SystemExit):
      self.removeFile(ana["fullname"])
      raise
    except Exception as e:
      self.removeFile(ana["fullname"])
      raise(PE.PyADownloadError("Could not download data from URL: " + str(url) + ".\n" + \
            "Error message: " + str(e), \
            solution="Check whether URL exists and is spelled correctly."))
    if verbose:
      print "  - Downloaded " + str(os.path.getsize(filename) / 1000.0) + " kb"
      print "    to file: " + filename

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
    except IOError, ioe:
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
