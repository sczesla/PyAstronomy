from __future__ import print_function, division
import numpy as np
import os
import gzip
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyasl import _ic 
import six.moves as smo
from PyAstronomy.pyaC import pyaPermanent as pp
from astropy.table import Table as _atable


class MMSCETSTable(pp.PyAUpdateCycle):
  """
    Access to the "Modern Mean Stellar Color and Effective Temperature Sequence" table. 
    
    A compilation of data in the form of a table entiled `"A Modern Mean Stellar Color and Effective Temperature
    Sequence for O9V-Y0V Dwarf Stars" <http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt>`_
    is made available by E. Mamajek.
    
    In particular, the following columns are provided:
    SpT, Teff, logT, BCv, Mv, logL, B-V, Bt-Vt, U-B, V-Rc, V-Ic, V-Ks, J-H, H-K, Ks-W1, Msun, logAge, b-y,
    #SpT, M_J, M_Ks, Mbol".
    
    Much of the material (but not all) is published in Pecaut & Mamajek (2013, ApJS, 208, 9), which should
    be cited when using the data, and another reference for a fraction of the data is
    Pecaut, Mamajek, & Bubar (2012, ApJ 756, 154). More information can be found on the E. Mamajek's
    website and within the data file itself (use the getContent method to get the entire file).
    
    Entries without a valid entry are represented by NaN in this class.
    
    .. note:: The author of the table, E. Mamajek, states "Please email me if you use the table in your research
              and/or have any questions. - EEM". 
  """
  
  def __init__(self, skipUpdate=False, forceUpdate=False):
    
    if not _ic.check["astropy"]:
      raise(PE.PyARequiredImport("The 'astropy' package is not installed. astropy is required to read VO tables.", \
                                 solution="Please install 'astropy'."))

    
    configFilename = os.path.join("pyasl", "pyasl", "mamatable.cfg")
    pp.PyAUpdateCycle.__init__(self, configFilename, "MamajekTableUpdate", updateCycle=30)
    self.dataFileName = os.path.join("pyasl", "pyasl", "mamatable.gz")
    self._fs = pp.PyAFS()
    if forceUpdate:
      self._update(self._download)
    elif (self.needsUpdate() or (not self._fs.fileExists(self.dataFileName))) and (not skipUpdate):
      # Download data if data file does not exist or
      # regular update is indicated
      self._update(self._download)
    self._readData()
  
  def getTable(self):
    """
      Get table as astropy table object.
    """
    return self._tdata.copy()

  def getContent(self):
    """
      Get content of data file.
      
      The file contains a lot more information about
      the sources of the data, a change log etc..
      
      Returns
      -------
      content : list of strings
          The content of the file
    """
    return self._lines[:]
  
  def availableColumns(self):
    """
      Returns a list of available column names.
    """
    return self._tdata.columns[:]
  
  def getColumn(self, col, asarray=True):
    """
      Get data of specific column.
      
      Parameter
      ---------
      col : string
          Column specifier
      asarray : boolean, optional
          If True (default), the data will be converted into
          a numpy array. Otherwise, an astropy Column object
          is returned.
          
      Returns
      -------
      column data : array or column object
          Data for the specified column.
    """
    if asarray:
      return np.array(self._tdata[col])
    return self._tdata[col]
  
  def _readData(self):
    """
      Read the data from the file
    """
    self._lines = self._fs.requestFile(self.dataFileName, 'rt', gzip.open).readlines()
    # Table lines
    headerFound = False
    headerLine = None
    for n, l in enumerate(self._lines):
      if l.startswith("#SpT") and (not headerFound):
        # The header line was found
        headerLine = n
        headerFound = True
      elif l.startswith("#SpT") and headerFound:
        # This is the line after the table ends
        break
    tablelines = self._lines[headerLine:n]
    # Replace dot-placeholders by NaN
    for i in smo.range(len(tablelines)):
      for j in smo.range(5,2,-1):
        tablelines[i] = tablelines[i].replace("."*j, "NaN")
    # Read data
    self._tdata = _atable.read(tablelines, format="ascii")
      
  def _download(self):
    """
      Download data.
    """   
    self._fs.downloadToFile("http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt", self.dataFileName, clobber=True,
                            verbose=False, openMethod=gzip.open)
    