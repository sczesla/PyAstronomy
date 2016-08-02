from __future__ import print_function, division
from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy import pyaC
import os
import gzip
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE 
import six

class Baraffe98Tracks:
  """
    Provide access to the evolutionary tracks of Baraffe et al. 98.
    
    Downloads the table (tab1-3.dat) pertaining to the paper
    Baraffe et al. 98, A&A 337, 403 and offers some helper
    method to access the individual tracks.
  """
  
  def _identifyModels(self):
    """
      Find all independent models in the data.
      
      Sets up the `models` attribute, which maps
      a model descriptor (metallicity, Y, L (mixing length),
      mass) to an array of indices, which specifies the
      fraction of the data containing that model.
    """
    self._models = {}
    metals = np.unique(self.dat[::,0])
    ys = np.unique(self.dat[::,1])
    ls = np.unique(self.dat[::,2])
    masses = np.unique(self.dat[::,3])
    nl = pyaC.NestedLoop([len(metals), len(ys), len(ls), len(masses)])
    for i in nl:
      indi = np.where(np.logical_and(self.dat[::,0] == metals[i[0]], \
                      np.logical_and(self.dat[::,1] == ys[i[1]], \
                      np.logical_and(self.dat[::,2] == ls[i[2]], self.dat[::,3] == masses[i[3]]))))[0]
      if len(indi) > 0:
        self._models[(metals[i[0]], ys[i[1]], ls[i[2]], masses[i[3]])] = indi.copy()
  
  def __init__(self):
    pfs = pp.pyaFS.PyAFS()
    dfn = os.path.join("pyasl", "resBased", "baraffe98tracks.dat.gz")
    pfs.downloadToFile("ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/337/403/tab1-3.dat.gz", \
                       dfn, verbose=True, clobber=False)
    self.dat = np.loadtxt(pfs.requestFile(dfn, openMethod=gzip.open))
    self._identifyModels()
    
    self.colToName = {0:"Met", 1:"Y", 2:"Lmix", 3:"Mass", 4:"Age", 5:"Teff", 6:"logg", 7:"Mbol", \
                      8:"MV", 9:"MR", 10:"MI", 11:"MJ", 12:"MH", 13:"MK"}
    self.nameToCol = {}
    for k, v in six.iteritems(self.colToName):
      self.nameToCol[v] = k
    
  def getUniqueValues(self, colName):
    """
      Find unique values in a column.
      
      Parameters
      ----------
      colName : string, {"Met","Y","Lmix","Mass","Age","Teff","logg","Mbol","MV","MR","MI","MJ","MH","MK"}
          The name of the column.
      
      Returns
      -------
      Unique values : array
          All unique values in that column.
    """
    if not colName in six.iterkeys(self.nameToCol):
      raise(PE.PyAValError("No column named '" + str(colName) + "'.", \
            solution="Choose one of: " + str(self.nameToCol.keys())))
    
    return np.unique(self.dat[::,self.nameToCol[colName]])
  
  def getModelData(self, model, columns=None):
    """
      Find data for a specific model.
      
      Parameters
      ----------
      model : tuple of float
          Specifies the model as a tuple containing four
          entries, which define: metallicity, initial
          helium mass fraction (Y), initial mixing length
          parameter (Lmix), and mass.
      columns : list of strings, optional
          A list of column names, which shall be included
          in the result. The default is including all
          columns.
      
      Returns
      -------
      Model : recarray
          The data pertaining to that particular model
          as a numpy recarray.
    """
    if not model in self._models:
      raise(PE.PyAValError("No model with parameters: " + str(model), \
            solution="Use, e.g., `getAvailableValues` to check what is available."))
    if columns is None:
      columns = list(self.nameToCol.keys())
    # Check validity of column names
    for col in columns:
      if not col in self.nameToCol.keys():
        raise(PE.PyAValError("No column named '" + str(col) + "'.", \
              solution="Choose one of: " + str(self.nameToCol.keys()), \
              where="getModel"))
    # There is a valid model, return it as a recarray
    dt = []
    for col in columns:
      dt.append((col, np.float))
    result = np.recarray(shape=(len(self._models[model]),), \
                         dtype=dt)
    for col in columns:
      result[col] = self.dat[self._models[model], self.nameToCol[col]]
    return result

  def findModels(self, Met=None, Y=None, Lmix=None, Mass=None):
    """
      Find models with given parameters.
      
      Parameters
      ----------
      Met : float, optional
          The metallicity.
      Y : float, optional
          The initial helium mass fraction.
      Lmix : float, optional
          The initial mixing length parameter.
      Mass : float, optional
          The mass [solar masses]
      
      Returns
      -------
      Models : list of tuples
          A list of models specified as tuples
          defining: metallicity, Y, Lmix, and mass.
    """
    result = []
    for m in self._models:
      if Met is not None:
        if Met != m[0]: continue
      if Y is not None:
        if Y != m[1]: continue
      if Lmix is not None:
        if Lmix != m[2]: continue
      if Mass is not None:
        if Mass != m[3]: continue
      result.append(m)
    return result
