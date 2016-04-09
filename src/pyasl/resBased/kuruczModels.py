from __future__ import print_function, division
import os
import re
import pickle
import gzip
import numpy as np
from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves.urllib as urllib


class KuruczMT:
  """
    Provides access to individual models in a model grid.
    
    Parameters
    ----------
    fn : string
        Name of the Kurucz model file.
    
    Attributes
    ----------
    teffs : array
        Available effective temperatures [K] in ascending order.
    loggs : array
        Available loggs in ascending order [cgs].
    models : dictionary
        The key is a tuple of the form: (teff, logg, met)
        and the value is a list of strings representing the
        model.
    met : float
        The logarithmic metallicity of the model grid.
    
  """
  
  def _extractModelProperties(self, model):
    """
      Extract Teff, logg, and logg from a model.
      
      Returns
      -------
      Teff, logg, met : float
          Effective temperatures, logg, and log10 of metallicity
          as encoded in the model header.
    """
    r = re.match("\s*TEFF\s+([^\s]+)\s+GRAVITY\s+([^\s]+).*", model[0])
    teff = float(r.group(1))
    logg = float(r.group(2))
    r = re.match("\s*ABUNDANCE SCALE\s+([^\s]+).*", model[4])
    met = np.log10(float(r.group(1)))
    met = np.floor(met*10.0+0.5) / 10.0
    return teff, logg, met
  
  def _readModelTable(self, fn):
    """
      Read file, extract models, and store in class.
      
      Parameters
      ----------
      fn : string
          Name of the Kurucz model-grid file. 
    """
    lines = gzip.open(fn, 'rt').readlines()
    if len(lines) == 1:
      lines = lines[0].split("\r")
    models = []
    model = []
    for l in lines:
      if (l.find("TEFF") != -1) and (l.find("GRAVITY") != -1):
        # A new model starts here
        if len(model) > 0:
          # Add model to list
          models.append(model[::])
        model = []
      model.append(l.rstrip("\n"))
    # Build models dictionary
    teffs = set()
    loggs = set()
    self.met = None
    self.models = {}
    for m in models:
      props = self._extractModelProperties(m)
      self.models[props] = m[::]
      teffs.add(props[0])
      loggs.add(props[1])
      if self.met is None:
        self.met = props[2]
      elif self.met != props[2]:
        raise(PE.PyAValError("Abundance of models differ. This should not happen!\n" + \
                             "The values are: " + str(props[2]) +" and " + str(self.met), \
                             where="KuruczMT::_readModelTable"))
    # Build array of available models    
    self._modelsAvailable = np.zeros( (len(teffs), len(loggs)), dtype=np.bool)
    self.teffs = np.array(sorted(teffs))
    self.loggs = np.array(sorted(loggs))
    for i, teff in enumerate(self.teffs):
      for j, logg in enumerate(self.loggs):
        props = (teff, logg, self.met)
        if props in self.models:
          self._modelsAvailable[i, j] = True

  def availableTeffs(self):
    """
      Get available effective temperatures.
      
      Returns
      -------
      Teffs : array
          Array of available effective temperatures
          sorted in ascending order. 
    """
    return self.teffs.copy()
  
  def availableLoggs(self):
    """
      Get available logg values.
      
      Note that not all models are available
      for all logg values.
      
      Returns
      -------
      Loggs : array
          Array of available loggs sorted in
          ascending order.
    """
    return self.loggs.copy()
  
  def metallicity(self):
    """
      Get metallicity of model grid.
      
      Returns
      -------
      met : float
          Log10 of the metallicity.
    """
    return self.met
  
  def modelAvailable(self, teff, logg, met=None):
    """
      Determine whether model is available.
      
      Parameters
      ----------
      teff : float
          Effective temperatures [K]
      logg : float
          Logg [cgs]
      met : float, optional
          Logarithmic metallicity. If not given, the metallicity
          of the model grid is used.
      
      Returns
      -------
      Availability flag : boolean
          True, if the model is available.
    """
    if met is None:
      met = self.met
    return ((teff, logg, met) in self.models)
  
  def getModel(self, teff, logg, met=None):
    """
      Get a model.
      
      Parameters
      ----------
      teff : float
          Effective temperatures [K]
      logg : float
          Logg [cgs]
      met : float, optional
          Logarithmic metallicity. If not given, the metallicity
          of the model grid is used.
      
      Returns
      -------
      Model : list of strings
          The model as found on the file.
    """
    if not self.modelAvailable(teff, logg, met):
      raise(PE.PyAValError(("No model for the combination: Teff = %6.3e, logg = %6.3e, met = " + str(met)) % (teff, logg)))
    if met is None:
      met = self.met
    return self.models[(teff, logg, met)][::]
  
  def __init__(self, fn):
    if not os.path.isfile(fn):
      raise(PE.PyAValError("No such file: " + str(fn)))
    self._readModelTable(fn)
    

class KuruczModels:
  """
    Provides access to the Kurucz models.
  """
  
  def _getListOfModelGrids(self):
    """
      Get list of available grids from http://kurucz.harvard.edu/grids.html
    """
    gridFile = os.path.join("pyasl", "resBased", "kuruczMG.dat.gz")
    if not self._fs.fileExists(gridFile):
      # Grid file has to be created
      try:
        with self._fs.requestFile(gridFile, mode='w', openMethod=gzip.open) as f:
          try:
            d = urllib.request.urlopen("http://kurucz.harvard.edu/grids.html").read()
            d = d.decode("utf8")
          except urllib.error.URLError as e:
            raise(PE.PyADownloadError("Could not access URL: http://kurucz.harvard.edu/grids.html\n" + \
                                      "Error: " + str(e), \
                                      solution="Are you online?"))
          pattern = '<A HREF=\"(.*)\">(.*)</A>'
          self.grids = {}
          for grid in re.findall(pattern, d):
            if grid[0].find("grids") == -1:
              continue
            # Find the grid file
            print("Checking: ", grid[0])
            try:
              gfs = urllib.request.urlopen(grid[0]).read().decode("utf8")
            except urllib.error.HTTPError:
              # Ignore dead links
              print("  Ignoring: ", grid[0])
              continue
            except urllib.error.URLError as e: 
              print("  Error opening url: ", grid[0], ". Error message: " + str(e) + ".\n Ignoring...")
              continue
            # <img src="/icons/blank.gif" alt="[   ]" width="20" height="20"> <a href="am40ak2odfnew.dat">am40ak2odfnew.dat</a>          13-Apr-2011 15:07  4.2M  
            fns = re.findall('<a href=\"(.*)\">', gfs)
            # First check for .datcd file
            gfn = None
            for fn in fns:
              r = re.match("^a.*k2.*\.datcd$", fn)
              if r is not None:
                gfn = fn
                break
            if gfn is None:
              # Check for .dat file (with k2)
              for fn in fns:
                r = re.match("^a.*k2.*\.dat$", fn)
                if r is not None:
                  gfn = fn
                  break
            if gfn is None:
              # Check for .dat file (with any k)
              for fn in fns:
                r = re.match("^a.*k?.*\.dat$", fn)
                if r is not None:
                  gfn = fn
                  break
            if gfn is None:
              # Ignore this grid, there is nothing here
              continue
            self.grids[grid[1]] = (grid[0]+gfn)
            
          pickle.dump(self.grids, f)
      except:
        # Delete the half-completed file...
        self._fs.removeFile(gridFile)
        raise
        
    else:
      # File does already exist
      try:
        self.grids = pickle.load(self._fs.requestFile(gridFile, 'r', openMethod=gzip.open))
      except ValueError as ve:
        ffn = self._fs.composeFilename(gridFile)
        per = PE.PyAValError("Reading the pickle file: " + str(ffn) + ", the error: '" + str(ve) + "' occurred.", \
                             solution="Likely, the file was written using Python 3 and you try to read it using Python 2.x. Try to delete the file.")
        raise(per)
  
  def _abundToStr(self, met):
    """
      Convert metallicity into a grid-naming compatible string
      
      Parameters
      ----------
      met : float
          Log10 of metallicity
      
      Returns
      -------
      Naming string : string
          E.g., M01, P00 etc.
    """
    if met < 0:
      result = "M"
    else:
      result = "P"
    if (met * 10) - int(abs(met) * 10) > 1e-14:
      raise(PE.PyAValError("Inappropriate met value: " + str(met)))
    met = int(abs(met) * 10)
    if met < 10:
      result += "0"
    result += str(met)
    return result
  
  def _downloadGrid(self, name):
    """
      Download a model-grid file from the Kurucz page.
      
      Only downloads, if the file does not already exist
      in PyA's data path.
      
      Parameters
      ----------
      name : string
          The name of the grid.
      
      Returns
      -------
      gfn : string
          The (relative) name of the (downloaded) file.
    """
    gfn = os.path.join("pyasl", "resBased", name + ".dat.gz")
    if self._fs.fileExists(gfn):
      return gfn
    print("Downloading model data...")
    print("  Writing data to file: ", gfn, " in PyA data path.")
    self._fs.downloadToFile(self.grids[name], gfn, clobber=True, verbose=False, \
                            openMethod=gzip.open)
    return gfn
    
  def requestModelGrid(self, met, add=""):
    """
      Get a model grid.
      
      Parameters
      ----------
      met : float
          Log10 of the metallicity. For instance use: 0.0, +0.1, +0.5, or -0.5.
      
      Returns
      -------
      Model grid : KuruczMT
          The model grid enclosed in a class instance, which allows
          easy access to the models.
    """
    if not self.gridAvailable(met, add):
      raise(PE.PyAValError("No appropriate grid available", \
                           solution="Use 'availableGrids' to see what is available."))
    gfn = self._downloadGrid(self._gridName(met, add))
    ffn = self._fs.composeFilename(gfn)
    return KuruczMT(ffn)
    
  def _gridName(self, met, add):
    """
      Compose metallicity and name addition to obtain name of model grid.
      
      Parameters
      ----------
      met : float
          Log10 of metallicity.
      add : string, optional
          An additional to the name such as "ODFNEW" or "NOVER"
    """
    return "GRID" + self._abundToStr(met) + add
  
  def gridAvailable(self, met, add=""):
    """
      Check whether model grid is available.
      
      Parameters
      ----------
      met : float
          Log10 of metallicity.
      
      Returns
      -------
      Availability flag : boolen
          True, if model grid is available.
    """
    return (self._gridName(met, add) in self.grids)      
  
  def availableGrids(self):
    """
      All available model grids.
      
      Returns
      -------
      Available grids : list of strings
          The names of all available model grids.
    """
    return list(self.grids.keys())
  
  def __init__(self):
    self._fs = pp.PyAFS()
    self._getListOfModelGrids()


def getKuruczModel(teff, logg, met, nameadd=""):
  """
    Obtain a Kurucz model
    
    Parameters
    ----------
    teff : float
        Effective temperature [K]
    logg : float
        Logarithmic surface gravity [cgs]
    met : float
        Logarithmic metallicity, e.g., +0.1.
    nameadd : string, optional
        Name extension of the model grid; for instance,
        "NOVER".
    
    Returns
    -------
    Model : list of strings
        The requested model.
  """
  km = KuruczModels()
  if not km.gridAvailable(met, nameadd):
    raise(PE.PyAValError("No model grid with logarithmic metallicity of " + str(met) + \
                         " and name addition " + str(nameadd) + " available."))
  mt = km.requestModelGrid(met, nameadd)
  if not mt.modelAvailable(teff, logg, met):
    raise(PE.PyAValError("No model available in grid for parameters: teff = %6.3e, logg = %6.3e, met = %6.3e " \
                         % (teff, logg, met), \
                         solution=["Available Teffs: " + str(list(mt.availableTeffs())), \
                                   "Available Loggs: " + str(list(mt.availableLoggs()))]))
  return mt.getModel(teff, logg, met)
  
    