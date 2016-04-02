from __future__ import print_function, division
from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy.pyaC import pyaErrors as PE 
import PyAstronomy.pyaC as pyaC
import os
import gzip
import csv
import numpy as np
import six

class NasaExoplanetArchive(pp.PyAUpdateCycle):
  """
    Easy access to NASA's exoplanet archive.
    
    This class downloads a table of exoplanet data from
    the "NASA Exoplanet Archive"
    (http://exoplanetarchive.ipac.caltech.edu/index.html)
    and provides access to these data. By default, the
    data will be re-downloaded every seven days.
    
    The following data are provided:
    
    ===========   ===================================  ======
    Column        Description                          Unit 
    -----------   -----------------------------------  ------
    pl_hostname   Name of host star                         
    pl_name       Name of the planet                        
    pl_letter     Planet letter (e.g., b, c, d, etc.)       
    ra            Right ascension                      deg  
    dec           Declination                          deg  
    pl_orbper     Planetary orbital period             d    
    pl_massj      Planetary mass                       MJ   
    pl_radj       Planetary radius                     RJ   
    pl_trandep    Central depth of transit             %    
    pl_trandur    Transit duration                     d    
    pl_tranmid    Transit midpoint                     BJD  
    pl_orbsmax    Semi-major-axis                      AU   
    pl_orbincl    Orbital inclination of planet        deg  
    st_rad        Stellar radii                        Solar
    st_dist       Distance to star                     pc
    st_mass       Stellar mass                         Solar
    st_teff       Effective temperature of star        K    
    st_vsini      Stellar vsin(i)                      km/s 
    st_logg       Stellar surface gravity              cm/s**2
    st_acts       Stellar S-Index                           
    st_vj         Stellar V-band brightness            mag
    ===========   ===================================  ======
  """
  
  def _downloadData(self):
    """
      Download data and store it to file in PyA's data directory.
    """
    urlRoot = "http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?"
    table = "&table=exoplanets"
    select = "&select="
    for v in six.itervalues(self._columns):
      select = ''.join([select, ',', v[0]])
    outformat = "&format=csv"
    
    url = ''.join((urlRoot, table+select, outformat))
    self._fs.downloadToFile(url, self.dataFileName, clobber=True, verbose=False, \
                            openMethod=gzip.open)
  
  def __init__(self):
    self.data = None
    self.dataFileName = os.path.join("pyasl", "resBased", "NEXA.csv.gz")
    configFileName = os.path.join("pyasl", "resBased", "NEXA.cfg")
    pp.PyAUpdateCycle.__init__(self, configFileName, "NEXA")
    # Define columns to select
    # Column name, Description, Unit
    self._columns = {}
    self._columns[0] = ["pl_hostname", "Name of host star", "", "S15"]
    self._columns[1] = ["pl_name", "Name of the planet", "", "S15"]
    self._columns[2] = ["pl_letter", "Planet letter (e.g., b, c, d, etc.)", "", "S2"]
    self._columns[3] = ["ra", "Right ascension", "deg", np.float]
    self._columns[4] = ["dec", "Declination", "deg", np.float]
    self._columns[5] = ["pl_orbper", "Planetary orbital period", "d", np.float]
    self._columns[6] = ["pl_massj", "Planetary mass", "MJ", np.float]
    self._columns[7] = ["pl_radj", "Planetary radius", "RJ", np.float]
    self._columns[8] = ["pl_trandep", "Central depth of transit", "%", np.float]
    self._columns[9] = ["pl_trandur", "Transit duration", "d", np.float]
    self._columns[10] = ["pl_tranmid", "Transit midpoint", "BJD", np.float]
    self._columns[11] = ["pl_orbsmax", "Semi-major-axis", "AU", np.float]
    self._columns[12] = ["pl_orbincl", "Orbital inclination of planet", "deg", np.float]
    self._columns[13] = ["st_rad", "Stellar radii", "Solar", np.float]
    self._columns[14] = ["st_dist", "Distance to star", "pc", np.float]
    self._columns[15] = ["st_mass", "Stellar mass", "Solar", np.float]
    self._columns[16] = ["st_teff", "Effective temperature of star", "K", np.float]
    self._columns[17] = ["st_vsini", "Stellar vsin(i)", "km/s", np.float]
    self._columns[18] = ["st_logg", "Stellar surface gravity", "cm/s**2", np.float]
    self._columns[19] = ["st_acts", "Stellar S-Index", "", np.float]
    self._columns[20] = ["st_vj", "Stellar V-band brightness", "mag", np.float]
    # Check whether data file exists
    self._fs = pp.PyAFS()
    
    if self.needsUpdate():
      # Data needs update
      print("Downloading exoplanet data from NASA exoplanet archive")
      self._update(self._downloadData)
      print("Saved data to file: ", self.dataFileName, " in data directory,")
      print("  which has been configured as: ", self._fs.dpath)
      print("By default, the data will be downloaded anew every 7 days.")
      print("You can use the `changeDownloadCycle` to change this behavior.")
    self._readData()

  def downloadData(self):
    """
      Trigger download of data.
    """
    self._update(self._downloadData)
      
  def _readData(self):
    """
      Read the data from local file into numpy recarray.
    """
    r = csv.DictReader(self._fs.requestFile(self.dataFileName, 'rt', gzip.open), delimiter=',')
    for nplanets, x in enumerate(r):
      pass
    r = csv.DictReader(self._fs.requestFile(self.dataFileName, 'rt', gzip.open), delimiter=',')
    dtype = [(self._columns[x][0], self._columns[x][3]) for x in range(len(self._columns))]
    self.data = np.recarray((nplanets+1,), dtype=dtype)
    for i, x in enumerate(r):
      for k, v in six.iteritems(x):
        if len(v) == 0:
          v = None
        self.data[k][i] = v
  
  def availableColumns(self, verbose=True):
    """
      Shows a list of available data columns.
      
      Parameters
      ----------
      verbose : boolean, optional
          If True (default), prints information
          to screen.
      
      Returns
      -------
      columns : list of strings
          The names of the available data columns.
    """
    if self.data is None:
      return None
    if verbose:
      print("{0:12s}  {1:35s}  {2:5s}".format("Column", "Description", "Unit"))
      print("-"*56)
      for v in six.itervalues(self._columns):
        print("{0:12s}  {1:35s}  {2:5s}".format(*v))
    
    return [self._columns[x][0] for x in range(len(self._columns))]
    
  def selectByPlanetName(self, planetName, caseSensitive=False):
    """
      Get entry by planet name.
      
      Parameters
      ----------
      planetName : string
          The name of the planet (includes planet letter,
          e.g., "corot-2 b"
      caseSensitive : boolean, optional
          If False (default), the search will be case-insensitive.
      
      Returns
      -------
      Data entry : dictionary
          A dictionary with a key for every data column holding
          the associated value from the data table.
    """
    pnames = [name.decode("utf8") for name in self.data.pl_name]
    r = pyaC.fuzzyMatch(planetName, pnames, caseSensitive=caseSensitive, raises=True)
    result = {}
    for c in six.itervalues(self._columns):
      result[c[0]] = self.data[c[0]][r["index"]]
    return result

  def getAllData(self):
    """
      Get all available data.
      
      Returns
      -------
      Data : numpy recarray
          All data stored in the table as a numpy recarray.
    """
    return self.data.copy()

    