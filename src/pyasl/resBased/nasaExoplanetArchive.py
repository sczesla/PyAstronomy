from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy.pyaC import pyaErrors as PE 
import PyAstronomy.pyaC as pyaC
import ConfigParser as CP
import datetime as DT
import os
import urllib2
import gzip
import csv
import numpy as np

class NasaExoplanetArchive:
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
    ===========   ===================================  ======
  """
  
  def _downloadData(self):
    """
      Download data and store it to file in PyA's data directory.
    """
    urlRoot = "http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?"
    table = "&table=exoplanets"
    select = "&select="
    for v in self._columns.itervalues():
      select = select.join([select, ',', v[0]])
    outformat = "&format=csv"
    
    response = urllib2.urlopen(''.join((urlRoot, table+select, outformat)))
    data = response.read()
    self._fs.requestFile(self.dataFileName, 'w', gzip.open).write(data)
    
    # Write download date to file
    self.config.set("NEXA", "DATA_DOWNLOAD_DATE", \
                   DT.datetime.now().strftime("%Y-%m-%d %H:%M"))
    self.config.write(self._fs.requestFile(self.configFileName, 'w'))
  
  def __init__(self):
    self.data = None
    self.dataFileName = os.path.join("pyasl", "resBased", "NEXA.csv.gz")
    self.configFileName = os.path.join("pyasl", "resBased", "NEXA.cfg")
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
    # Check whether data file exists
    self._fs = pp.PyAFS()
    self.config = CP.RawConfigParser()
    if not self._fs.fileExists(self.dataFileName):
      # It does not yet exist
      # Write config stub
      self.config.add_section('NEXA')
      # By default, reload every 7 days
      self.config.set('NEXA', 'DATA_UPDATE_CYCLE_DAYS', 7)
      self.config.write(self._fs.requestFile(self.configFileName, 'w'))
      
      print "Downloading exoplanet data from NASA exoplanet archive"
      self._downloadData()
      print "Saved data to file: ", self.dataFileName, " in data directory,"
      print "  which has been configured as: ", self._fs.dpath
      print "By default, the data will be downloaded anew every 7 days."
      print "You can use the `changeDownloadCycle` to change this behavior."
    else:
      # Data file exists, check whether it is too old
      self.config.readfp(self._fs.requestFile(self.configFileName, 'r'))
      oldTimeStr = self.config.get("NEXA", "DATA_DOWNLOAD_DATE")
      oldTime = DT.datetime.strptime(oldTimeStr, "%Y-%m-%d %H:%M")
      delta = DT.datetime.now() - oldTime
      # Convert to days
      ddays = (delta.total_seconds()/86400.)
      try:
        ddaysLimit = self.config.getfloat("NEXA", "DATA_UPDATE_CYCLE_DAYS")
      except ValueError:
        ddaysLimit = None
      if ddaysLimit is not None:
        if ddays > ddaysLimit:
          print "Download of exoplanet data. Last download was more than " \
                + str(ddaysLimit) + " days ago."
          self._downloadData()
    self._readData()

  def downloadData(self):
    """
      Trigger download of data.
    """
    self._downloadData()

  def changeDownloadCycle(self, c):
    """
      Change the time after which the data are re-downloaded.
      
      By default, the data will be downloaded if they are older
      than seven days.
      This method allows you to change that cycle.
      
      Parameters
      ----------
      c : float or None
          The new download cycle in days. If `None` is
          provided, re-downloading is switched off.
    """
    if c is not None:
      if c < 0.0:
        raise(PE.PyAValError("Update cycle needs to be a positive number of days."))
      self.config.set("NEXA", "DATA_UPDATE_CYCLE_DAYS", c)
    else:
      self.config.set("NEXA", "DATA_UPDATE_CYCLE_DAYS")
    self.config.write(self._fs.requestFile(self.configFileName, 'w'))
      
  def _readData(self):
    """
      Read the data from local file into numpy recarray.
    """
    r = csv.DictReader(self._fs.requestFile(self.dataFileName, 'r', gzip.open), delimiter=',')
    for nplanets, x in enumerate(r):
      pass
    r = csv.DictReader(self._fs.requestFile(self.dataFileName, 'r', gzip.open), delimiter=',')
    dtype = map(lambda x: (self._columns[x][0], self._columns[x][3]), range(len(self._columns)))
    self.data = np.recarray((nplanets+1,), dtype=dtype)
    for i, x in enumerate(r):
      for k, v in x.iteritems():
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
      print "{0:12s}  {1:35s}  {2:5s}".format("Column", "Description", "Unit")
      print "-"*56
      for v in self._columns.itervalues():
        print "{0:12s}  {1:35s}  {2:5s}".format(*v)
    
    return map(lambda x:self._columns[x][0], range(len(self._columns)))
    
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
    r = pyaC.fuzzyMatch(planetName, self.data.pl_name, caseSensitive=caseSensitive, raises=True)
    result = {}
    for c in self._columns.itervalues():
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

    