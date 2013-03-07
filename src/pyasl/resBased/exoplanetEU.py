from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy.pyaC import pyaErrors as PE
import os
import urllib2
import gzip
import csv
import numpy as np

class ExoplanetEU(pp.PyAUpdateCycle):
  """
    Provides access to exoplanet.eu data base.
    
    This class downloads the data base as a csv
    file and converts it into a numpy recarray.
    By default, the data are re-downloaded
    every 7 days.
    
    The available columns are:
    
    ============  ==============================  =====
     Column Name                     Description   Unit

          plName                  Name of planet       
          plMass                  Mass of planet     MJ
        plRadius                Radius of planet     RJ
          period                  Orbital period      d
             sma                 Semi-major axis     AU
    eccentricity            Orbital eccentricity       
     inclination             Orbital inclination    deg
     angDistance                Angular Distance arcsec
       pubStatus              Publication status       
      discovered               Year of discovery     yr
         updated             Date of data update       
           omega          Argument of Periastron    deg
           tperi             Epoch of Periastron      d
         detType                  Detection type       
       molecules      List of detected molecules       
          stName                    Name of star       
              ra         Right ascension (J2000)    hms
             dec             Declination (J2000)    dms
           mag_v      V magnitude of a host star    mag
           mag_i      I magnitude of a host star    mag
           mag_j      J magnitude of a host star    mag
           mag_h      H magnitude of a host star    mag
           mag_k      K magnitude of a host star    mag
            dist           Distance to host star     pc
              mh        Metallicity of host star    dex
          stMass                    Stellar mass  solar
        stRadius                  Radius of star  solar
             SpT      Spectral type of host star       
           stAge                     Stellar age     Ga
          stTeff   Stellar effective temperature      K
    ============  ==============================  =====
    
    Parameters
    ----------
    skipUpdate : boolean, optional
        If True, no re-download of the data
        will be initiated no matter how old
        they are.
    
  """
  
  def _download(self):
    """
      Download data.
    """
    response = urllib2.urlopen("http://exoplanet.eu/catalog/csv")
    data = response.read()
    self._fs.requestFile(self.dataFileName, 'w', gzip.open).write(data)
  
  def _readData(self):
    """
    """
    # Determine number of planets in the csv file
    r = csv.DictReader(self._fs.requestFile(self.dataFileName, 'r', gzip.open), delimiter=',')
    for nplanets, x in enumerate(r):
      pass
    # Reinitialize csv file
    r = csv.DictReader(self._fs.requestFile(self.dataFileName, 'r', gzip.open), delimiter=',')
    # Determine data types for numpy recarray from columns
    # and initialize array
    dtype = map(lambda x: (self._columns[x][0], self._columns[x][3]), range(len(self._columns)))
    self.data = np.recarray((nplanets+1,), dtype=dtype)
    for i, x in enumerate(r):
      for k, v in x.iteritems():
        # Remove hash and white spaces from column names
        k = k.strip('#')
        k = k.strip()
        # Translate csv column name into internal column name
        if k in self._ident:
          key = self._ident[k]
        else:
          key = k
        if len(v) == 0:
          v = None
        self.data[key][i] = v

  def availableColumns(self):
    """
      Show a summary of the available columns.
    """
    print "-"*51
    print "%12s  %30s  %5s" % ("Column Name", "Description", "Unit")
    print "-"*51
    for k, v in self._columns.iteritems():
      print "%12s  %30s  %5s" % tuple(v[0:3])
    print "-"*51
  
  def getAllData(self):
    """
      Provides all data as a numpy recarray.
      
      Returns
      -------
      Data : numpy recarray
          All available data. Use, e.g.,
          `availableColumns` to get an
          overview of the available
          information.
    """
    return self.data.copy()

  def forceUpdate(self):
    """
      Force a fresh download of the data.
      
      By default, the data will be updated every
      7 days. 
    """
    self._update(self._download)
  
  def __init__(self, skipUpdate=False):
    configFilename = os.path.join("pyasl", "resBased", "epeu.cfg")
    pp.PyAUpdateCycle.__init__(self, configFilename, "ExoUpdate")
    self.dataFileName = os.path.join("pyasl", "resBased", "epeu.csv.gz")
    self._fs = pp.PyAFS()
    if self.needsUpdate() and (not skipUpdate):
      self._update(self._download)

    # Define internal column names and data types
    self._columns = {}
    self._columns[0] = ["plName", "Name of planet", "", "S15"]
    self._columns[1] = ["plMass", "Mass of planet", "MJ", np.float]
    self._columns[2] = ["plRadius", "Radius of planet", "RJ", np.float]
    self._columns[3] = ["period", "Orbital period", "d", np.float]
    self._columns[4] = ["sma", "Semi-major axis", "AU", np.float]
    self._columns[5] = ["eccentricity", "Orbital eccentricity", "", np.float]
    self._columns[6] = ["inclination", "Orbital inclination", "deg", np.float]
    self._columns[7] = ["angDistance", "Angular Distance", "arcsec", np.float]
    self._columns[8] = ["pubStatus", "Publication status", "", "S2"]
    self._columns[9] = ["discovered", "Year of discovery", "yr", np.float]
    self._columns[10] = ["updated", "Date of data update", "", "S10"]
    self._columns[11] = ["omega", "Argument of Periastron", "deg", np.float]
    self._columns[12] = ["tperi", "Epoch of Periastron", "d", np.float]
    self._columns[13] = ["detType", "Detection type", "", "S2"]
    self._columns[14] = ["molecules", "List of detected molecules", "", "S10"]
    self._columns[15] = ["stName", "Name of star", "", "S15"]
    self._columns[16] = ["ra", "Right ascension (J2000)", "hms", "S12"]
    self._columns[17] = ["dec", "Declination (J2000)", "dms", "S12"]
    self._columns[18] = ["mag_v", "V magnitude of a host star", "mag", np.float]
    self._columns[19] = ["mag_i", "I magnitude of a host star", "mag", np.float]
    self._columns[20] = ["mag_j", "J magnitude of a host star", "mag", np.float]
    self._columns[21] = ["mag_h", "H magnitude of a host star", "mag", np.float]
    self._columns[22] = ["mag_k", "K magnitude of a host star", "mag", np.float]
    self._columns[23] = ["dist", "Distance to host star", "pc", np.float]
    self._columns[24] = ["mh", "Metallicity of host star", "dex", np.float]
    self._columns[25] = ["stMass", "Stellar mass", "solar", np.float]
    self._columns[26] = ["stRadius", "Radius of star", "solar", np.float]
    self._columns[27] = ["SpT", "Spectral type of host star", "", "S5"]
    self._columns[28] = ["stAge", "Stellar age", "Ga", np.float]
    self._columns[29] = ["stTeff", "Stellar effective temperature", "K", np.float]
    # Identify exoplanet.eu csv column names with internal column names
    self._ident = {"name":"plName", "mass":"plMass", "radius":"plRadius", \
                   "axis":"sma", "angular_distance":"angDistance", "publication_status":"pubStatus", \
                   "detection_type":"detType", "star.name":"stName", "star.magnitude_v":"mag_v", \
                   "star.magnitude_i":"mag_i", "star.magnitude_j":"mag_j", "star.magnitude_h":"mag_h", \
                   "star.magnitude_k":"mag_k", "star.distance":"dist", "star.metallicity":"mh", \
                   "star.mass":"stMass", "star.radius":"stRadius", "star.spec_type":"SpT", \
                   "star.age":"stAge", "star.teff":"stTeff"}
    
    self._readData()