from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy.pyaC import pyaErrors as PE
import os
import urllib2
import gzip
import numpy as np
import pandas as pd


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
    # Read the data structure into a pandas DataFrame
    self.data = pd.read_csv(self._fs.requestFile(self.dataFileName, 'r', gzip.open))
    # Remove # from the keys, and trailing whitespaces
    for key in self.data.keys():
        if '#' in key:
            self.data.rename(columns={key: key.replace('#', '')}, inplace=True)
            key = key.replace('#', '')
        if ' ' in key:
            # Replace whitespaces inside names with underscore, otherwise strip it
            if not key.startswith(' ') or not key.endswith(' '):
                self.data.rename(columns={key: key.replace(' ', '_')}, inplace=True)
            else:
                self.data.rename(columns={key: key.strip(' ')}, inplace=True)


  def availableColumns(self):
    """
      Show a summary of the available columns.

      Returns
      -------
      Column names : list of strings
          The names of the columns.
    """
    print "-"*51
    print "%12s  %30s  %5s" % ("Column Name", "Description", "Unit")
    print "-"*51
    cols = []
    for k, v in self._columns.iteritems():
      print "%12s  %30s  %5s" % tuple(v[0:3])
      cols.append(v[0])
    print "-"*51
    return cols


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
                   "semi_major_axis":"sma", "angular_distance":"angDistance", "publication_status":"pubStatus", \
                   "detection_type":"detType", "star_name":"stName", "mag_v":"mag_v", \
                   "mag_i":"mag_i", "mag_j":"mag_j", "mag_h":"mag_h", \
                   "mag_k":"mag_k", "star_distance":"dist", "star_metallicity":"mh", \
                   "star_mass":"stMass", "star_radius":"stRadius", "star_sp_type":"SpT", \
                   "star_age":"stAge", "star_teff":"stTeff"}

    self._readData()
