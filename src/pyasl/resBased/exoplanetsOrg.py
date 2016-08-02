from __future__ import print_function, division
from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy.pyaC import pyaErrors as PE 
import PyAstronomy.pyaC as pyaC
import os
import gzip
import csv
import numpy as np
import six.moves as smo
import six

class ExoplanetsOrg(pp.PyAUpdateCycle):
  """
    Easy access to explanets.org exoplanet archive.
    
    This class downloads a table of exoplanet data from
    the explanets.org
    (http://exoplanets.org/csv-files/)
    and provides access to these data. By default, the
    data will be re-downloaded every seven days.
    
    The following data are provided:
    
    ===========   ===================================  ======
    Column        Description                          Unit 
    -----------   -----------------------------------  ------                       
    pl_name       Name of the planet                             
    pl_orbper     Planetary orbital period             d    
    pl_massj      Planetary mass                       MJ  
    pl_msini      Minimum planetary mass               MJ
    pl_radj       Planetary radius                     RJ   
    pl_trandep    Central depth of transit             %
    pl_impact     Impact Parameter                     Stellar Radii  
    pl_trandur    Transit duration                     d    
    pl_tranmid    Transit midpoint                     BJD
    pl_tperi      Time of Periastron passage           BJD
    pl_orbsmax    Semi-major-axis                      AU 
    pl_orbsmaxr   Semi-major-axis / R_star             Stellar Radii
    pl_orbincl    Orbital inclination of planet        deg  
    pl_missal     Orbital misalignment of planet       deg
    pl_omega      Argument of Periastron               deg
    pl_ecc        Planetary orbital eccentricity       %
    pl_grav       Planetary surface gravity            log10(cm/s^2)
    pl_dens       Planetary Density                    g/cm^3
    pl_dtype      Detection type                       %
    KOI           Kepler ID (if available)
    
    pl_hostname   Name of host star
    st_binary     Binary Flag
    st_rad        Stellar radii                        Solar
    st_dist       Distance to star                     pc
    st_par        Stellar parallax                     mas
    st_mass       Stellar mass                         Solar
    st_teff       Effective temperature of star        K    
    st_vsini      Stellar vsin(i)                      km/s 
    st_logg       Stellar surface gravity              cm/s**2
    st_acts       Stellar S-Index
    st_actc       Stellar chromospheric activity                     
    st_vj         Stellar V-band brightness            mag
    st_fe         Stellar metallicity
    st_radv       System radial velocity               km/s
    st_dens       Density of star                      g/cm**3
    K             Velocity Semi-amplitude              m/s
    
    dec           Declination (J2000)                  dms  
    ra            Right ascension (J2000)              hms 
    ===========   ===================================  ======
  """
  
  def _downloadData(self):
    """
      Download data.
    """
    url = "http://exoplanets.org/csv-files/exoplanets.csv"
    self._fs.downloadToFile(url, self.dataFileName, clobber=True, verbose=False,
                            openMethod=gzip.open)
 
  def __init__(self, skipUpdate=False):
    self.data = None
    self.dataFileName = os.path.join("pyasl", "resBased", "eporg.csv.gz")
    configFilename = os.path.join("pyasl", "resBased", "eporg.cfg")
    pp.PyAUpdateCycle.__init__(self, configFilename, "eporgupdate")
    # Define columns to select
    # Column name, Description, Unit
    self._columns = {}
    # planetary parameters
    self._columns[0] = ["pl_name", "Name of the planet", "", "S15"]
    self._columns[1] = ["pl_orbper", "Planetary orbital period", "d", np.float]
    self._columns[2] = ["pl_massj", "Planetary mass", "MJ", np.float]
    self._columns[3] = ["pl_msini", "Minimum planetary mass", "MJ", np.float]
    self._columns[4] = ["pl_radj", "Planetary radius", "RJ", np.float]
    self._columns[5] = ["pl_trandep", "Central depth of transit", "(r_p/r_star)^2", np.float]
    self._columns[6] = ["pl_impact", "Impact Parameter", "Stellar Radii", np.float]
    self._columns[7] = ["pl_trandur", "Duration of transit", "d", np.float]
    self._columns[8] = ["pl_tranmid", "Transit midpoint", "BJD", np.float]
    self._columns[9] = ["pl_tperi", "Time of Periastron passage", "d", np.float] 
    self._columns[10] = ["pl_orbsmax", "Semi-major-axis", "AU", np.float]    
    self._columns[11] = ["pl_orbsmaxr", "Ratio sma to R_star", "Stellar Radii", np.float]
    self._columns[12] = ["pl_orbincl", "Orbital inclination of planet", "deg", np.float]
    self._columns[13] = ["pl_missal", "Orbital misalignment of planet ", "deg", np.float]
    self._columns[14] = ["pl_omega", "Argument of Periastron", "deg", np.float]
    self._columns[15] = ["pl_ecc", "Planetary orbital eccentricity", "", np.float]
    self._columns[16] = ["pl_grav", "Planetary surface gravity", "log10(cm/s^2)", np.float]
    self._columns[17] = ["pl_dens", "Planetary Density", "g/cm^3", np.float]
    self._columns[18] = ["pl_dtype", "Detection type", "", "S27"]
    self._columns[19] = ["KOI", "Kepler ID (if available)", "", np.float]
    # stellar pars
    self._columns[20] = ["pl_hostname","Name of host star", "", "S15"]
    self._columns[21] = ["st_binary", "Binary Flag", "", np.float]
    self._columns[22] = ["st_rad", "Radius of star", "Solar", np.float]
    self._columns[23] = ["st_dist", "Distance to Star", "pc", np.float]
    self._columns[24] = ["st_par", "Stellar parallax", "mas", np.float]
    self._columns[25] = ["st_mass", "Stellar mass", "Solar", np.float]
    self._columns[26] = ["st_teff", "Stellar effective temperature", "K", np.float]
    self._columns[27] = ["st_vsini", "Stellar vsin(i)", "km/s", np.float]
    self._columns[28] = ["st_logg", "Stellar surface gravity", "cm/s**2", np.float]    
    self._columns[29] = ["st_acts", "Stellar S-Index ", "", np.float]
    self._columns[30] = ["st_actc", "Stellar chromospheric activity", "", np.float]
    self._columns[31] = ["st_vj", "Stellar V-band brightness", "mag", np.float]
    self._columns[32] = ["st_fe", "Stellar metallicity", "", np.float]
    self._columns[33] = ["st_radv", "System radial velocity", "km/s", np.float]
    self._columns[34] = ["st_dens", "Density of star", "g/cm^3", np.float]
    self._columns[35] = ["K", "Velocity Semi-amplitude", "m/s", np.float] 
    self._columns[36] = ["dec", "Declination (J2000)", "dms", "S12"]
    self._columns[37] = ["ra", "Right ascension (J2000)", "hms", "S12"]
    
    self._ident = {"STAR":"pl_hostname", "NAME":"pl_name", "PER":"pl_orbper", "MASS":"pl_massj", "MSINI":"pl_msini", "R":"pl_radj", \
                   "DEPTH":"pl_trandep","B": "pl_impact","T14":"pl_trandur", "TT":"pl_tranmid", "T0":"pl_tperi", "A":"pl_orbsmax", "AR":"pl_orbsmaxr", \
                   "I":"pl_orbincl","LAMBDA":"pl_missal", "OM":"pl_omega", "ECC":"pl_ecc","GRAVITY": "pl_grav", "DENSITY":"pl_dens", \
                   "PLANETDISCMETH": "pl_dtype",\
                   "BINARY":"st_binary","RSTAR":"st_rad", "DIST":"st_dist", "PAR":"st_par" ,"MSTAR":"st_mass","TEFF":"st_teff","LOGG":"st_logg", "VSINI":"st_vsini", \
                   "SHK":"st_acts", "RHK":"st_actc","V":"st_vj", "FE":"st_fe", "GAMMA":"st_radv", "RHOSTAR":"st_dens",\
                   "DEC_STRING":"dec", "RA_STRING":"ra", "KOI":"KOI", "K":"K"}
    # Check whether data file exists
    self._fs = pp.PyAFS()
    if not skipUpdate:
      if self.needsUpdate() or (not self._fs.fileExists(self.dataFileName)):
        # Data needs update
        print("Downloading exoplanet data from explanets.org exoplanet archive")
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
    
    # List of keys which shall be extracted
    desiredKeys = [x for x in six.iterkeys(self._ident)]
    
    # Assigned data
    for i, x in enumerate(r):
      for k in desiredKeys:
        v = x[k]
        if len(v) == 0:
          v = None
        self.data[self._ident[k]][i] = v

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
      print("-"*63)
      for v in six.itervalues(self._columns):
        print("{0:12s}  {1:35s}  {2:5s}".format(*v))
      print("-"*63)
    
    return [self._columns[x][0] for x in smo.range(len(self._columns))]
    
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
    plnames = [name.decode("utf8") for name in self.data.pl_name]
    r = pyaC.fuzzyMatch(planetName, plnames, caseSensitive=caseSensitive, raises=True)
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

    