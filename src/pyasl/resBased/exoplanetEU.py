from __future__ import print_function, division
from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import pyaC
from PyAstronomy.pyasl import _ic 
import os
import gzip
import csv
import numpy as np
import warnings
import six.moves as smo
import six
import copy

if _ic.check["astropy"]:
    from astropy.io import votable

class ExoplanetEU(pp.PyAUpdateCycle):
  """
    Provides access to exoplanet.eu data base.
    
    This class downloads the data base as a csv
    file and converts it into a numpy recarray.
    By default, the data are re-downloaded
    every 7 days.
    
    The available columns are:
    
    ============  ==============================  ======
     Column Name                     Description    Unit

          plName                  Name of planet       
          plMass                  Mass of planet      MJ
        plRadius                Radius of planet      RJ
          period                  Orbital period       d
             sma                 Semi-major axis      AU
    eccentricity            Orbital eccentricity       
     inclination             Orbital inclination     deg
     angDistance                Angular Distance  arcsec
       pubStatus              Publication status       
      discovered               Year of discovery      yr
         updated             Date of data update       
           omega          Argument of Periastron     deg
           tperi             Epoch of Periastron       d
         detType                  Detection type       
       molecules      List of detected molecules       
          stName                    Name of star       
              ra         Right ascension (J2000)     hms
             dec             Declination (J2000)     dms
           mag_v      V magnitude of a host star     mag
           mag_i      I magnitude of a host star     mag
           mag_j      J magnitude of a host star     mag
           mag_h      H magnitude of a host star     mag
           mag_k      K magnitude of a host star     mag
            dist           Distance to host star      pc
              mh        Metallicity of host star     dex
          stMass                    Stellar mass   solar
        stRadius                  Radius of star   solar
             SpT      Spectral type of host star       
           stAge                     Stellar age      Ga
          stTeff   Stellar effective temperature       K
        plRadMM          Measuring method of Rpl
    ============  ==============================  ======
    
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
    self._fs.downloadToFile("http://exoplanet.eu/catalog/csv", self.dataFileName, clobber=True,
                            verbose=False, openMethod=gzip.open)
  
  def _readData(self):
    """
    """
    # Determine number of planets in the csv file
    r = csv.DictReader(self._fs.requestFile(self.dataFileName, 'rt', gzip.open), delimiter=',')
    for nplanets, x in enumerate(r):
      pass
    # Reinitialize csv file
    r = csv.DictReader(self._fs.requestFile(self.dataFileName, 'rt', gzip.open), delimiter=',')
    # Determine data types for numpy recarray from columns
    # and initialize 
    dtype = [(self._columns[x][0], self._columns[x][3]) for x in range(len(self._columns))]
    self.data = np.recarray((nplanets+1,), dtype=dtype)
    colnotfilled = [self._columns[x][0] for x in six.iterkeys(self._columns)]
    for i, x in enumerate(r):
      for k, v in six.iteritems(x):
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
        # Accept only expected fields 
        if not key in self.data.dtype.names:
          continue
        try:
          colnotfilled.remove(key)
        except ValueError:
          # Ignoring already removed value
          pass
        self.data[key][i] = v
    if len(colnotfilled) > 0:
      PE.warn(PE.PyAAlgorithmFailure("Not all columns could be filled with data. The following columns must not be used: " + ", ".join(colnotfilled), \
                                     where="ExoplanetEU", \
                                     solution="The format of the data base must be checked. Please consider issuing a bug report via github."))

  def availableColumns(self):
    """
      Show a summary of the available columns.
      
      Returns
      -------
      Column names : list of strings
          The names of the columns.
    """
    print("-"*51)
    print("%12s  %30s  %5s" % ("Column Name", "Description", "Unit"))
    print("-"*51)
    cols = []
    for k, v in six.iteritems(self._columns):
      print("%12s  %30s  %5s" % tuple(v[0:3]))
      cols.append(v[0])
    print("-"*51)
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
    self._columns[30] = ["plRadMM", "Measuring method of Rpl", "", "S15"]
    # Identify exoplanet.eu csv column names with internal column names
    self._ident = {"name":"plName", "mass":"plMass", "radius":"plRadius", \
                   "semi_major_axis":"sma", "angular_distance":"angDistance", "publication_status":"pubStatus", \
                   "detection_type":"detType", "star_name":"stName", "mag_v":"mag_v", \
                   "mag_i":"mag_i", "mag_j":"mag_j", "mag_h":"mag_h", \
                   "mag_k":"mag_k", "star_distance":"dist", "star_metallicity":"mh", \
                   "star_mass":"stMass", "star_radius":"stRadius", "star_sp_type":"SpT", \
                   "star_age":"stAge", "star_teff":"stTeff", "orbital_period":"period", \
                   "radius_detection_type":"plRadMM"}
    
    self._readData()
    
    
class ExoplanetEU2(pp.PyAUpdateCycle):
  """
    Provides access to exoplanet.eu data base.
    
    This class downloads the data base as a VO table.
    By default, the data are re-downloaded
    every 7 days.
    
    The class provides access to the entire data base
    offered by exoplanet.eu as shown below. Use `showAvailableData`
    to view the table of available data, including potential
    updates.
    
    ==========================  =======  =======  ============================================================================================  
                          name    dtype     unit  description                                                                                   
    ==========================  =======  =======  ============================================================================================  
                          name   object           Name of a planet                                                                              
                          mass  float64  jovMass  Planetary Mass                                                                                
                  mass_err_min  float64  jovMass  Planetary Mass error                                                                          
                  mass_err_max  float64  jovMass  Planetary Mass error                                                                          
                     mass_sini  float64  jovMass  Planetary Mass*sin(i)                                                                         
           mass_sini_error_min  float64  jovMass  Planetary Mass*sin(i) error                                                                   
           mass_sini_error_max  float64  jovMass  Planetary Mass*sin(i) error                                                                   
                        radius  float64     Rjup  Planetary Radius                                                                              
              radius_error_min  float64     Rjup  Planetary Radius error                                                                        
              radius_error_max  float64     Rjup  Planetary Radius error                                                                        
                orbital_period  float64        d  Orbital Period                                                                                
        orbital_period_err_min  float64        d  Orbital Period error                                                                          
        orbital_period_err_max  float64        d  Orbital Period error                                                                          
               semi_major_axis  float64       AU  Semi-Major Axis                                                                               
     semi_major_axis_error_min  float64       AU  Semi-Major Axis error                                                                         
     semi_major_axis_error_max  float64       AU  Semi-Major Axis error                                                                         
                  eccentricity  float64           Orbital Eccentricity                                                                          
        eccentricity_error_min  float64           Orbital Eccentricity error                                                                    
        eccentricity_error_max  float64           Orbital Eccentricity error                                                                    
                   inclination  float64      deg  Orbital Inclination                                                                           
         inclination_error_min  float64      deg  Orbital Inclination error                                                                     
         inclination_error_max  float64      deg  Orbital Inclination error                                                                     
              angular_distance  float64     arcs  Angular Distance                                                                              
                    discovered    int32       yr  Year of Discovery                                                                             
                       updated   object           Last Update                                                                                   
                         omega  float64      deg  Argument of Periastron                                                                        
               omega_error_min  float64      deg  Argument of Periastron error                                                                  
               omega_error_max  float64      deg  Argument of Periastron error                                                                  
                         tperi  float64        d  Epoch of Periastron                                                                           
               tperi_error_min  float64        d  Epoch of Periastron error                                                                     
               tperi_error_max  float64        d  Epoch of Periastron error                                                                     
                         tconj  float64        d  Conjonction Date                                                                              
               tconj_error_min  float64        d  Conjonction Date error                                                                        
               tconj_error_max  float64        d  Conjonction Date error                                                                        
                      tzero_tr  float64        d  Primary Transit                                                                               
            tzero_tr_error_min  float64        d  Primary Transit error                                                                         
            tzero_tr_error_max  float64        d  Primary Transit error                                                                         
                  tzero_tr_sec  float64        d  Secondary Transit                                                                             
        tzero_tr_sec_error_min  float64        d  Secondary transit error                                                                       
        tzero_tr_sec_error_max  float64        d  Secondary transit error                                                                       
                  lambda_angle  float64      deg  Sky-projected angle between the planetary orbital spin and the stellar rotational spin        
        lambda_angle_error_min  float64      deg  Sky-projected angle between the planetary orbital spin and the stellar rotational spin error  
        lambda_angle_error_max  float64      deg  Sky-projected angle between the planetary orbital spin and the stellar rotational spin error  
              impact_parameter  float64           Impact Parameter b                                                                            
    impact_parameter_error_min  float64           Impact Parameter b error                                                                      
    impact_parameter_error_max  float64           Impact Parameter b error                                                                      
                      tzero_vr  float64        d  Zero Radial Speed time                                                                        
            tzero_vr_error_min  float64        d  Zero Radial Speed time error                                                                  
            tzero_vr_error_max  float64        d  Zero Radial Speed time error                                                                  
                             k  float64    m / s  Velocity Semiamplitude K                                                                      
                   k_error_min  float64    m / s  Velocity Semiamplitude K error                                                                
                   k_error_max  float64    m / s  Velocity Semiamplitude K error                                                                
               temp_calculated  float64        K  Calculated temperature                                                                        
                 temp_measured  float64        K  Measured temperature                                                                          
                 hot_point_lon  float64      deg  Hottest point longitude                                                                       
              geometric_albedo  float64           Geometric albedo                                                                              
    geometric_albedo_error_min  float64           Geometric albedo error                                                                        
    geometric_albedo_error_max  float64           Geometric albedo error                                                                        
                         log_g  float64           log(g)                                                                                        
            publication_status   object           Publication Status                                                                            
                detection_type   object           Detection type                                                                                
           mass_detection_type   object           Mass Measurement Method                                                                       
         radius_detection_type   object           Radius Measurement Method                                                                     
               alternate_names   object           List of planet alternative names                                                              
                     molecules   object           List of detected molecules                                                                    
                     star_name   object           Name of a host star                                                                           
                            ra  float64      deg  RA (J2000) of a star                                                                          
                           dec  float64      deg  Dec (J2000) of a star                                                                         
                         mag_v  float64      mag  V magnitude of a host star                                                                    
                         mag_i  float64      mag  I magnitude of a host star                                                                    
                         mag_j  float64      mag  J magnitude of a host star                                                                    
                         mag_h  float64      mag  H magnitude of a host star                                                                    
                         mag_k  float64      mag  K magnitude of a host star                                                                    
                 star_distance  float64       pc  Distance to a host star                                                                       
       star_distance_error_min  float64       pc  Distance to a host star error                                                                 
       star_distance_error_max  float64       pc  Distance to a host star error                                                                 
              star_metallicity  float64           Metallicity of a host star                                                                    
    star_metallicity_error_min  float64           Metallicity of a host star error                                                              
    star_metallicity_error_max  float64           Metallicity of a host star error                                                              
                     star_mass  float64     Msun  Mass of a host star                                                                           
           star_mass_error_min  float64     Msun  Mass of a host star error                                                                     
           star_mass_error_max  float64     Msun  Mass of a host star error                                                                     
                   star_radius  float64     Rsun  Radius of a host star                                                                         
         star_radius_error_min  float64     Rsun  Radius of a host star error                                                                   
         star_radius_error_max  float64     Rsun  Radius of a host star error                                                                   
                  star_sp_type   object           Spectral type of a host star                                                                  
                      star_age  float64      Gyr  Age of a host star                                                                            
            star_age_error_min  float64      Gyr  Age of a host star error                                                                      
            star_age_error_max  float64      Gyr  Age of a host star error                                                                      
                     star_teff  float64        K  Effective temperature of a host star                                                          
           star_teff_error_min  float64        K  Effective temperature of a host star error                                                    
           star_teff_error_max  float64        K  Effective temperature of a host star error                                                    
            star_detected_disc   object           Star Detected Disc                                                                            
           star_magnetic_field     bool           Star magnetic field                                                                           
          star_alternate_names   object           List of star alternative names                                                                
    ==========================  =======  =======  ============================================================================================  
    
    
    Parameters
    ----------
    skipUpdate : Boolean, optional
        If True, update of data will be skipped (default is False).
    forceUpdate : Boolean, optional
        If True, udpate of data will be forced (independent of value
        if `skipUpdate`, default is False)
  """
  
  
  def _readData(self):
    """
      Read data from VO table
    """
    # Temporarily suppress warnings (astropy issues one on reading this table)
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      self.vot = votable.parse(self._fs.requestFile(self.dataFileName, 'r', gzip.open))
    # Use 'name' over ID field to specify column names
    self.vot = self.vot.get_first_table().to_table(use_names_over_ids=True)
  
  def forceUpdate(self):
    """
      Force a fresh download of the data and read the data.
      
      By default, the data will be updated every
      7 days. 
    """
    self._update(self._download)
    self._readData()
  
  def showAvailableData(self):
    """
      Show available data.
      
      Displays the available column names and associated
      data types, units, and descriptions.
    """
    # Information to be used in output
    fs = ["name", "dtype", "unit", "description"]
    dat = {}
    cols = list(self.vot.columns.values())
    # Collect data
    for i in smo.range(len(cols)):
      dat[i] = {}
      for v in fs:
        dat[i][v] = cols[i].info('attributes', out=None)[v]
    # Determine maximum length of table elements 
    maxlens = {v:-1 for v in fs}
    for i in smo.range(len(cols)):
      for v in fs:
        maxlens[v] = max(maxlens[v], len(dat[i][v]))
    
    # Format the table header
    l = ""
    sep = ""
    for v in fs:
      sep += "="*maxlens[v] + "  "
      if v == "description":
        # Use left-aligned text for description
        l += ("%-" + str(maxlens[v]) + "s  ") % v
        continue
      l += ("%" + str(maxlens[v]) + "s  ") % v
    # Output header
    print(sep)
    print(l)
    print(sep)
    
    # Output information
    for i in smo.range(len(cols)):
      l = ""
      for v in fs:
        if v == "description":
          # Use left-aligned text for description
          l += ("%-" + str(maxlens[v]) + "s  ") % dat[i][v]
          continue
        l += ("%" + str(maxlens[v]) + "s  ") % dat[i][v]
      print(l)
    print(sep)
  
  def getUnitOf(self, col):
    """
      Get unit of column.
      
      Parameters
      ----------
      col : string
          The name of of column.
      
      Returns
      -------
      unit : string
          The unit used for this column.
    """
    if not col in self.vot.colnames:
      raise(PE.PyAValError("No such column: " + str(col), \
                           solution="Choose one of: " + ", ".join(self.vot.colnames)))
    return self.vot[col].info('attributes', out=None)["unit"]
  
  def _printHRInfo(self, dat):
    """
      Produce human-readable summary of available information.
      
      Parameters
      ----------
      dat : dictionary
          Output of `selectByPlanetName`.
    """
    # Search keys with/without errors
    kwe, kwoe = [], []
    for k in six.iterkeys(dat):
      if k.find("_error") != -1:
        # Ignore keys containing the '_error' phrase
        continue
      if k.find("_err_") != -1:
        # Ignore keys containing the '_error' phrase
        continue        
      if ((k + "_error_min") in dat) or ((k + "_err_min") in dat):
        # Key has error
        kwe.append(k)
      else:
        # Key has no error
        kwoe.append(k)
    
    # Maximal length of column name
    maxlen = max(list(smo.map(lambda x:len(x), self.getColnames())))
    # Get units
    units = {c:self.getUnitOf(c) for c in self.getColnames()}
    # Maximum length of unit
    mlu = max(list(smo.map(lambda x:len(x), units.values())))
    
    lines = []
    
    for k in self.getColnames():
      if (k.find("_error") != -1) or (k.find("_err_") != -1):
        # Ignore keys containing the '_error' phrase
        continue
      if k in kwoe:
        lines.append( ("%" + str(maxlen) + "s") % k + ("  [%" + str(mlu) + "s]  ") % units[k] + str(dat[k]))
      else:
        try:
          # Try _error_ to locate errors
          ep = "_error_"
          lines.append(("%" + str(maxlen) + "s") % k + ("  [%" + str(mlu) + "s]  ") % units[k] + str(dat[k]) + "(+" + str(dat[k+ep+"max"]) + ", -" + str(dat[k+ep+"min"]) + ")")
        except KeyError:
          # Try _err_ to locate errors
          ep = "_err_"
          lines.append(("%" + str(maxlen) + "s") % k + ("  [%" + str(mlu) + "s]  ") % units[k] + str(dat[k]) + "(+" + str(dat[k+ep+"max"]) + ", -" + str(dat[k+ep+"min"]) + ")")
        except:
          raise
    
    # Maximum length of output line
    mll = max(list(smo.map(lambda x:len(x), lines)))
    
    # Print to screen
    print("-"*mll)
    for l in lines :
      print(l)   
    print("-"*mll)    
    
  def selectByPlanetName(self, planetName, toScreen=True, caseSensitive=False):
    """
      Get entry by planet name.
      
      Parameters
      ----------
      planetName : string
          The name of the planet (includes planet letter,
          e.g., "corot-2 b"
      caseSensitive : boolean, optional
          If False (default), the search will be case-insensitive.
      toScreen : boolean, optional
          If True (default), the information on the system is printed
          to screen in human-readable format.
      
      Returns
      -------
      Data entry : dictionary
          A dictionary with a key for every data column holding
          the associated value from the data table.
    """
    names = [n.decode("utf8") for n in self.vot["name"]]
    r = pyaC.fuzzyMatch(planetName, names, caseSensitive=caseSensitive, raises=True)
    result = {cn:self.vot[r["index"]][cn] for cn in self.vot.colnames}
    if toScreen:
      self._printHRInfo(result)
    
    return result
  
  def getColnames(self):
    """
      Get all column names.
      
      Returns
      -------
      Columns names : list
          All column names.
    """
    return self.vot.colnames[:]
  
  def getAllDataAPT(self):
    """
      Get all data as astropy table object.
      
      Returns
      -------
      votable : astropy.table.table.Table
          All tabulated data as an astropy table
    """
    return copy.copy(self.vot)
  
  def getAllDataPandas(self):
    """
      Get all data as pandas DataFrame.
      
      Returns
      -------
      table : DataFrame
          All available data in pandas format.
    """
    if not _ic.check["pandas"]:
      raise(PE.PyARequiredImport("You need to install 'pandas' to use pandas DataFrames.", \
                                 solution="Install 'pandas' package."))
    return self.vot.to_pandas()
  
  def __init__(self, skipUpdate=False, forceUpdate=False):
    
    if not _ic.check["astropy"]:
      raise(PE.PyARequiredImport("The 'astropy' package is not installed. astropy is required to read VO tables.", \
                                 solution="Please install 'astropy'."))

    
    configFilename = os.path.join("pyasl", "resBased", "epeuvo.cfg")
    pp.PyAUpdateCycle.__init__(self, configFilename, "ExoUpdate")
    self.dataFileName = os.path.join("pyasl", "resBased", "epeu.vo.gz")
    self._fs = pp.PyAFS()
    if forceUpdate:
      self._update(self._download)
    elif (self.needsUpdate() or (not self._fs.fileExists(self.dataFileName))) and (not skipUpdate):
      # Download data if data file does not exist or
      # regular update is indicated
      self._update(self._download)
    self._readData()
      
  def _download(self):
    """
      Download data.
    """   
    self._fs.downloadToFile("http://exoplanet.eu/catalog/votable", self.dataFileName, clobber=True,
                            verbose=False, openMethod=gzip.open)
  
  