from __future__ import print_function, division
from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy.pyaC import pyaErrors as PE
import os
import gzip
from PyAstronomy.pyasl import _ic
import ssl


class SWEETCat(pp.PyAUpdateCycle):
  """
    Access the SWEET-Cat catalog.

    The SWEET-Cat catalog provides parameters for planet host stars.

    The following data are provided

    ===========   ===================================  ======
    Column        Description                          Unit
    -----------   -----------------------------------  ------
    star          Name of the star
    hd            The HD number (if available)
    ra            The right ascension                   hms
    dec           The declination                       dms
    vmag          V magnitude                           mag
    ervmag        Error on V magnitude                  mag
    par           Parallax                              mas
    erpar         Error on parallax                     mas
    parsource     If par is from Simbad or calculated
    teff          Effective temperature                 K
    erteff        Error on effective temperature        K
    logg          Surface gravity                       cgs
    erlogg        Error on surface gravity              cgs
    logglc        Surface gravity from LC               cgs
    erlogglc      Error on surface gravity from LC      cgs
    vt            Micro turbulence                      km/s
    ervt          Error on micro turbulence             km/s
    metal         Metallicity ([Fe/H])
    ermetal       Error on metallicity
    mass          Mass calculated from Torres et al.    Solar
    ermass        Error on calculated mass              Solar
    author        Author of source
    link          Link to paper in ADS
    source        1 means CAUP's method. 0 otherwise
    update        When the parameters were updated
    comment1      Special comment
    comment2      Blank
    ===========   ===================================  ======


    Detailed information can be found here:
    https://www.astro.up.pt/resources/sweet-cat/
    and in the associated publications.

    Attributes
    ----------
    data : pandas data frame
        The catalog data

  """


  def _downloadData(self):
    """
      Download SWEETCAT and write it to file
    """
    url = 'https://www.astro.up.pt/resources/sweet-cat/download.php'
    dfn = self._fs.composeFilename(self.dataFileName)
    try:
      # The next two lines by-pass certificate verification
      context = ssl._create_unverified_context()
      self._fs.downloadToFile(url, dfn, clobber=True, verbose=False, openMethod=gzip.open, \
                              context=context)
    except AttributeError:
      # Python version does not support ssl._create_unverified_context()
      self._fs.downloadToFile(url, dfn, clobber=True, verbose=False, openMethod=gzip.open)
    except Exception as e:
      # "Handle" unexpected error
      raise(PE.PyANetworkError("Could not download SWEET-Cat data. The following error was raised: " + str(e)))

  def _read_sweetcat(self):
    """
      Read SWEETCat into a pandas DataFrame
    """
    if not _ic.check["pandas"]:
      raise(PE.PyARequiredImport("Could not import pandas module.", \
                                  solution="Please install pandas (http://pandas.pydata.org). Please also check the dependencies."))
    else:
      import pandas as pd
    ffn = self._fs.requestFile(self.dataFileName, 'r', gzip.open)
    self.data = pd.read_csv(ffn, sep='\t', names=self.names, na_values=['~'])

  def __init__(self, skipUpdate=False):
    self.dataFileName = os.path.join("pyasl", "resBased", "sweetcat.csv.gz")
    configFilename = os.path.join("pyasl", "resBased", "sweetcat.cfg")
    pp.PyAUpdateCycle.__init__(self, configFilename, "sweetcatupdate")

    self.names = ['star', 'hd', 'ra', 'dec', 'vmag', 'ervmag', 'par', 'erpar',
             'parsource', 'teff', 'erteff', 'logg', 'erlogg', 'logglc',
             'erlogglc', 'vt', 'ervt', 'metal', 'ermetal', 'mass', 'ermass',
             'author', 'source', 'update', 'comment']

    # Check whether data file exists
    self._fs = pp.PyAFS()
    if (self.needsUpdate() or (not self._fs.fileExists(self.dataFileName))) and (not skipUpdate) :
      # Data needs update
      print("Downloading exoplanet data from SWEET-Cat archive")
      self._update(self._downloadData)
      print("Saved data to file: ", self.dataFileName, " in data directory,")
      print("  which has been configured as: ", self._fs.dpath)
      print("By default, the data will be downloaded anew every 7 days.")
      print("You can use the `changeDownloadCycle` to change this behavior.")
    self._read_sweetcat()


  def downloadData(self):
    """
      Trigger download of data.
    """
    self._update(self._downloadData)
    self._read_sweetcat()
