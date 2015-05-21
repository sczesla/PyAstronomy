from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy.pyaC import pyaErrors as PE 
import os
import urllib
import gzip
from PyAstronomy.pyasl import _ic
import ssl


class SWEETCat(pp.PyAUpdateCycle):
  """
    Access the SWEET-Cat catalog.
    
    The SWEET-Cat catalog provides parameters for planet host stars.
    
    TBD
    
    Detailed information can be found here:
    https://www.astro.up.pt/resources/sweet-cat/
    and in the associated publications.
    
    Attributes
    ----------
    df : pandas data frame
        The catalog data
    
  """
  
  def _downloadData(self):
    """
      Download SWEETCAT and write it to file
    """
    url = 'https://www.astro.up.pt/resources/sweet-cat/download.php'
    dfn = self._fs.composeFilename(self.dataFileName)
    # The next two lines by-pass certificate verification
    context = ssl._create_unverified_context()
    urllib.urlretrieve(url, dfn, context=context)
    d = self._fs.requestFile(dfn, 'r').readlines()
    self._fs.requestFile(dfn, 'w', gzip.open).writelines(d)
 
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
    self.df = pd.read_csv(ffn, sep='\t', names=self.names, na_values=['~'])
    # Adding luminosity to the DataFrame
    self.df['lum'] = (self.df.teff/5777)**4 * self.df.mass 
  
  def __init__(self, skipUpdate=False):
    self.dataFileName = os.path.join("pyasl", "resBased", "sweetcat.csv.gz")
    configFilename = os.path.join("pyasl", "resBased", "sweetcat.cfg")
    pp.PyAUpdateCycle.__init__(self, configFilename, "sweetcatupdate")

    self.names = ['star', 'hd', 'ra', 'dec', 'vmag', 'ervmag', 'par', 'erpar',
             'parsource', 'teff', 'erteff', 'logg', 'erlogg', 'logglc',
             'erlogglc', 'vt', 'ervt', 'metal', 'ermetal', 'mass', 'ermass',
             'author', 'link', 'source', 'update', 'comment1', 'comment2']
                    

    # Check whether data file exists
    self._fs = pp.PyAFS()
    if self.needsUpdate():
      # Data needs update
      print "Downloading exoplanet data from SWEET-Cat archive"
      self._update(self._downloadData)
      print "Saved data to file: ", self.dataFileName, " in data directory,"
      print "  which has been configured as: ", self._fs.dpath
      print "By default, the data will be downloaded anew every 7 days."
      print "You can use the `changeDownloadCycle` to change this behavior."
    self._read_sweetcat()

  def downloadData(self):
    """
      Trigger download of data.
    """
    self._update(self._downloadData)

