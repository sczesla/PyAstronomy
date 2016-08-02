from __future__ import division
from .pyaFS import PyAFS
import six.moves.configparser as CP
import datetime as DT
from PyAstronomy.pyaC import pyaErrors as PE
import io

class PyAUpdateCycle:
  """
    A simple data update cycle.
    
    This class provides a simple way to realize an
    update cycle for data, which need to be re-initialized
    periodically.
    
    Parameters
    ----------
    fn : string
        The name of the file (within PyAFS), which
        is used to store the settings.
    section : string
        The section within the config file used
        to store the settings.
    updateCycle : float or int
        The number of days defining the update cycle.
  """
  
  def __init__(self, fn, section, updateCycle=7):
    self._config = CP.RawConfigParser()
    self._fs = PyAFS()
    self._fn = fn
    self._section = section
    if not self._fs.fileExists(fn):
      # Config does not yet exist
      # Write the default version
      self._config.add_section(section)
      self._config.set(section, 'DATA_UPDATE_CYCLE_DAYS', updateCycle)
      self._config.set(section, "DATA_DOWNLOAD_DATE", "NEVER")
      self._config.write(self._fs.requestFile(fn, 'wt'))
    else:
      # File exists
      self._config.readfp(self._fs.requestFile(self._fn, 'rt'))
      if not self._config.has_section(section):
        # Only the section does not yet exists
        # Add the content
        self._config.add_section(section)
        self._config.set(section, 'DATA_UPDATE_CYCLE_DAYS', updateCycle)
        self._config.set(section, "DATA_DOWNLOAD_DATE", "NEVER")
        self._config.write(self._fs.requestFile(fn, 'wt'))
    
    # Read the current data set from file
    self._config.readfp(self._fs.requestFile(self._fn, 'rt'))
    self._updateCycle = self._config.getfloat(self._section, "DATA_UPDATE_CYCLE_DAYS")
    self._timeString = self._config.get(self._section, "DATA_DOWNLOAD_DATE")
    
  def dataAge(self):
    """
      Determine the "age" of the data.
      
      Returns
      -------
      age : float
          The time since last data update
          in days. None, if no age can be
          determined, e.g., if data have
          never been downloaded.
    """
    if self._timeString == "NEVER":
      return None
    downloadTime = DT.datetime.strptime(self._timeString, "%Y-%m-%d %H:%M")
    delta = DT.datetime.now() - downloadTime
    # Convert into days
    ddays = (delta.total_seconds()/86400.)
    return ddays
  
  def needsUpdate(self):
    """
      Determine whether data need to be updated.
      
      Returns
      -------
      Update flag : boolean
          True if data need update and False
          otherwise.
    """
    if self._updateCycle < 0:
      # Updating switched off
      return False
    age = self.dataAge()
    if age is None:
      return True
    if age > self._updateCycle:
      return True
    return False
  
  def _setDownloadDate(self, ddate=None):
    """
      Set the download date.
      
      Parameters
      ----------
      ddate : datetime object, optional
          The time of last update. If not
          specified, the current time is
          used.
    """
    if ddate is None:
      ddate = DT.datetime.now()
    self._config.set(self._section, 'DATA_DOWNLOAD_DATE', \
                     ddate.strftime("%Y-%m-%d %H:%M"))
    self._config.write(self._fs.requestFile(self._fn, 'wt'))
    self._timeString = self._config.get(self._section, "DATA_DOWNLOAD_DATE")
  
  def changeDownloadCycle(self, c):
    """
      Change the time after which the data are updated.
      
      By default, the data will be updated if they are older
      than the given update cycle.
      This method allows you to change that cycle.
      
      Parameters
      ----------
      c : float or None
          The new update cycle in days. If `None` is
          provided, updating is switched off.
    """
    if c is not None:
      if c < 0.0:
        raise(PE.PyAValError("Update cycle needs to be a positive number of days."))
      self._updateCycle = c
    else:
      self._updateCycle = -1.0
    self._config.set(self._section, "DATA_UPDATE_CYCLE_DAYS", self._updateCycle)
    self._config.write(self._fs.requestFile(self._fn, 'wt'))

  def _update(self, func, *args, **kwargs):
    """
      Initiate data update
      
      The download date is updated automatically.
      
      Parameters
      ----------
      func : callable
          Function to call to carry out the update.
          The function will be called with *args and
          **kwargs as parameters.
    """
    func(*args, **kwargs)
    self._setDownloadDate() 