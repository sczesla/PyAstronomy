# -*- coding: utf-8 -*-
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE


def localTime(utc, lon, diff=True):
  """
    Computes the Local Time for a given UTC at a given geographic longitude.
    
    The local time is computed as UTC + LONGITUDE/15.
    
    Parameters
    ----------
    utc : float or array
        The time in UTC in hours.
    lon : float or array
       The geographic (East) longitude in DEGREES for which
       local time should be calculated.
    diff : boolean, optional
         If True (default), returns the difference in HOURS between
         UTC and local time.

    Returns
    -------
    Time : float or array
        Local time in HOURS (0 - 24) for given geographic
        longitude and UTC.
    Time difference : float or array
        The difference between local and UTC time in
        hours (only returned if `diff` is True)
  """

  utc = np.array(utc, ndmin=1)
  lon = np.array(lon, ndmin=1)
  if lon.size != utc.size:
    raise(PE.PyAValError("You need to specify the same number of longitudes and times", \
          solution="Make `lon` and `utc` arrays have the same length.", \
          where="localTime"))
  indi = np.where(np.logical_or(lon<0.0, lon>360.))[0]
  if len(indi) > 0:
    raise(PE.PyAValError("Longitude needs to be in the range 0-360.", \
          solution="Change the input.", \
          where="localTime"))
    
  localtime = utc + lon/15.
  localtime = localtime % 24.0

  if diff == True:
    return localtime, lon/15.
  else:
    return localtime

