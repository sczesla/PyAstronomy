# -*- coding: utf-8 -*-
from __future__ import print_function, division
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo
import six
from six.moves import configparser as ConfigParser
import os

def observatory(obsName):
  """
    Get the location of an observatory.
    
    This function emulates the behavior of
    the IDL Astrolib's `observatory.pro`
    procedure. Use `showObservatories`
    to get a list of known observatories.
    
    Parameters
    ----------
    obsName : string
        The observatory code. Use `listObservatories`
        to get a list of available observatories
        and their codes.
    
    Returns
    -------
    Observatory location : dictionary
        The keys of this dictionary are: "longitude",
        "latitude", "altitude", and "tz". Longitude
        and latitude are given in degrees and the
        observatory's altitude in meter. Finally, tz
        defines the time zone (hours West of Greenwich).
  """
  config = ConfigParser.RawConfigParser()
  config.read(os.path.join(os.path.dirname(__file__), 'observatory.cfg'))
  
  
  if not config.has_section(obsName):
    raise(PE.PyAValError("No such observatory: '" + str(obsName) + "'", \
          solution="Choose existing observatory (use listObservatories" + \
          "to get a list)", \
          where="observatory"))
  
  result = {}
  result["name"] = config.get(obsName, "NAME")
  result["longitude"] = config.getfloat(obsName, "LONGITUDE")
  result["latitude"] = config.getfloat(obsName, "LATITUDE")
  result["altitude"] = config.getfloat(obsName, "ALTITUDE")
  result["tz"] = config.getfloat(obsName, "TZ")
  
  return result


def listObservatories(show=True):
  """
    Get a list of available observatories.
    
    Parameters
    ----------
    show : boolean, optional
        If True (default), the observatory data
        will be written to screen.
    
    Returns
    -------
    Observatory data : dictionary
        For every observatory code (key), the dictionary
        holds another dictionary with "longitude",
        "latitude", "altitude", "name", and "tz". Longitude
        and latitude are given in degrees and the
        observatory's altitude in meter. Finally, tz
        defines the time zone (hours West of Greenwich)
        and "name" gives the full observatory name.
  """
  config = ConfigParser.RawConfigParser()
  config.read(os.path.join(os.path.dirname(__file__), 'observatory.cfg'))
  codes = config.sections()
  
  result = {}
  for code in codes:
    result[code] = observatory(code)
  
  if show:
    print("List of available observatories")
    print()
    maxCodeLen = max([len(k) for k in six.iterkeys(result)])
    print(("{0:"+str(maxCodeLen)+"s}     ").format("Code") + "Observatory name")
    print("-" * (21+maxCodeLen))
    for k in sorted(list(result), key=lambda s: s.lower()):
      print(("{0:"+str(maxCodeLen)+"s} --- ").format(k) + result[k]["name"])
    print()
    print("Observatory location")
    print("--------------------")
    for k in sorted(list(result), key=lambda s: s.lower()):
      print(("Code: {code:"+str(maxCodeLen)+"s}, " + \
             "Longitude = {longitude:7.3f}°, Latitude = {latitude:+7.3f}°, " + \
            "Altitude = {altitude:5.0f} m, TZ = {tz: 5.1f}").format(code=k, **result[k]))
  return result 