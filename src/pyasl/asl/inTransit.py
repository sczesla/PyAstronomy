# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
from .astroTimeLegacy import helio_jd, daycnv
from . import observatory as pyaobs
from . import eq2hor
from . import sunpos
from . import twilight
from .moonpos import moonpos
from .moonphase import moonphase
from .angularDistance import getAngDist
from .cardinalPoint import getCardinalPoint
import sys
from . import airmass
from . import localtime
import six
import six.moves as smo

def isInTransit(time, T0, period, halfDuration, boolOutput=False):
  """
    Check whether time is inclosed by transit interval.
    
    This function uses the given ephemerides (T0, period, and
    halfDuration) to check whether the time point(s)
    specified by `time` are within a transit window or not.
    The edges of the window are counted as transit times.
    
    Parameters
    ----------
    time : float or array like
        The time(s) which are to be checked.
    T0 : float
        The time reference point (center of transit).
    period : float
        The orbital period.
    halfDuration : float
        The half-duration of the event.
        Must have same units as `time`.
    boolOutput : boolean, optional
        If set True and `time` is an array, the function will
        return a bool array holding True for time points in-
        and False for time points out-of-transit.
    
    Returns
    -------
    inTransit : boolean or array of int
        If `time` was a float, the return value is a boolean,
        which is True if the give time falls into a transit
        interval and False otherwise.
        If `time` was given as an array, the return value is
        an array holding the indices of those time points,
        which fall into a transit window. The `boolOutput`
        option may be used to obtain a boolean array holding
        True for in-transit points. 
    
  """
  if halfDuration > period/2.:
    raise(PE.PyAValError("The half-duration is longer than half the period. This cannot be true.", \
                         where="isInTransit"))
  if (period <= 0.0) or (halfDuration <= 0.0):
    raise(PE.PyAValError("Both period and half-duration must be larger 0.", \
                         where="isInTransit"))
  absPhase = np.array(np.abs((np.array(time) - T0)/period), ndmin=1)
  absPhase -= np.floor(absPhase)
  dPhase = halfDuration/period
  isIn = np.logical_or(absPhase <= dPhase, absPhase >= (1.-dPhase))
  indi = np.where(isIn)[0]
  if isinstance(time, float):
    return (len(indi) == 1)
  if boolOutput:
    return isIn
  return indi


def transitDuration(sma, rp, rs, inc, period):
  """
    Calculate the transit duration.
    
    The routine calculates the transit duration assuming
    a circular star and planet and a circular planetary orbit.
    
    In particular, it evaluates:
    
    .. math::
    
        T_D = \\frac{P}{\\pi} \\arcsin\\left(\\frac{\\sqrt{(R_s+R_p)^2 - b^2}}{a} \\right)
    
    where P is the orbital period, b the impact parameter (:math:`b=a/R_s\\cos(i)`),
    and a the semi-major axis.
    
    .. note:: The units of the transit duration are the same as the units
              of the input orbital period.
    
    Parameters
    ----------
    sma : float
        The semi-major axis in AU.
    rp : float
        The planetary radius in Jovian radii.
    rs : float
        The stellar radius in solar radii.
    inc : float
        The orbital inclination in degrees.
    period : float
        The orbital period.
    
    Returns
    -------
    Transit duration : float
        The duration of the transit (same units as
        `period`).
  """
  from PyAstronomy import constants as PC
  
  c = PC.PyAConstants()
  # Calculate the impact parameter
  impact = (sma*c.AU) * np.cos(inc/180.*np.pi)
  # Calculate the geometric transit duration
  dur = (period/np.pi) * \
        np.arcsin(np.sqrt((rs*c.RSun + rp*c.RJ)**2 - impact**2) / (sma*c.AU))
  return dur


def transitTimes(tmin, tmax, planetData, obsOffset=0., hjd=True, \
                 observatory=None, lon=None, lat=None, alt=None, minAltitude=None, \
                 showTwilight="all", moonDist=None, nexaInput=False, fileOutput=None):
  """
    Calculate transit times for a given planet and a given period of time.
    
    The `planetData` dictionary must contain the following information:
    
    =======  =================================
    Key      Value
    -------  ---------------------------------
    ra       Right ascension of object [deg]
    dec      Declination of object [deg]
    T0       Time reference point (HJD)
    orbPer   Orbital period [d]
    orbInc   Orbital inclination [deg]
    SMA      Semi-major axis [AU]
    RpJ      Planetary radius [Jovian radii]
    RsSun    Stellar Radius [solar]
    Tdur     OPTIONAL, Transit duration [d]
    =======  =================================
    
    If the transit duration (Tdur) is not given, the duration will
    be estimated using pyasl's `transitDuration` function. 
    
    .. note:: The input times (`tmin` and `tmax`) are expected in JD (UT).
              Input time will be calculated to HJD.
              Time output is in HJD.
    
    Parameters
    ----------
    tmin : float
        Start of time interval in Julian days (UT).
    tmax : float
        End of time interval in Julian days (UT).
    planetData: dictionary
        A dictionary containing the parameters of the exoplanet
        for which the transit times should be calculated.
        The required keys are specified above.
    obs_offset : float, optional
        Specifies additional time before AND after the transit in DAYS.
        This is useful if the observation should start and end
        some time before and after the actual transit.
    hjd : boolean, optional
        If True (default), the given Julian dates specifying the time
        interval (`tmin` and `tmax`) are automatically
        converted into the heliocentric frame (HJD).
    observatory : string, optional
        If given, pyasl's `observatory` function will be used to automatically
        resolve the name and obtain longitude, latitude, and altitude
        of the observatory.
        If `observatory` is given, `lon`, `lat`, and `alt` must not be specified.
    lon : float, optional
        East longitude of the observatory given in DEGREES.
         
        Longitude is positive in EASTWARD direction.
        If LON is not given, transitTimes will only return beginning and end
        of the observation and the time of mid transit.
    lat : float, optional
        Latitude of the observatory given in DEGREES
        (positive in NORTHWARD direction).
    alt : float, optional
        Altitude of the observatory given in METER.
    minAltitude : float, optional
        Minimum altitude of the object in DEGREES.
         
        If a minimum altitude is given, only transits for which the
        object is above the given altitude during the ENTIRE
        observation are included in the list
        created by `transitTimes`. Note that `minAltitude` can
        only be used if the observer's location has been specified
        either via `observatory` or via `lon`, `lat`, and `alt`.
    showTwilight : string, optional, {"all", "civil", "nautical", "astronomical", "night"}
        Specifies the twilight acceptable during the observation.
        By default all twilight conditions are acceptable.
         
        Only the transits for which the ENTIRE observation
        occurs during the specified or darker twilight conditions
        are listed.
         
        The choices are:
          - "all": all transits are shown (even during day)
          - "civil": only transits during civil twilight and better are shown
          - "nautical": only transits during nautical twilight and better are shown
          - "astronomical": only transits during astronomical twilight and better are shown
          - "night": only transits during night are shown

        Note that this can only have an effect, if the observer's location is
        specified.
    moonDist : float
        Minimum distance between the Moon and the target in DEGREES.
        By default all Moon distances are acceptable (moonDist=0.0).
         
        Only observations are listed for which the angular distance between
        the Moon and the target is larger
        than `moonDist` during the ENTIRE observation.

        Note that this can only have an effect, if the observer's location is
        specified.
    fileOutput : string or file, optional
        If a string is given, a file with the name will be created
        and the output will be written to that file. If a (writable)
        file object is given, the output will be written to that
        file. In both cases, no output will be given on screen. 

    Returns
    -------
    Transit times : dictionary
        Returns a dictionary containing the transit details. The dictionary key
        is a running number (starting with one), which is equivalent to that 
        listed in the first column of the table.
        
        For each transit, the function returns a dictionary with the transit
        details.
        
        If the observer's location was not specified, the dictionary has the
        following keys:
        
          ============    ====================================================
          Key             Description
          ------------    ----------------------------------------------------
          Planet name     Name of the planet
          Tmid            HJD of transit center
          Transit jd      Array giving JD of start, mid-time, and end of
                          transit.
          Obs jd          Array specifying the HJD of the start, center and
                          end of the observation.
          Obs cal         Equivalent to 'Obs jd', but in the form of the
                          calendar date. In particular, for each date, a list
                          containing [Year, month, day, fractional hours]
                          is given.
                          
                          **Below follows optional output only present**
                          **if the observer's location is known**
                          
          Obs coord       East longitude [deg], latitude [deg], and
                          altitude [m] of the observatory.
          Sun ra          Right ascension of the Sun at center of
                          observation.
          Sun dec         Declination of the Sun at center of
                          observation.
          Sun alt         Altitude of the Sun [deg] at begin, center, and
                          end of the observation.
          Sun az          Azimuth if the Sun [deg] at begin, center, and
                          end of the observation.
          Moon phase      Array giving lunar phase (in percent) at start,
                          center, and end of the observation.
          Moon AD         Angular distance between the target and the Moon
                          at begin, center, and end of the observation [deg].
          Moon ra         Right ascension of the Moon at begin, center, and
                          end of the observation [deg].
          Moon dec        Declination of the Moon at begin, center, and
                          end of the observation [deg].
          Star ra         Right ascension of the star [deg].
          Star dec        Declination of the star [deg].
          Star CP         Cardinal point of the star at begin, center, and
                          end of the observation.
          Star alt        Altitude of the star [deg] at begin, center, and
                          end of the observation.
          Star az         Azimuth of the star [deg] at begin, center, and
                          end of the observation.
          Twilight        The worst, i.e., brightest type of twilight
                          encountered during the observation.
          ============    ====================================================
                  
  """
    
  if fileOutput is not None:
    oldStdout = sys.stdout
    if isinstance(fileOutput, six.string_types):
      sys.stdout = open(fileOutput, 'w')
    else:
      sys.stdout = fileOutput
  
  try:
  
    if tmin >= tmax:
      raise(PE.PyAValError("The given time range is inconsistent (tmin >= tmax)", \
            where="transitTimes", \
            solution="Adapt tmin and tmax."))      
    
    # Copy input dictionary, because it may be changed
    planetData = planetData.copy()
    
    if nexaInput:
      pdin = planetData.copy()
      planetData = {}
      planetData["ra"] = pdin["ra"]
      planetData["dec"] = pdin["dec"]
      planetData["orbPer"] = pdin["pl_orbper"]
      planetData["T0"] = pdin["pl_tranmid"]
      planetData["orbInc"] = pdin["pl_orbincl"]
      planetData["SMA"] = pdin["pl_orbsmax"]
      planetData["RpJ"] = pdin["pl_radj"]
      planetData["RsSun"] = pdin["st_rad"]
      planetData["Tdur"] = pdin["pl_trandur"]
      planetData["plName"] = pdin["pl_name"]
      if np.isnan(planetData["Tdur"]):
        del planetData["Tdur"]
  
    # Check whether required keys are present
    reke = ["ra", "dec", "orbPer", "T0", "orbInc", "SMA", "RpJ", "RsSun", "plName"]
    msg = ""
    fail = False
    for key in reke:
      if (not key in planetData): 
        msg += "The required key '" + key + "' is missing in the input data!\n"
        fail = True
        continue
      if isinstance(planetData[key], (tuple(six.integer_types) + (float,))):
        if np.isnan(planetData[key]):
          msg += "The required key '" + key + "' has NaN value in the input data\n!"
          fail = True
    if fail:
      raise(PE.PyAValError("The input `planetData` is inappropriate:\n" + msg, \
                           where="transitTimes", \
                           solution="Specify all required input values."))
    
    # Object position [degrees]
    ra = planetData["ra"]
    dec = planetData["dec"]
  
    if hjd:
      # Convert input times into heliocentric frame
      # Using 'reduced' JD in calculation
      tmin = helio_jd(tmin-2.4e6, ra, dec) + 2.4e6
      tmax = helio_jd(tmax-2.4e6, ra, dec) + 2.4e6
  
    print("Specified time span")
    print("Start date (DDDD-MM-YY and fractional hours): {0:4d}-{1:02d}-{2:02d} {3:6.3f}".format(*daycnv(tmin)))
    print("End date (DDDD-MM-YY and fractional hours): {0:4d}-{1:02d}-{2:02d} {3:6.3f}".format(*daycnv(tmax)))
    print()
    
    # Transit parameters
    # Orbital period in days
    period = planetData["orbPer"]
    # Transit reference time (should be HJD)
    T0 = planetData["T0"]
    
    if not "Tdur" in planetData:
      # No duration specified in the data
      inc = planetData["orbInc"] # deg
      sma = planetData["SMA"] # au
      rp = planetData["RpJ"] # Mjup
      rs = planetData["RsSun"] # Msun
      dur = transitDuration(sma, rp, rs, inc, period)
      print("Estimating transit duration using orbital inclination, semi-major axis,")
      print("  planetary radius, and stellar radius")
    else:
      dur = planetData["Tdur"]
  
    print("Transit duration: ", dur*24.*60., " minutes")
    print("Off-transit time before and after transit: ", obsOffset*24.*60., " minutes")
  
    # First and last epoch contained in specified range
    trnum_start = np.floor( (tmin - T0)/period )
    trnum_end = np.ceil( (tmax - T0)/period )
    # Relevant transit epochs
    tr = np.arange( trnum_start, trnum_end, 1)
  
    if (observatory is not None) and \
                   ((lon is not None) or (lat is not None) or (alt is not None)):
      raise(PE.PyAParameterConflict("You must either specify `observatory` OR `lon`, `lat`, and `alt`.", \
            where="transitTimes", \
            solution="Adapt function call."))
      
    if observatory is not None:
      # Resolve observatory string
      observatory_data = pyaobs.observatory(observatory)
      lon = observatory_data["longitude"]
      lat = observatory_data["latitude"]
      alt = observatory_data["altitude"]
  
    # Check if the observatory data are complete
    obsCompl = (lon is None) + (lat is None) + (alt is None)
    if (obsCompl == 1) or (obsCompl == 2):
      raise(PE.PyAValError("Observatory data is incomplete. `lon`, `lat`, and `alt` must all be specified.\n" + \
            "Current values are: lon = " + str(lon) + ", lat = " + str(lat) + ", alt = " + str(alt), \
            where="transitTimes", \
            solution="Provide complete observatory information."))
    
    if (minAltitude is not None) and (lon is None):
      # Observer's location not given so minAltitude cannot have any effect
      raise(PE.PyAParameterConflict("The observer's location is not specified, but `minAltitude` is given.\n" + \
            "This parameter can only be used, if the observer's location is known.", \
            where="transitTimes", \
            solution="Either specify the observer's location or set `minAltitude` to None."))
  
    if (showTwilight != "all") and (lon is None):
      # Observer's location not given so showTwilight cannot have any effect
      raise(PE.PyAParameterConflict("The observer's location is not specified, but `showTwilight` is given.\n" + \
            "This parameter can only be used, if the observer's location is known.", \
            where="transitTimes", \
            solution="Either specify the observer's location or set `showTwilight` to \"all\"."))
    
    if (showTwilight != "all") and (showTwilight != "civil") and (showTwilight != "nautical") and \
       (showTwilight != "astronomical") and (showTwilight != "night"):
      # None of the possible choices for showTwilight have been used.
      raise(PE.PyAValError("Wrong keyword given for showTwilight.\n" + \
            "Current keyword is " + showTwilight, \
            where="transitTimes", \
            solution="Select a valid keyword for showTwilight: `all', `civil', `nautical', `astronomical', or `night'."))
    
    if moonDist is None:
      # No limit on Moon distance
      moonDist = 0.0
    if (moonDist != 0.0) and (lon is None):
      # Observer's location not given so showTwilight cannot have any effect
      raise(PE.PyAParameterConflict("The observer's location is not specified, but `moonDist` is given.\n" + \
            "This parameter can only be used, if the observer's location is known.", \
            where="transitTimes", \
            solution="Either specify the observer's location or set `moonDist` to 0.0 or None."))

    if moonDist < 0.0:
      # Moon distance below zero does not make sense
      PE.warn("The specified `moonDist' is below zero ("+str(moonDist)+") which does not make sense.\n"+\
              "It was changed to 0.0.\n"+\
              "Please use a value >= 0.0 or None if specifying `moonDist'.")

    print()
    if np.logical_and(lon != None, lat != None):
      print("No. Tmid [HJD]      Obs. start [UT] [ALT, DIR(AZI)]     Transit mid [UT] [ALT, DIR(AZI)]     Obs. end [UT] [ALT, DIR(AZI)]   twilight"+\
            " (SUN ALT)                   moon distance     moon phase")
    else:
      print("No. Tmid [HJD]      Obs. start [UT]    Transit mid [UT]   Obs. end [UT]")    
    
    allData = {}
    trcounter = 1
    for i in tr:
      trData = {}
      # Get times
      Tmid = T0 + float(i)*period
      obs_start_hjd = Tmid - (dur/2.0) - obsOffset
      obs_start = daycnv(obs_start_hjd)
      obs_mid = daycnv(Tmid)
      obs_end_hjd = Tmid + (dur/2.0) + obsOffset
      obs_end = daycnv(obs_end_hjd)
      time_temp = np.array([obs_start_hjd, Tmid, obs_end_hjd])
      transit_only = np.array([Tmid - (dur/2.0), Tmid, Tmid + (dur/2.0)])
      
      # Get visibility
      if (lon is not None) and (lat is not None):
        # Get alt/az of object for current transit
        altaz = eq2hor.eq2hor(time_temp, np.ones(time_temp.size)*ra, \
                              np.ones(time_temp.size)*dec, lon=lon, lat=lat, alt=alt)      
        # If minimum altitude is not fulfilled during observation,
        # do not show transit
        if minAltitude is not None:
          minalt = np.where(altaz[0] >= minAltitude)[0]
          if len(minalt) < time_temp.size:
            # Skip this transit
            continue
        # Get Sun position for current transit
        sunpos_radec = sunpos.sunpos(time_temp[1])
        sunpos_altaz = eq2hor.eq2hor(time_temp, np.ones(time_temp.size)*sunpos_radec[1], \
                                     np.ones(time_temp.size)*sunpos_radec[2], \
                                     lon=lon, lat=lat, alt=alt)
        twi = twilight.twilightName(max(sunpos_altaz[0]))
        # Check type of twilight -> if requirement not fulfilled, don't show transit
        if showTwilight == "civil":
          # Show civil or better
          if twi == "day": continue
        if showTwilight == "nautical":
          # Show nautical or better
          if twi == "day": continue
          if twi == "civil twilight": continue
        if showTwilight == "astronomical":
          # Show astronomical or better
          if twi == "day": continue
          if twi == "civil twilight": continue
          if twi == "nautical twilight": continue
        if showTwilight == "night":
          # Only show night
          if twi != "night": continue
        
        # Get Moon position for current transit
        mpos = moonpos(time_temp)
        mdists = []
        for i in range(time_temp.size):
          mdists.append(getAngDist(mpos[0][i], mpos[1][i], ra, dec))
        mdist = min(mdists)
        # Check Moon distance, if not fulfilled, neglect the transit
        if mdist < moonDist: continue
        # Get lunar phase in percent
        moonpha = moonphase(time_temp) * 100.
        print("%3d %10.5f   %2d.%2d. %2d:%02d    [%3d°,%s(%3d°)]      %2d.%2d. %2d:%02d     [%3d°,%s(%3d°)]      %2d.%2d. %2d:%02d  [%3d°,%s(%3d°)]   %18s (%3d°,%3d°,%3d°)   (%3d°,%3d°,%3d°)  %3d%%" \
              %(trcounter, Tmid, obs_start[2], obs_start[1], np.floor(obs_start[3]), (obs_start[3]-np.floor(obs_start[3]))*60., \
                      altaz[0][0], getCardinalPoint(altaz[1][0]), altaz[1][0], \
                      obs_mid[2], obs_mid[1], np.floor(obs_mid[3]), (obs_mid[3]-np.floor(obs_mid[3]))*60., \
                      altaz[0][1], getCardinalPoint(altaz[1][1]), altaz[1][1], \
                      obs_end[2], obs_end[1], np.floor(obs_end[3]), (obs_end[3]-np.floor(obs_end[3]))*60., \
                      altaz[0][2], getCardinalPoint(altaz[1][2]), altaz[1][2], twi, sunpos_altaz[0][0], sunpos_altaz[0][1], sunpos_altaz[0][2], \
                      mdists[0], mdists[1], mdists[2], np.max(moonpha)))
        # Save transit data
        trData["Tmid"] = Tmid
        trData["Obs jd"] = time_temp
        trData["Obs cal"] = [obs_start, obs_mid, obs_end]
        trData["Star ra"] = ra
        trData["Star dec"] = dec
        trData["Star alt"] = altaz[0]
        trData["Star az"] = altaz[1]
        trData["Sun ra"] = sunpos_radec[1]
        trData["Sun dec"] = sunpos_radec[2]
        trData["Sun alt"] = sunpos_altaz[0]
        trData["Sun az"] = sunpos_altaz[1]
        trData["Twilight"] = twi
        trData["Moon ra"] = mpos[0]
        trData["Moon dec"] = mpos[1]
        trData["Moon AD"] = mdist
        trData["Moon phase"] = moonpha
        trData["Star CP"] = [getCardinalPoint(altaz[1][0]), getCardinalPoint(altaz[1][1]), getCardinalPoint(altaz[1][2])]
        
        trData["Obs coord"] = [lon, lat, alt]
      else:
        # If you do not specify the observer's location, return all transits of the object  
        print("%3d %10.5f   %2d.%2d. %2d:%02d       %2d.%2d. %2d:%02d       %2d.%2d. %2d:%02d" \
              %(trcounter, Tmid, obs_start[2], obs_start[1], np.floor(obs_start[3]), (obs_start[3]-np.floor(obs_start[3]))*60., \
                      obs_mid[2], obs_mid[1], np.floor(obs_mid[3]), (obs_mid[3]-np.floor(obs_mid[3]))*60., \
                      obs_end[2], obs_end[1], np.floor(obs_end[3]), (obs_end[3]-np.floor(obs_end[3]))*60. ))
        trData["Tmid"] = Tmid
        trData["Obs jd"] = time_temp
        trData["Obs cal"] = [obs_start, obs_mid, obs_end]
      
      trData["Transit jd"] = transit_only
      trData["Planet name"] = planetData["plName"]
      allData[trcounter] = trData
      trcounter += 1
  
    if len(allData) == 0:
      print()
      print("------------------------------------------------------")
      print("!!! No transits found for the given restrictions. !!!")
      print("------------------------------------------------------")
      print()
  
  except:
    raise
  finally:
    if fileOutput is not None:
      if isinstance(fileOutput, six.string_types):
        sys.stdout.close()
      sys.stdout = oldStdout
  
  return allData


def transitVisibilityPlot(allData, markTransit=False, plotLegend=True, showMoonDist=True, print2file=False):
  """
    Plot the visibility of transits.
    
    This function can conveniently be used with the output of
    the transitTimes function.
    
    Parameters
    ----------
    allData : dictionary
        Essentially the output of `transitTimes`.
        A dictionary mapping consecutive numbers (one per transit) to
        another dictionary providing the following keys:
        
          ============    ====================================================
          Key             Description
          ------------    ----------------------------------------------------
          Planet name     Name of the planet
          Transit jd      (Only if `markTransit is True)
                          Array giving JD of start, mid-time, and end of
                          transit.
          Obs jd          Array specifying the HJD of the start, center and
                          end of the observation.
          Obs cal         Equivalent to 'Obs jd', but in the form of the
                          calendar date. In particular, for each date, a list
                          containing [Year, month, day, fractional hours]
                          is given.
          Obs coord       East longitude [deg], latitude [deg], and
                          altitude [m] of the observatory.
          Star ra         Right ascension of the star [deg].
          Star dec        Declination of the star [deg].
          ============    ====================================================

        .. note:: To use the list created by transitTimes, the LONGITUDE and LATITUDE
                  of the observatory location must have been specified.
    markTransit : boolean, optional
        If True (default is False), the in-transit times will
        be clearly indicated in the plot.
        Note that this would not be the case otherwise, which is particularly
        important if extra off-transit time before and after the transit has been
        requested. 
    showMoonDist : boolean, optional
        If True (default), the Moon distance will be shown.
    print2file : boolean, optional
        If True, the plot will be dumped to a png-file named:
        "transitVis-"[planetName].png. The default is False.
  """
  from PyAstronomy.pyasl import _ic
  if not _ic.check["matplotlib"]:
    raise(PE.PyARequiredImport("matplotlib is not installed.", \
          where="transitVisibilityPlot", \
          solution="Install matplotlib (http://matplotlib.org/)"))
  
  import matplotlib
  import matplotlib.pylab as plt
  from mpl_toolkits.axes_grid1 import host_subplot
  from matplotlib.ticker import MultipleLocator
  from matplotlib.font_manager import FontProperties
  from matplotlib import rcParams
  
  rcParams['xtick.major.pad'] = 12 
  
  if len(allData) == 0:
    raise(PE.PyAValError("Input dictionary is empty", \
          where="transitVisibilityPlot", \
          solution=["Use `transitTimes` to generate input dictionary",
                    "Did you forget to supply observer's location?", \
                    "If you used `transitTime`, you might need to change the call argument (e.g., times)"]))
  
  # Check whether all relevant data have been specified
  reqK = ["Obs jd", "Obs coord", "Star ra", "Star dec", "Obs cal", "Planet name"]
  if markTransit:
    reqK.append("Transit jd")
  missingK = []
  for k in reqK:
    if not k in allData[1]:
      missingK.append(k)
  if len(missingK) > 0:
    raise(PE.PyAValError("The following keys are missing in the input dictionary: " + ', '.join(missingK), \
                         where="transitVisibilityPlot", \
                         solution="Did you specify observer's location in `transitTimes`?"))
  
  fig = plt.figure(figsize=(15,10))
  fig.subplots_adjust(left=0.07, right=0.8, bottom=0.15, top=0.88)
  ax = host_subplot(111)

  font0 = FontProperties()
  font1 = font0.copy()
  font0.set_family('sans-serif')
  font0.set_weight('light')
  font1.set_family('sans-serif')
  font1.set_weight('medium')

  for n in six.iterkeys(allData):
    # JD array
    jdbinsize = 1.0/24./10.
    jds = np.arange(allData[n]["Obs jd"][0], allData[n]["Obs jd"][2], jdbinsize)
    # Get JD floating point
    jdsub = jds - np.floor(jds[0])
    # Get alt/az of object
    altaz = eq2hor.eq2hor(jds, np.ones(jds.size)*allData[n]["Star ra"], np.ones(jds.size)*allData[n]["Star dec"], \
                        lon=allData[n]["Obs coord"][0], lat=allData[n]["Obs coord"][1], \
                        alt=allData[n]["Obs coord"][2])
    # Get alt/az of Sun
    sunpos_altaz = eq2hor.eq2hor(jds, np.ones(jds.size)*allData[n]["Sun ra"], np.ones(jds.size)*allData[n]["Sun dec"], \
                                lon=allData[n]["Obs coord"][0], lat=allData[n]["Obs coord"][1], \
                                alt=allData[n]["Obs coord"][2])
    
    # Define plot label
    plabel = "[%02d]  %02d.%02d.%4d" % (n, allData[n]["Obs cal"][0][2], \
                                        allData[n]["Obs cal"][0][1], allData[n]["Obs cal"][0][0])
    
    # Find periods of: day, twilight, and night
    day = np.where( sunpos_altaz[0] >= 0. )[0]
    twi = np.where( np.logical_and(sunpos_altaz[0] > -18., sunpos_altaz[0] < 0.) )[0]
    night = np.where( sunpos_altaz[0] <= -18. )[0]
    
    if (len(day) == 0) and (len(twi) == 0) and (len(night) == 0):
      print()
      print("transitVisibilityPlot - no points to draw for date %2d.%2d.%4d" \
            % (allData[n]["Obs cal"][0][2], allData[n]["Obs cal"][0][1], allData[n]["Obs cal"][0][0]))
      print("Skip transit and continue with next")
      print()
      continue

    mpos = moonpos(jds)
    mpha = moonphase(jds)
    mpos_altaz = eq2hor.eq2hor(jds, mpos[0], mpos[1], lon=allData[n]["Obs coord"][0], \
                               lat=allData[n]["Obs coord"][1], alt=allData[n]["Obs coord"][2])
    moonind = np.where( mpos_altaz[0] > 0. )[0]

    if showMoonDist:
      mdist = getAngDist(mpos[0], mpos[1], np.ones(jds.size)*allData[n]["Star ra"], \
                         np.ones(jds.size)*allData[n]["Star dec"])
      bindist = int((2.0/24.)/jdbinsize)
      firstbin = np.random.randint(0,bindist)
      for mp in range(0, int(len(jds)/bindist)):
        bind = firstbin+float(mp)*bindist
        ax.text(jdsub[bind], altaz[0][bind]-1., str(int(mdist[bind]))+r"$^\circ$", ha="center", va="top", \
                fontsize=8, stretch='ultra-condensed', fontproperties=font0, alpha=1.)

    if markTransit:
      # Mark points within transit. These may differ from that pertaining to the
      # observation if an extra offset was given to provide off-transit time.
      transit_only_ind = np.where( np.logical_and(jds >= allData[n]["Transit jd"][0], \
                                                  jds <= allData[n]["Transit jd"][2]) )[0]
      ax.plot( jdsub[transit_only_ind], altaz[0][transit_only_ind], 'g', linewidth=6, alpha=.3)

    if len(twi) > 1:
      # There are points in twilight
      linebreak = np.where( (jdsub[twi][1:]-jdsub[twi][:-1]) > 2.0*jdbinsize)[0]
      if len(linebreak) > 0:
        plotrjd = np.insert(jdsub[twi], linebreak+1, np.nan)
        plotdat = np.insert(altaz[0][twi], linebreak+1, np.nan)
        ax.plot( plotrjd, plotdat, "-", color='#BEBEBE', linewidth=1.5)
      else:
        ax.plot( jdsub[twi], altaz[0][twi], "-", color='#BEBEBE', linewidth=1.5)

    ax.plot( jdsub[night], altaz[0][night], 'k', linewidth=1.5, label=plabel)
    ax.plot( jdsub[day], altaz[0][day], color='#FDB813', linewidth=1.5)
    
    altmax = np.argmax(altaz[0])
    ax.text( jdsub[altmax], altaz[0][altmax], str(n), color="b", fontsize=14, \
             fontproperties=font1, va="bottom", ha="center")

    if n == 29:
      ax.text( 1.1, 1.0-float(n)*0.04, "too many transits", ha="left", va="top", transform=ax.transAxes, \
              fontsize=10, fontproperties=font0, color="r")      
    else:
      ax.text( 1.1, 1.0-float(n)*0.04, plabel, ha="left", va="top", transform=ax.transAxes, \
              fontsize=12, fontproperties=font0, color="b")

  ax.text( 1.1, 1.03, "Start of observation", ha="left", va="top", transform=ax.transAxes, \
          fontsize=12, fontproperties=font0, color="b")
  ax.text( 1.1, 1.0, "[No.]  Date", ha="left", va="top", transform=ax.transAxes, \
          fontsize=12, fontproperties=font0, color="b")
  
  axrange = ax.get_xlim()
  ax.set_xlabel("UT [hours]")
  
  if axrange[1]-axrange[0] <= 1.0:
    jdhours = np.arange(0,3,1.0/24.)
    utchours = (np.arange(0,72,dtype=int)+12)%24
  else:
    jdhours = np.arange(0,3,1.0/12.)
    utchours = (np.arange(0,72, 2, dtype=int)+12)%24
  ax.set_xticks(jdhours)
  ax.set_xlim(axrange)
  ax.set_xticklabels(utchours, fontsize=18)  
  
  # Make ax2 responsible for "top" axis and "right" axis
  ax2 = ax.twin()
  # Set upper x ticks
  ax2.set_xticks(jdhours)
  ax2.set_xticklabels(utchours, fontsize=18)
  ax2.set_xlabel("UT [hours]")

  # Horizon angle for airmass
  airmass_ang = np.arange(5.,90.,5.)
  geo_airmass = airmass.airmassPP(90.-airmass_ang) 
  ax2.set_yticks(airmass_ang)
  airmassformat = []
  for t in range(geo_airmass.size):
    airmassformat.append("%2.2f" % geo_airmass[t])
  ax2.set_yticklabels(airmassformat, rotation=90)
  ax2.set_ylabel("Relative airmass", labelpad=32)
  ax2.tick_params(axis="y", pad=10, labelsize=10)
  plt.text(1.015,-0.04, "Plane-parallel", transform=ax.transAxes, ha='left', \
           va='top', fontsize=10, rotation=90)

  ax22 = ax.twin()
  ax22.set_xticklabels([])  
  ax22.set_frame_on(True)
  ax22.patch.set_visible(False)
  ax22.yaxis.set_ticks_position('right')
  ax22.yaxis.set_label_position('right')
  ax22.spines['right'].set_position(('outward', 25))
  ax22.spines['right'].set_color('k')
  ax22.spines['right'].set_visible(True)
  airmass2 = np.array([airmass.airmassSpherical(90. - ang, allData[n]["Obs coord"][2]) for ang in airmass_ang])
  ax22.set_yticks(airmass_ang)
  airmassformat = []
  for t in range(airmass2.size): airmassformat.append("%2.2f" % airmass2[t])
  ax22.set_yticklabels(airmassformat, rotation=90)
  ax22.tick_params(axis="y", pad=10, labelsize=10)
  plt.text(1.045,-0.04, "Spherical+Alt", transform=ax.transAxes, ha='left', va='top', \
           fontsize=10, rotation=90)

  ax3 = ax.twiny()
  ax3.set_frame_on(True)
  ax3.patch.set_visible(False)
  ax3.xaxis.set_ticks_position('bottom')
  ax3.xaxis.set_label_position('bottom')
  ax3.spines['bottom'].set_position(('outward', 50))
  ax3.spines['bottom'].set_color('k')
  ax3.spines['bottom'].set_visible(True)

  ltime, ldiff = localtime.localTime(utchours, np.repeat(allData[n]["Obs coord"][0], len(utchours)))
  jdltime = jdhours - ldiff/24.
  ax3.set_xticks(jdltime)
  ax3.set_xticklabels(utchours)
  ax3.set_xlim([axrange[0],axrange[1]])
  ax3.set_xlabel("Local time [hours]")

  ax.yaxis.set_major_locator(MultipleLocator(15))
  ax.yaxis.set_minor_locator(MultipleLocator(5))
  yticks = ax.get_yticks()
  ytickformat = []
  for t in range(yticks.size): ytickformat.append(str(int(yticks[t]))+r"$^\circ$")
  ax.set_yticklabels(ytickformat, fontsize=20)
  ax.set_ylabel("Altitude", fontsize=18)
  yticksminor = ax.get_yticks(minor=True)
  ymind = np.where( yticksminor % 15. != 0. )[0]
  yticksminor = yticksminor[ymind]
  ax.set_yticks(yticksminor, minor=True)
  m_ytickformat = []
  for t in range(yticksminor.size): m_ytickformat.append(str(int(yticksminor[t]))+r"$^\circ$")
  ax.set_yticklabels(m_ytickformat, minor=True)
  
  ax.yaxis.grid(color='gray', linestyle='dashed')
  ax.yaxis.grid(color='gray', which="minor", linestyle='dotted')
  ax2.xaxis.grid(color='gray', linestyle='dotted')

  plt.text(0.5,0.95,"Transit visibility of "+allData[n]["Planet name"].decode("utf8"), \
           transform=fig.transFigure, ha='center', va='bottom', fontsize=20)

  if plotLegend:
    line1 = matplotlib.lines.Line2D((0,0),(1,1), color='#FDB813', linestyle="-", linewidth=2)
    line2 = matplotlib.lines.Line2D((0,0),(1,1), color='#BEBEBE', linestyle="-", linewidth=2)
    line3 = matplotlib.lines.Line2D((0,0),(1,1), color='k', linestyle="-", linewidth=2)
    line4 = matplotlib.lines.Line2D((0,0),(1,1), color='g', linestyle="-", linewidth=6, alpha=.3)

    if markTransit:
      lgd2 = plt.legend((line1,line2,line3, line4),("day","twilight","night","transit",), \
                        bbox_to_anchor=(0.88, 0.15), loc=2, borderaxespad=0.,prop={'size':12}, fancybox=True)
    else:
      lgd2 = plt.legend((line1,line2,line3),("day","twilight","night",), \
                        bbox_to_anchor=(0.88, 0.13), loc=2, borderaxespad=0.,prop={'size':12}, fancybox=True)
    lgd2.get_frame().set_alpha(.5)

  targetco = r"Target coordinates: (%8.4f$^\circ$, %8.4f$^\circ$)" % \
            (allData[n]["Star ra"], allData[n]["Star dec"])
  obsco = "Obs coord.: (%8.4f$^\circ$, %8.4f$^\circ$, %4d m)" % \
          (allData[n]["Obs coord"][0], allData[n]["Obs coord"][1], allData[n]["Obs coord"][2])
  plt.text(0.01,0.97, targetco, transform=fig.transFigure, ha='left', va='center', fontsize=10)
  plt.text(0.01,0.95, obsco, transform=fig.transFigure, ha='left', va='center', fontsize=10)
  
  if print2file:
    outfile = "transVis-"+allData[n]["Planet name"].replace(" ", "")+".png"
    plt.savefig(outfile, format="png", dpi=300)
  else:
    plt.show()