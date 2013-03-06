# -*- coding: utf-8 -*-
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
from astroTimeLegacy import helio_jd, daycnv
from PyAstronomy import constants as PC
import observatory as pyaobs
import eq2hor
import sunpos
import twilight
from moonpos import moonpos
from moonphase import moonphase
from angularDistance import getAngDist
from cardinalPoint import getCardinalPoint
import sys

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
  absPhase = np.abs((np.array(time) - T0)/period)
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
  c = PC.PyAConstants()
  # Calculate the impact parameter
  impact = (sma*c.AU)/(rs*c.RSun) * np.cos(inc/180.*np.pi)
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
        Specifies additional time before AND after the transit.
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
    if isinstance(fileOutput, basestring):
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
      planetData["dec"] = pdin["ra"]
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
      if isinstance(planetData[key], (int, long, float)):
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
      tmin = helio_jd(tmin, ra, dec)
      tmax = helio_jd(tmax, ra, dec)
  
    print "Specified time span"
    print "Start date (DDDD-MM-YY and fractional hours): {0:4d}-{1:02d}-{2:02d} {3:6.3f}".format(*daycnv(tmin))
    print "End date (DDDD-MM-YY and fractional hours): {0:4d}-{1:02d}-{2:02d} {3:6.3f}".format(*daycnv(tmax))
    print
    
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
      print "Estimating transit duration using orbital inclination, semi-major axis,"
      print "  planetary radius, and stellar radius"
    else:
      dur = planetData["Tdur"]
  
    print "Transit duration: ", dur*24.*60., " minutes"
    print "Off-transit time before and after transit: ", obsOffset*24.*60., " minutes"
  
    # First and last epoch contain in specified range
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
  
    # Check the observatory data are complete
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
    
    if moonDist is None:
      # No limit on Moon distance
      moonDist = 0.0
    if (moonDist != 0.0) and (lon is None):
      # Observer's location not given so showTwilight cannot have any effect
      raise(PE.PyAParameterConflict("The observer's location is not specified, but `moonDist` is given.\n" + \
            "This parameter can only be used, if the observer's location is known.", \
            where="transitTimes", \
            solution="Either specify the observer's location or set `moonDist` to 0.0 or None."))
  
    print
    if np.logical_and(lon != None, lat != None):
      print "No. Tmid [HJD]      Obs. start [UT] [ALT, DIR(AZI)]     Transit mid [UT] [ALT, DIR(AZI)]     Obs. end [UT] [ALT, DIR(AZI)]   twilight"+\
            " (SUN ALT)                   moon distance     moon phase"
    else:
      print "No. Tmid [HJD]      Obs. start [UT]    Transit mid [UT]   Obs. end [UT]"    
    
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
      
      # Get visibility
      if (lon is not None) and (lat is not None):
        # Get alt/az of object for current transit
        altaz = eq2hor.eq2hor(time_temp, np.ones(time_temp.size)*ra, \
                              np.ones(time_temp.size)*dec, lon=lon, lat=lat, alt=alt)      
        # If minimum altitude is not fulfilled during observation,
        # do not show transit
        if minAltitude is not None:
          minalt = np.where(altaz[0] >= minAltitude)[0]
          if len(minalt) < 3:
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
        for i in range(3):
          mdists.append(getAngDist(mpos[0][i], mpos[1][i], ra, dec))
        mdist = min(mdists)
        # Check Moon distance, if not fulfilled, neglect the transit
        if moonDist is not None:
          moonind = np.where(mdist < moonDist)[0]
          if len(moonind) > 0:
            # Neglect transit for it violates the Moon distance limit
            continue
        # Get lunar phase in percent
        moonpha = moonphase(time_temp) * 100.
        print "%3d %10.5f   %2d.%2d. %2d:%2d    [%3d°,%s(%3d°)]      %2d.%2d. %2d:%2d     [%3d°,%s(%3d°)]      %2d.%2d. %2d:%2d  [%3d°,%s(%3d°)]   %18s (%3d°,%3d°,%3d°)   (%3d°,%3d°,%3d°)  %3d%%" \
              %(trcounter, Tmid, obs_start[2], obs_start[1], np.floor(obs_start[3]), (obs_start[3]-np.floor(obs_start[3]))*60., \
                      altaz[0][0], getCardinalPoint(altaz[1][0]), altaz[1][0], \
                      obs_mid[2], obs_mid[1], np.floor(obs_mid[3]), (obs_mid[3]-np.floor(obs_mid[3]))*60., \
                      altaz[0][1], getCardinalPoint(altaz[1][1]), altaz[1][1], \
                      obs_end[2], obs_end[1], np.floor(obs_end[3]), (obs_end[3]-np.floor(obs_end[3]))*60., \
                      altaz[0][2], getCardinalPoint(altaz[1][2]), altaz[1][2], twi, sunpos_altaz[0][0], sunpos_altaz[0][1], sunpos_altaz[0][2], \
                      mdists[0], mdists[1], mdists[2], np.max(moonpha))
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
        print "%3d %10.5f   %2d.%2d. %2d:%2d       %2d.%2d. %2d:%2d       %2d.%2d. %2d:%2d" \
              %(trcounter, Tmid, obs_start[2], obs_start[1], np.floor(obs_start[3]), (obs_start[3]-np.floor(obs_start[3]))*60., \
                      obs_mid[2], obs_mid[1], np.floor(obs_mid[3]), (obs_mid[3]-np.floor(obs_mid[3]))*60., \
                      obs_end[2], obs_end[1], np.floor(obs_end[3]), (obs_end[3]-np.floor(obs_end[3]))*60. )
        trData["Tmid"] = Tmid
        trData["Obs jd"] = time_temp
        trData["Obs cal"] = [obs_start, obs_mid, obs_end]
      
      trData["Planet name"] = planetData["plName"]
      allData[trcounter] = trData
      trcounter += 1
  
    if len(allData.keys()) == 0:
      print
      print "------------------------------------------------------"
      print "!!! No transits found for the given restrictions. !!!"
      print "------------------------------------------------------"
      print
  
  except:
    raise
  finally:
    if fileOutput is not None:
      if isinstance(fileOutput, basestring):
        sys.stdout.close()
      sys.stdout = oldStdout
  
  return allData