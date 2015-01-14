from PyAstronomy.pyaC import pyaErrors as PE
from datetime import datetime as dt
import time

def decimalYear(date):
    """
    Calculates the decimal representation of a date, e.g., 2013.12.
    
    The code uses Python's datetime package to determine the fractional
    year. Thus, leap years are taken into account. There may still be
    issues with time zones or daylight saving times etc..
    
    Code from:
    http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years

    Parameters
    ----------
    date : python date instance
        The input date (and time).
      
    Returns
    -------
    Decimal year : float
        Decimal representation of the date.
    """
    def s(date):
      # returns seconds since epoch
      return time.mktime(date.timetuple())

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


def decimalYearGregorianDate(date, form="datetime"):
  """
    Convert decimal year into gregorian date.
    
    Formally, the precision of the result is one microsecond.
    
    Parameters
    ----------
    date : float
        The input date (and time).
    form : str, optional
        Output format for the date. Either one of the following strings
        defining a format: "dd-mm-yyyy [hh:mm:ss]", "yyyy-mm-dd [hh:mm:ss]",
        where the term in braces is optional, or 
        "tuple" or "datetime". If 'tuple' is specified, the result will be a
        tuple holding (year, month, day, hour, minute, second, microseconds).
        In the case of "datetime" (default), the result will be a
        datetime object.
      
    Returns
    -------
    Gregorian date : str, tuple, or datetime instance
        The gregorian representation of the input in the specified format.
        In case of an invalid format, None is returned.
  """
  def s(date):
    # returns seconds since epoch
    return time.mktime(date.timetuple())
  
  # Shift the input of 1e-2 microseconds
  # This step accounts for rounding issues.
  date += 1e-5/(365.*86400.0)
  
  year = int(date)
  yearFraction = float(date) - int(date)
  startOfThisYear = dt(year=year, month=1, day=1)
  startOfNextYear = dt(year=year+1, month=1, day=1)
  secondsInYear = (s(startOfNextYear) - s(startOfThisYear)  ) * yearFraction
  # Find the month
  m = 1
  while m<=12 and s(dt(year=year, month=m, day=1)) - s(startOfThisYear) <= secondsInYear: m+=1
  m-=1
  # Find the day
  d = 1
  tdt = dt(year=year, month=m, day=d)
  while s(tdt) - s(startOfThisYear) <= secondsInYear: 
    d+=1
    try: tdt=dt(year=year, month=m, day=d)
    except: break
  d-=1 
  # Find the time
  secondsRemaining = secondsInYear + s(startOfThisYear) - s(dt(year=year, month=m, day=d))
  hh = int(secondsRemaining/3600.)
  mm = int((secondsRemaining - hh*3600) / 60.)
  ss = int(secondsRemaining - hh*3600 - mm * 60) 
  ff = secondsRemaining - hh*3600 - mm * 60 - ss
  
  # Output formating
  if "tuple" == form:
    r = (year, m, d, hh, mm, ss, int(ff*1000))
  elif "datetime" == form:
    r = dt(year, m, d, hh, mm, ss, int(ff*1000))
  elif "dd-mm-yyyy" in form:
    r = str("%02i-%02i-%04i" % (d,m,year))
    if "hh:mm:ss" in form:
      r+=str(" %02i:%02i:%02i" % (hh,mm,ss))
  elif "yyyy-mm-dd" in form:
    r = str("%04i-%02i-%02i" % (year,m,d))
    if "hh:mm:ss" in form:
      r+=str(" %02i:%02i:%02i" % (hh,mm,ss))
  else:
    raise(PE.PyAValError("Invalid input form of `form` parameter.", \
                         where="gregorianDate"))
    return None
  return r  