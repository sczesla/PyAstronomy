import numpy as np
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
