from __future__ import print_function, division
import numpy as np
from .idlMod import idlMod
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyasl import _ic


def sunpos(jd, end_jd=None, jd_steps=None, outfile=None, radian=False, plot=False, full_output=False):
  """
    Compute right ascension and declination of the Sun at a given time.
    
    Parameters
    ----------
    jd : float
         The Julian date
    end_jd : float, optional
         The end of the time period as Julian date. If given,
         `sunpos` computes RA and DEC at `jd_steps` time points
         between `jd` and ending at `end_jd`.
    jd_steps : integer, optional
         The number of steps between `jd` and `end_jd`
         for which RA and DEC are to be calculated.
    outfile : string, optional
         If given, the output will be written to a file named according
         to `outfile`.
    radian : boolean, optional
         Results are returned in radian instead of in degrees.
         Default is False.
    plot : boolean, optional
         If True, the result is plotted.
    full_output: boolean, optional
         If True, `sunpos`, additionally, returns the elongation and
         obliquity of the Sun.

    Returns
    -------
    Time : array
        The JDs for which calculations where carried out.
    Ra : array
        Right ascension of the Sun.
    Dec : array
        Declination of the Sun.
    Elongation : array, optional
        Elongation of the Sun (only of `full_output`
        is set to True).
    Obliquity : array, optional
        Obliquity of the Sun (only of `full_output`
        is set to True).

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:

    NAME:
          SUNPOS
    PURPOSE:
          To compute the RA and Dec of the Sun at a given date.
    
    CALLING SEQUENCE:
          SUNPOS, jd, ra, dec, [elong, obliquity, /RADIAN ]
    INPUTS:
          jd    - The Julian date of the day (and time), scalar or vector
                  usually double precision
    OUTPUTS:
          ra    - The right ascension of the sun at that date in DEGREES
                  double precision, same number of elements as jd
          dec   - The declination of the sun at that date in DEGREES
    
    OPTIONAL OUTPUTS:
          elong - Ecliptic longitude of the sun at that date in DEGREES.
          obliquity - the obliquity of the ecliptic, in DEGREES
    
    OPTIONAL INPUT KEYWORD:
          /RADIAN - If this keyword is set and non-zero, then all output variables 
                  are given in Radians rather than Degrees
    
    NOTES:
          Patrick Wallace (Rutherford Appleton Laboratory, UK) has tested the
          accuracy of a C adaptation of the sunpos.pro code and found the 
          following results.   From 1900-2100 SUNPOS  gave 7.3 arcsec maximum 
          error, 2.6 arcsec RMS.  Over the shorter interval 1950-2050 the figures
          were 6.4 arcsec max, 2.2 arcsec RMS.  
    
          The returned RA and Dec are in the given date's equinox.
    
          Procedure was extensively revised in May 1996, and the new calling
          sequence is incompatible with the old one.
    METHOD:
          Uses a truncated version of Newcomb's Sun.    Adapted from the IDL
          routine SUN_POS by CD Pike, which was adapted from a FORTRAN routine
          by B. Emerson (RGO).
    EXAMPLE:
          (1) Find the apparent RA and Dec of the Sun on May 1, 1982
          
          IDL> jdcnv, 1982, 5, 1,0 ,jd      ;Find Julian date jd = 2445090.5   
          IDL> sunpos, jd, ra, dec
          IDL> print,adstring(ra,dec,2)
                   02 31 32.61  +14 54 34.9
    
          The Astronomical Almanac gives 02 31 32.58 +14 54 34.9 so the error
                  in SUNPOS for this case is < 0.5".      
    
          (2) Find the apparent RA and Dec of the Sun for every day in 1997
    
          IDL> jdcnv, 1997,1,1,0, jd                ;Julian date on Jan 1, 1997
          IDL> sunpos, jd+ dindgen(365), ra, dec    ;RA and Dec for each day 
    
    MODIFICATION HISTORY:
          Written by Michael R. Greason, STX, 28 October 1988.
          Accept vector arguments, W. Landsman     April,1989
          Eliminated negative right ascensions.  MRG, Hughes STX, 6 May 1992.
          Rewritten using the 1993 Almanac.  Keywords added.  MRG, HSTX, 
                  10 February 1994.
          Major rewrite, improved accuracy, always return values in degrees
          W. Landsman  May, 1996 
          Added /RADIAN keyword,    W. Landsman       August, 1997
          Converted to IDL V5.0   W. Landsman   September 1997
  """
    
  if end_jd is None:
    # Form time in Julian centuries from 1900.0
    start_jd = (jd - 2415020.0)/36525.0
    # Zime array
    time = np.array([start_jd])
  else:
    if jd >= end_jd:
      raise(PE.PyAValError("`end_jd` needs to be larger than `jd`.", \
            where="sunpos", \
            solution="Modify the parameters."))
      if jd_steps is None:
        raise(PE.PyAValError("You specified `end_jd`, but no value for `jd_steps`.", \
              where="sunpos", \
              solution="Specify `jd_steps`, e.g., given jd_steps=10"))
    # Form time in Julian centuries from 1900.0
    start_jd = (jd - 2415020.0)/36525.0
    end_jd = (end_jd - 2415020.0)/36525.0    
    # Time array
    timestep = (end_jd-start_jd)/float(jd_steps)
    time = np.arange(start_jd, end_jd, timestep)

  # Mean solar longitude
  sunlon = (279.696678 + idlMod( (36000.768925*time), 360.0) )*3600.0

  # Allow for ellipticity of the orbit (equation of center)
  # using the Earth's mean anomaly ME
  me = 358.475844 + idlMod( (35999.049750*time) , 360.0 )
  ellcor  = ( 6910.1 - 17.2*time ) * np.sin(me*np.pi/180.) + 72.3 * np.sin(2.0*me*np.pi/180.)
  sunlon += ellcor

  # Allow for the Venus perturbations using the mean anomaly of Venus MV
  mv = 212.603219 + idlMod( (58517.803875*time) , 360.0 )
  vencorr = 4.8 * np.cos( (299.1017 + mv - me)*np.pi/180. ) + \
            5.5 * np.cos( (148.3133 +  2.0 * mv  -  2.0 * me )*np.pi/180. ) + \
            2.5 * np.cos( (315.9433 +  2.0 * mv  -  3.0 * me )*np.pi/180. ) + \
            1.6 * np.cos( (345.2533 +  3.0 * mv  -  4.0 * me )*np.pi/180. ) + \
            1.0 * np.cos( (318.15   +  3.0 * mv  -  5.0 * me )*np.pi/180. )
  sunlon += vencorr

  # Allow for the Mars perturbations using the mean anomaly of Mars MM
  mm = 319.529425  + idlMod( (19139.858500*time)  , 360.0 )
  marscorr = 2.0 * np.cos( (343.8883 -  2.0 * mm  +  2.0 * me)*np.pi/180. ) + \
            1.8 * np.cos( (200.4017 -  2.0 * mm  + me)*np.pi/180. )
  sunlon += marscorr

  # Allow for the Jupiter perturbations using the mean anomaly of Jupiter MJ
  mj = 225.328328  + idlMod( (3034.6920239*time) , 360.0 )
  jupcorr = 7.2 * np.cos( (179.5317 - mj + me )*np.pi/180. ) + \
            2.6 * np.cos( (263.2167 -  mj )*np.pi/180. ) + \
            2.7 * np.cos( ( 87.1450 -  2.0 * mj  +  2.0 * me )*np.pi/180. ) + \
            1.6 * np.cos( (109.4933 -  2.0 * mj  +  me )*np.pi/180. )
  sunlon += jupcorr

  # Allow for the Moon's perturbations using the mean elongation of
  # the Moon from the Sun D
  d = 350.7376814 + idlMod( (445267.11422*time) , 360.0 )
  mooncorr  = 6.5*np.sin(d*np.pi/180.)
  sunlon += mooncorr

  # Allow for long period terms
  longterm  = 6.4*np.sin( (231.19 + 20.20*time)*np.pi/180. )
  sunlon += longterm
  sunlon = idlMod( ( sunlon + 2592000.0 ) , 1296000.0 )
  longmed = sunlon/3600.0

  # Allow for Aberration
  sunlon -=  20.5

  # Allow for Nutation using the longitude of the Moons mean node OMEGA
  omega = 259.183275 - idlMod( (1934.142008*time) , 360.0 )
  sunlon = sunlon - 17.2*np.sin(omega*np.pi/180.)

  # Calculate the True Obliquity
  oblt = 23.452294 - 0.0130125*time + ( 9.2*np.cos(omega*np.pi/180.) )/3600.0

  # Right Ascension and Declination
  sunlon /= 3600.0
  ra = np.arctan2( np.sin(sunlon*np.pi/180.) * np.cos(oblt*np.pi/180.), np.cos(sunlon*np.pi/180.) )

  neg = np.where(ra < 0.0)[0]
  nneg = len(neg)
  if nneg > 0: ra[neg] += 2.0*np.pi

  dec = np.arcsin( np.sin(sunlon*np.pi/180.) * np.sin(oblt*np.pi/180.) )

  if radian:
    oblt *= (np.pi/180.)
    longmed *= (np.pi/180.)
  else:
    ra /= (np.pi/180.)
    dec /= (np.pi/180.)

  jd = time*36525.0 + 2415020.0

  if outfile is not None:
    # Write results to a file
    of = open(outfile, 'w')
    of.write("# File created by 'sunpos'\n")
    if not full_output:
      of.write("# 1) JD, 2) ra, 3) dec\n")
      np.savetxt(of, np.transpose(np.vstack((jd, ra, dec))))
    else:
      of.write("# 1) JD, 2) ra, 3) dec, 4) longitude, 5) obliquity\n")
      np.savetxt(of, np.transpose(np.vstack((jd, ra, dec, longmed, oblt))))
    of.close()

  if plot:
    if not _ic.check["matplotlib"]:
      raise(PE.PyARequiredImport("Could not import matplotlib.", \
            where="sunpos", \
            solution=["Install matplotlib", "Switch `plot` flag to False."]))
    import matplotlib.pylab as plt
    plt.plot(jd, ra, 'k-', label="RA")
    plt.plot(jd, dec, 'g-', label="DEC")
    plt.legend()
    plt.show()

  if full_output:
    return jd, ra, dec, longmed, oblt
  else:
    return jd, ra, dec
