# -*- coding: utf-8 -*-
from __future__ import print_function, division
import datetime
import numpy
from numpy import sin, cos, tan, sqrt, arcsin
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo


def daycnv(xjd, mode="idl"):
  """
    Converts Julian dates to Gregorian calendar dates. 
    
    Handles both individual floats as xjd and iterables such as
    lists and arrays. In the latter case, the result is returned
    in the form of a list.
    
    Parameters
    ----------
    xjd : float, list, array
        The Julian date
    mode : string, {idl, dtlist, dt}, optional
        Determines format of output. If 'idl' is given (default),
        a list holding [year, month, day, (fractional) hours] is
        returned; this mimics the behavior of the IDL astrolib function.
        If 'dtlist' is given, a list holding
        [year, month, day, hours, minutes, seconds, microseconds] is
        returned. Finally, if 'dt' is specified, a Python
        datetime object will be returned. If the input is an iterable,
        the mode determines the format of the individual items in the
        result list.
    
    Returns
    -------
    Calendar date : list or datetime object
        A list holding [year, month, day, (fractional) hours] (default)
        or [year, month, day, hours, minutes, seconds, microseconds].
        Alternatively, a Python datetime object is returned. The format
        depends on the 'mode' specified. If the input is an iterable of
        Julian dates, the output is a list.
        
    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.
    :IDL - Documentation:
    NAME:
          DAYCNV
    PURPOSE:
          Converts Julian dates to Gregorian calendar dates
    CALLING SEQUENCE:
          DAYCNV, XJD, YR, MN, DAY, HR
    INPUTS:
          XJD = Julian date, positive double precision scalar or vector
    OUTPUTS:
          YR = Year (Integer)
          MN = Month (Integer)
          DAY = Day (Integer)
          HR = Hours and fractional hours (Real).   If XJD is a vector,
                  then YR,MN,DAY and HR will be vectors of the same length.
    EXAMPLE:
          IDL> DAYCNV, 2440000.D, yr, mn, day, hr    
          yields yr = 1968, mn =5, day = 23, hr =12.   
    WARNING:
          Be sure that the Julian date is specified as double precision to
          maintain accuracy at the fractional hour level.
    METHOD:
          Uses the algorithm of Fliegel and Van Flandern (1968) as reported in
          the "Explanatory Supplement to the Astronomical Almanac" (1992), p. 604
          Works for all Gregorian calendar dates with XJD > 0, i.e., dates after
          -4713 November 23.
    REVISION HISTORY:
          Converted to IDL from Yeoman's Comet Ephemeris Generator, 
          B. Pfarr, STX, 6/16/88
          Converted to IDL V5.0   W. Landsman   September 1997
  """
  
  if not mode in ('idl', 'dtlist', 'dt'):
    raise(PE.PyAValError("Unknown mode: " + str(mode), \
                         where="daycnv", \
                         solution="Use any of 'idl', 'dtlist', or 'dt'."))

  # Adjustment needed because Julian day starts at noon, calendar day at midnight
  
  iterable = hasattr(xjd, "__iter__")
  
  # Use iterable throughout calculations
  if not iterable:
    xjd = [xjd]
  
  jd = numpy.array(xjd).astype(int)                         #Truncate to integral day
  frac = numpy.array(xjd).astype(float) - jd + 0.5          #Fractional part of calendar day
  gi = numpy.where(frac >= 1.0)
  frac[gi] -= 1.0
  jd[gi] += 1
   
  hr = frac*24.0
  l = jd + 68569
  n = 4*l // 146097
  l = l - (146097*n + 3) // 4
  yr = 4000*(l+1) // 1461001
  l = l - 1461*yr // 4 + 31        #1461 = 365.25 * 4
  mn = 80*l // 2447
  day = l - 2447*mn // 80
  l = mn//11
  mn = mn + 2 - 12*l
  yr = 100*(n-49) + yr + l
  if mode in ('dt', 'dtlist'):
    # [year, month, day, hours, minutes, seconds, microseconds] requested
    hour = numpy.floor(hr).astype(int)
    minute = numpy.floor((hr - numpy.floor(hr))*60).astype(int)
    sec = numpy.floor((hr - hour - minute/60.)*3600.).astype(int)
    msec = (3600*1e6*(hr - hour - minute/60. - sec/3600.)).astype(int)
    if mode == 'dtlist':
      if not iterable:
        return [yr[0], mn[0], day[0], hour[0], minute[0], sec[0], msec[0]]      
      return [[yr[i], mn[i], day[i], hour[i], minute[i], sec[i], msec[i]] for i in smo.range(len(yr))]
    # Return datetime object
    dts = [datetime.datetime(*(yr[i], mn[i], day[i], hour[i], minute[i], sec[i], msec[i])) for i in smo.range(len(yr))]
    if not iterable: 
      return dts[0]
    return dts
  if not iterable:
    return [yr[0], mn[0], day[0], hr[0]]
  return [[yr[i], mn[i], day[i], hr[i]] for i in smo.range(len(yr))]





def bprecess(ra, dec, mu_radec = None,  \
                  parallax = 0.0,  rad_vel = 0.0, epoch = 2000.0):
  """
    Precess positions from J2000.0 (FK5) to B1950.0 (FK4).
  
    Parameters
    ----------
    ra : float
        Right ascension [deg]
    dec : float
        Declination [deg]
    mu_radec : list
        List of two float entries, optional
        Proper motion [arcsec per tropical CENTURY]
    parallax : float
        The parallax of the target
    rad_vel : float
        Radial velocity [km/s]
        
    Returns
    -------
    Precessed position : list
        [ra_1950, dec_1950, MU_RADEC, PARALLAX, RAD_VEL] referring to 1950

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:

    NAME:
          BPRECESS
    PURPOSE:
          Precess positions from J2000.0 (FK5) to B1950.0 (FK4)
    EXPLANATION:
          Calculates the mean place of a star at B1950.0 on the FK4 system from
          the mean place at J2000.0 on the FK5 system.    

    CALLING SEQUENCE:
          bprecess, ra, dec, ra_1950, dec_1950, [ MU_RADEC = , PARALLAX = 
                                          RAD_VEL =, EPOCH =   ]

    INPUTS:
          RA,DEC - Input J2000 right ascension and declination in *degrees*.
                  Scalar or N element vector

    OUTPUTS:
          RA_1950, DEC_1950 - The corresponding B1950 right ascension and 
                  declination in *degrees*.    Same number of elements as
                  RA,DEC but always double precision.

    OPTIONAL INPUT-OUTPUT KEYWORDS
          MU_RADEC - 2xN element double precision vector containing the proper 
                    motion in seconds of arc per tropical *century* in right 
                    ascension and declination.
          PARALLAX - N_element vector giving stellar parallax (seconds of arc)
          RAD_VEL  - N_element vector giving radial velocity in km/s

          The values of MU_RADEC, PARALLAX, and RADVEL will all be modified
          upon output to contain the values of these quantities in the
          B1950 system.  The parallax and radial velocity will have a very 
          minor influence on the B1950 position.   

          EPOCH - scalar giving epoch of original observations, default 2000.0d
              This keyword value is only used if the MU_RADEC keyword is not set.
    NOTES:
          The algorithm is taken from the Explanatory Supplement to the 
          Astronomical Almanac 1992, page 186.
          Also see Aoki et al (1983), A&A, 128,263

          BPRECESS distinguishes between the following two cases:
          (1) The proper motion is known and non-zero
          (2) the proper motion is unknown or known to be exactly zero (i.e.
                  extragalactic radio sources).   In this case, the reverse of 
                  the algorithm in Appendix 2 of Aoki et al. (1983) is used to 
                  ensure that the output proper motion is  exactly zero. Better 
                  precision can be achieved in this case by inputting the EPOCH 
                  of the original observations.

          The error in using the IDL procedure PRECESS for converting between
          B1950 and J1950 can be up to 12", mainly in right ascension.   If
          better accuracy than this is needed then BPRECESS should be used.

          An unsystematic comparison of BPRECESS with the IPAC precession 
          routine (http://nedwww.ipac.caltech.edu/forms/calculator.html) always 
          gives differences less than 0.15".
    EXAMPLE:
          The SAO2000 catalogue gives the J2000 position and proper motion for
          the star HD 119288.   Find the B1950 position. 

          RA(2000) = 13h 42m 12.740s      Dec(2000) = 8d 23' 17.69''  
          Mu(RA) = -.0257 s/yr      Mu(Dec) = -.090 ''/yr

          IDL> mu_radec = 100D* [ -15D*.0257, -0.090 ]
          IDL> ra = ten(13, 42, 12.740)*15.D 
          IDL> dec = ten(8, 23, 17.69)
          IDL> bprecess, ra, dec, ra1950, dec1950, mu_radec = mu_radec
          IDL> print, adstring(ra1950, dec1950,2)
                  ===> 13h 39m 44.526s    +08d 38' 28.63"

    REVISION HISTORY:
          Written,    W. Landsman                October, 1992
          Vectorized, W. Landsman                February, 1994
          Treat case where proper motion not known or exactly zero  November 1994
          Handling of arrays larger than 32767   Lars L. Christensen, march, 1995
          Converted to IDL V5.0   W. Landsman   September 1997
          Fixed bug where A term not initialized for vector input 
              W. Landsman        February 2000
  """

  radeg = 180.0/numpy.pi
  sec_to_radian = 1.0/radeg/3600.0

  M = numpy.array( [ [+0.9999256795, -0.0111814828, -0.0048590040,  \
        -0.000551,  -0.238560,     +0.435730     ], \
      [ +0.0111814828, +0.9999374849, -0.0000271557,  \
        +0.238509,     -0.002667,      -0.008541     ], \
      [ +0.0048590039, -0.0000271771, +0.9999881946 , \
        -0.435614,      +0.012254,      +0.002117      ], \
      [ -0.00000242389840, +0.00000002710544, +0.00000001177742, \
        +0.99990432,    -0.01118145,    -0.00485852    ], \
      [ -0.00000002710544, -0.00000242392702, +0.00000000006585, \
        +0.01118145,     +0.99991613,    -0.00002716    ], \
      [ -0.00000001177742, +0.00000000006585,-0.00000242404995, \
        +0.00485852,   -0.00002717,    +0.99996684] ] )

  A_dot = 1e-3*numpy.array([1.244, -1.579, -0.660 ])           # in arc seconds per century

  ra_rad = ra/radeg       ;      dec_rad = dec/radeg
  cosra =  cos( ra_rad )  ;       sinra = sin( ra_rad )
  cosdec = cos( dec_rad ) ;      sindec = sin( dec_rad )

  dec_1950 = 0.0
  ra_1950 = 0.0


  # Following statement moved inside loop in Feb 2000.
  A = 1e-6*numpy.array([ -1.62557, -0.31919, -0.13843])        # in radians

  r0 = numpy.array([ cosra*cosdec, sinra*cosdec, sindec ])

  if mu_radec is not None:
    mu_a = mu_radec[0]
    mu_d = mu_radec[1]
    # Velocity vector
    r0_dot = numpy.array([ -mu_a*sinra*cosdec - mu_d*cosra*sindec , \
              mu_a*cosra*cosdec - mu_d*sinra*sindec , \
              mu_d*cosdec ]) + 21.095 * rad_vel * parallax * r0
  else:
    r0_dot = numpy.zeros(3)

  R_0 = numpy.concatenate([ r0, r0_dot ])
  R_1 =  numpy.dot(numpy.transpose(M), R_0)

  # Include the effects of the E-terms of aberration to form r and r_dot.

  r1 = R_1[0:3]
  r1_dot = R_1[3:6]

  if mu_radec is None:
    r1 = r1 + sec_to_radian * r1_dot * (epoch - 1950.0)/100.
    A = A + sec_to_radian * A_dot * (epoch - 1950.0)/100.

  x1 = R_1[0]   ;   y1 = R_1[1]    ;  z1 = R_1[2]
  rmag = sqrt( x1**2 + y1**2 + z1**2 )

  s1 = r1/rmag    ; s1_dot = r1_dot/rmag

  s = s1
  for j in smo.range(2):
    r = s1 + A - ((s * A).sum())*s
    s = r/rmag
 
  x = r[0]          ; y = r[1]     ;  z = r[2]
  r2 = x**2 + y**2 + z**2
  rmag = sqrt( r2 )
 
  if mu_radec is not None:
    r_dot = s1_dot + A_dot - ( ( s * A_dot).sum() )*s
    x_dot = r_dot[0]  ; y_dot= r_dot[1]  ;  z_dot = r_dot[2]
    mu_radec[0] = ( x*y_dot - y*x_dot) / ( x**2 + y**2)
    mu_radec[1] = ( z_dot* (x**2 + y**2) - z*(x*x_dot + y*y_dot) ) / ( r2*sqrt( x**2 + y**2) )

  dec_1950 = arcsin( z / rmag)
  ra_1950 = numpy.arctan2( y, x)

  if parallax > 0.0:
    rad_vel = ( x*x_dot + y*y_dot + z*z_dot )/ (21.095*parallax*rmag)
    parallax = parallax / rmag


  if ra_1950 < 0.0:
    ra_1950 += 2.0*numpy.pi

  ra_1950 = ra_1950*radeg ; dec_1950 = dec_1950*radeg

  return [ra_1950, dec_1950, mu_radec, parallax, rad_vel]





def premat( equinox1, equinox2, FK4 = False):
  """
    Return the precession matrix needed to go from EQUINOX1 to EQUINOX2.
  
    Parameters
    ----------
    equinox1, equinox2 : float
        EQUINOX1 and EQUINOX2
    FK4 : boolean
        Set this to True to obtain output in FK4 system

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:

    NAME:
          PREMAT
    PURPOSE:
          Return the precession matrix needed to go from EQUINOX1 to EQUINOX2.  
    EXPLANTION:
          This matrix is used by the procedures PRECESS and BARYVEL to precess 
          astronomical coordinates

    CALLING SEQUENCE:
          matrix = PREMAT( equinox1, equinox2, [ /FK4 ] )

    INPUTS:
          EQUINOX1 - Original equinox of coordinates, numeric scalar.  
          EQUINOX2 - Equinox of precessed coordinates.

    OUTPUT:
        matrix - double precision 3 x 3 precession matrix, used to precess
                  equatorial rectangular coordinates

    OPTIONAL INPUT KEYWORDS:
          /FK4   - If this keyword is set, the FK4 (B1950.0) system precession
                  angles are used to compute the precession matrix.   The 
                  default is to use FK5 (J2000.0) precession angles

    EXAMPLES:
          Return the precession matrix from 1950.0 to 1975.0 in the FK4 system

          IDL> matrix = PREMAT( 1950.0, 1975.0, /FK4)

    PROCEDURE:
          FK4 constants from "Computational Spherical Astronomy" by Taff (1983), 
          p. 24. (FK4). FK5 constants from "Astronomical Almanac Explanatory
          Supplement 1992, page 104 Table 3.211.1.

    REVISION HISTORY
          Written, Wayne Landsman, HSTX Corporation, June 1994
          Converted to IDL V5.0   W. Landsman   September 1997
  """

  deg_to_rad = numpy.pi/180.0
  sec_to_rad = deg_to_rad/3600.0

  t = 0.001*(equinox2 - equinox1)

  if not FK4:
    st = 0.001*(equinox1 - 2000.0)
    # Compute 3 rotation angles
    A = sec_to_rad * t * (23062.181 + st*(139.656 +0.0139*st) \
        + t*(30.188 - 0.344*st+17.998*t))

    B = sec_to_rad * t * t * (79.280 + 0.410*st + 0.205*t) + A

    C = sec_to_rad * t * (20043.109 - st*(85.33 + 0.217*st) \
        + t*(-42.665 - 0.217*st -41.833*t))

  else:
    st = 0.001*(equinox1 - 1900.0)
    # Compute 3 rotation angles

    A = sec_to_rad * t * (23042.53 + st*(139.75 +0.06*st) \
        + t * (30.23 - 0.27*st+18.0*t))

    B = sec_to_rad * t * t * (79.27 + 0.66*st + 0.32*t) + A

    C = sec_to_rad * t * (20046.85 - st*(85.33 + 0.37*st) \
        + t*(-42.67 - 0.37*st -41.8*t))

  sina = sin(A) ;  sinb = sin(B)  ; sinc = sin(C)
  cosa = cos(A) ;  cosb = cos(B)  ; cosc = cos(C)

  r = numpy.zeros((3,3))
  r[::,0] = numpy.array([ cosa*cosb*cosc-sina*sinb, sina*cosb+cosa*sinb*cosc,  cosa*sinc])
  r[::,1] = numpy.array([-cosa*sinb-sina*cosb*cosc, cosa*cosb-sina*sinb*cosc, -sina*sinc])
  r[::,2] = numpy.array([-cosb*sinc, -sinb*sinc, cosc])
  return r







def precess(ra, dec, equinox1, equinox2, FK4 = False, radian=False):
  """
    Precess coordinates from EQUINOX1 to EQUINOX2.
  
    Parameters
    ----------
    ra, dec, equinox1, equinox2 : float
        Position and equinox
    FK4 : boolean
        Set to True to obtain output in FK4 system.
    radian : boolean
        If True, `ra` and `dec` must be given in radian (degrees otherwise). 

    Returns
    -------
    Position : list of ra and dec
        A list with [ra, dec] precessed from equinox 1 to equinox 2.

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:

    NAME:
        PRECESS
    PURPOSE:
        Precess coordinates from EQUINOX1 to EQUINOX2.  
    EXPLANATION:
        For interactive display, one can use the procedure ASTRO which calls 
        PRECESS or use the /PRINT keyword.   The default (RA,DEC) system is 
        FK5 based on epoch J2000.0 but FK4 based on B1950.0 is available via 
        the /FK4 keyword.

        Use BPRECESS and JPRECESS to convert between FK4 and FK5 systems
    CALLING SEQUENCE:
        PRECESS, ra, dec, [ equinox1, equinox2, /PRINT, /FK4, /RADIAN ]

    INPUT - OUTPUT:
        RA - Input right ascension (scalar or vector) in DEGREES, unless the 
                /RADIAN keyword is set
        DEC - Input declination in DEGREES (scalar or vector), unless the 
                /RADIAN keyword is set

        The input RA and DEC are modified by PRECESS to give the 
        values after precession.

    OPTIONAL INPUTS:
        EQUINOX1 - Original equinox of coordinates, numeric scalar.  If 
                  omitted, then PRECESS will query for EQUINOX1 and EQUINOX2.
        EQUINOX2 - Equinox of precessed coordinates.

    OPTIONAL INPUT KEYWORDS:
        /PRINT - If this keyword is set and non-zero, then the precessed
                  coordinates are displayed at the terminal.    Cannot be used
                  with the /RADIAN keyword
        /FK4   - If this keyword is set and non-zero, the FK4 (B1950.0) system
                  will be used otherwise FK5 (J2000.0) will be used instead.
        /RADIAN - If this keyword is set and non-zero, then the input and 
                  output RA and DEC vectors are in radians rather than degrees

    RESTRICTIONS:
          Accuracy of precession decreases for declination values near 90 
          degrees.  PRECESS should not be used more than 2.5 centuries from
          2000 on the FK5 system (1950.0 on the FK4 system).

    EXAMPLES:
          (1) The Pole Star has J2000.0 coordinates (2h, 31m, 46.3s, 
                  89d 15' 50.6"); compute its coordinates at J1985.0

          IDL> precess, ten(2,31,46.3)*15, ten(89,15,50.6), 2000, 1985, /PRINT

                  ====> 2h 16m 22.73s, 89d 11' 47.3"

          (2) Precess the B1950 coordinates of Eps Ind (RA = 21h 59m,33.053s,
          DEC = (-56d, 59', 33.053") to equinox B1975.

          IDL> ra = ten(21, 59, 33.053)*15
          IDL> dec = ten(-56, 59, 33.053)
          IDL> precess, ra, dec ,1950, 1975, /fk4

    PROCEDURE:
          Algorithm from Computational Spherical Astronomy by Taff (1983), 
          p. 24. (FK4). FK5 constants from "Astronomical Almanac Explanatory
          Supplement 1992, page 104 Table 3.211.1.

    PROCEDURE CALLED:
          Function PREMAT - computes precession matrix 

    REVISION HISTORY
          Written, Wayne Landsman, STI Corporation  August 1986
          Correct negative output RA values   February 1989
          Added /PRINT keyword      W. Landsman   November, 1991
          Provided FK5 (J2000.0)  I. Freedman   January 1994
          Precession Matrix computation now in PREMAT   W. Landsman June 1994
          Added /RADIAN keyword                         W. Landsman June 1997
          Converted to IDL V5.0   W. Landsman   September 1997
          Correct negative output RA values when /RADIAN used    March 1999 
          Work for arrays, not just vectors  W. Landsman    September 2003 
  """
  deg_to_rad = numpy.pi/180.0

  if not radian:
    # ra, dec are given in degrees
    ra_rad = ra*deg_to_rad     # Convert to double precision if not already
    dec_rad = dec*deg_to_rad
  else:
    ra_rad = ra ; dec_rad = dec

  a = cos( dec_rad )


  x = [a*cos(ra_rad), a*sin(ra_rad), sin(dec_rad)] # input direction

  sec_to_rad = deg_to_rad/3600.0

  # Use PREMAT function to get precession matrix from Equinox1 to Equinox2

  r = premat(equinox1, equinox2, FK4 = FK4)

  x2 = numpy.dot(r,x)      # rotate to get output direction cosines

  ra_rad = numpy.arctan2(x2[1],x2[0])
  dec_rad = arcsin(x2[2])

  if not radian:
    ra = ra_rad/deg_to_rad
    ra = ra + int(ra < 0.0)*360.0            # RA between 0 and 360 degrees
    dec = dec_rad/deg_to_rad
  else:
    ra = ra_rad ; dec = dec_rad
    ra = ra + int(ra < 0.0)*2.0*numpy.pi

  return [ra, dec]








def precess_xyz(x,y,z,equinox1,equinox2):
  """
    Precess equatorial geocentric rectangular coordinates. 
  
    Parameters
    ----------
    x, y, z, equinox1, equinox2 : float
    
    Returns
    -------
    Precessed coordinates : list
        A list containing the updated `x`, `y`, and `z` values. 
    
    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.
    
    :IDL - Documentation:

    NAME:
          PRECESS_XYZ

    PURPOSE:
          Precess equatorial geocentric rectangular coordinates. 

    CALLING SEQUENCE:
          precess_xyz, x, y, z, equinox1, equinox2

    INPUT/OUTPUT:
          x,y,z: scalars or vectors giving heliocentric rectangular coordinates
                THESE ARE CHANGED UPON RETURNING.
    INPUT:
          EQUINOX1: equinox of input coordinates, numeric scalar
          EQUINOX2: equinox of output coordinates, numeric scalar

    OUTPUT:
          x,y,z are changed upon return

    NOTES:
      The equatorial geocentric rectangular coords are converted
        to RA and Dec, precessed in the normal way, then changed
        back to x, y and z using unit vectors.

    EXAMPLE:
          Precess 1950 equinox coords x, y and z to 2000.
          IDL> precess_xyz,x,y,z, 1950, 2000

    HISTORY:
          Written by P. Plait/ACC March 24 1999 
            (unit vectors provided by D. Lindler)
          Use /Radian call to PRECESS     W. Landsman     November 2000
          Use two parameter call to ATAN   W. Landsman    June 2001
  """

  # take input coords and convert to ra and dec (in radians)
  ra = numpy.arctan2(y,x)
  delp = sqrt(x*x + y*y + z*z)  # magnitude of distance to Sun
  dec = arcsin(z/delp)

  # precess the ra and dec
  ra, dec = precess(ra, dec, equinox1, equinox2, radian=True)

  # convert back to x, y, z
  xunit = cos(ra)*cos(dec)
  yunit = sin(ra)*cos(dec)
  zunit = sin(dec)

  x = xunit * delp
  y = yunit * delp
  z = zunit * delp
  
  return [x,y,z]








def xyz(date, velocity=False, equinox=1950.0):
  """
    Calculate geocentric X,Y, and Z  and velocity coordinates of the Sun.
  
    Parameters
    ----------
    date : float
        Julian date
    equinox : float
        Equinox of output. If None, Equinox will be 1950.
    velocity : boolean
        If False, the velocity of the Sun will not be calculated
    
    Returns
    -------
    Sun position and velocity : list
        A list of the from [X, Y, Z, XVEL, YVEL, ZVEL]. Last three values are None
        if `velocity` flag is set to False.

    Notes
    -----

    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:

    NAME:
          XYZ
    PURPOSE:
          Calculate geocentric X,Y, and Z  and velocity coordinates of the Sun
    EXPLANATION:
          Calculates geocentric X,Y, and Z vectors and velocity coordinates 
          (dx, dy and dz) of the Sun.   (The positive X axis is directed towards 
          the equinox, the y-axis, towards the point on the equator at right 
          ascension 6h, and the z axis toward the north pole of the equator).
          Typical position accuracy is <1e-4 AU (15000 km).

    CALLING SEQUENCE:
          XYZ, date, x, y, z, [ xvel, yvel, zvel, EQUINOX = ]

    INPUT:
          date: reduced julian date (=JD - 2400000), scalar or vector

    OUTPUT:
          x,y,z: scalars or vectors giving heliocentric rectangular coordinates
                    (in A.U) for each date supplied.    Note that sqrt(x^2 + y^2
                    + z^2) gives the Earth-Sun distance for the given date.
          xvel, yvel, zvel: velocity vectors corresponding to X, Y and Z.

    OPTIONAL KEYWORD INPUT:
          EQUINOX: equinox of output. Default is 1950.

    EXAMPLE:
          What were the rectangular coordinates and velocities of the Sun on 
          Jan 22, 1999 0h UT (= JD 2451200.5) in J2000 coords? NOTE:
          Astronomical Almanac (AA) is in TDT, so add 64 seconds to 
          UT to convert.

          IDL> xyz,51200.5+64.d/86400.d,x,y,z,xv,yv,zv,equinox = 2000

          Compare to Astronomical Almanac (1999 page C20)
                      X  (AU)        Y  (AU)     Z (AU)
          XYZ:      0.51456871   -0.76963263  -0.33376880
          AA:       0.51453130   -0.7697110   -0.3337152
          abs(err): 0.00003739    0.00007839   0.00005360
          abs(err)
              (km):   5609          11759         8040 

          NOTE: Velocities in AA are for Earth/Moon barycenter
                (a very minor offset) see AA 1999 page E3
                    X VEL (AU/DAY) YVEL (AU/DAY)   Z VEL (AU/DAY)
          XYZ:      -0.014947268   -0.0083148382    -0.0036068577
          AA:       -0.01494574    -0.00831185      -0.00360365
          abs(err):  0.000001583    0.0000029886     0.0000032077
          abs(err)
          (km/sec): 0.00265        0.00519          0.00557

    PROCEDURE CALLS:
          PRECESS_XYZ
    REVISION HISTORY
          Original algorithm from Almanac for Computers, Doggett et al. USNO 1978
          Adapted from the book Astronomical Photometry by A. Henden
          Written  W. Landsman   STX       June 1989
          Correct error in X coefficient   W. Landsman HSTX  January 1995
          Added velocities, more terms to positions and EQUINOX keyword,
            some minor adjustments to calculations 
            P. Plait/ACC March 24, 1999
  """

  picon = numpy.pi/180.0
  t = (date - 15020.0)/36525.0         # Relative Julian century from 1900

  # NOTE: longitude arguments below are given in *equinox* of date.
  #   Precess these to equinox 1950 to give everything an even footing.
  #   Compute argument of precession from equinox of date back to 1950
  pp = (1.396041 + 0.000308*(t + 0.5))*(t-0.499998)

  # Compute mean solar longitude, precessed back to 1950
  el = 279.696678 + 36000.76892*t + 0.000303*t*t - pp

  # Compute Mean longitude of the Moon
  c = 270.434164 + 480960.*t + 307.883142*t - 0.001133*t*t - pp

  # Compute longitude of Moon's ascending node
  n = 259.183275 - 1800.*t - 134.142008*t + 0.002078*t*t - pp

  # Compute mean solar anomaly
  g = 358.475833 + 35999.04975*t - 0.00015*t*t

  # Compute the mean jupiter anomaly
  j = 225.444651 + 2880.0*t + 154.906654*t*t

  # Compute mean anomaly of Venus
  v = 212.603219 + 58320.*t + 197.803875*t + 0.001286*t*t

  # Compute mean anomaly of Mars
  m = 319.529425 + 19080.*t + 59.8585*t + 0.000181*t*t

  # Convert degrees to radians for trig functions
  el = el*picon
  g  = g*picon
  j =  j*picon
  c  = c*picon
  v  = v*picon
  n  = n*picon
  m  = m*picon

  # Calculate X,Y,Z using trigonometric series
  X =   0.999860*cos(el)                          \
      - 0.025127*cos(g - el)                      \
      + 0.008374*cos(g + el)                      \
      + 0.000105*cos(g + g + el)                  \
      + 0.000063*t*cos(g - el)                    \
      + 0.000035*cos(g + g - el)                  \
      - 0.000026*sin(g - el - j)                  \
      - 0.000021*t*cos(g + el)                    \
      + 0.000018*sin(2.*g + el - 2.*v)            \
      + 0.000017*cos(c)                           \
      - 0.000014*cos(c - 2.*el)                   \
      + 0.000012*cos(4.*g + el - 8.*m + 3.*j)     \
      - 0.000012*cos(4.*g - el - 8.*m + 3.*j)     \
      - 0.000012*cos(g + el - v)                  \
      + 0.000011*cos(2.*g + el - 2.*v)            \
      + 0.000011*cos(2.*g - el - 2.*j)
  

  Y =   0.917308*sin(el)                             \
      + 0.023053*sin(g - el)                         \
      + 0.007683*sin(g + el)                         \
      + 0.000097*sin(g + g + el)                     \
      - 0.000057*t*sin(g - el)                       \
      - 0.000032*sin(g + g - el)                     \
      - 0.000024*cos(g - el - j)                     \
      - 0.000019*t*sin(g + el)                       \
      - 0.000017*cos(2.00*g + el - 2.00*v)           \
      + 0.000016*sin(c)                              \
      + 0.000013*sin(c - 2.00*el )                   \
      + 0.000011*sin(4.00*g + el - 8.00*m + 3.00*j)  \
      + 0.000011*sin(4.00*g - el - 8.00*m + 3.00*j)  \
      - 0.000011*sin(g + el - v)                     \
      + 0.000010*sin(2.00*g + el - 2.00*v )          \
      - 0.000010*sin(2.00*g - el - 2.00*j )


  Z =   0.397825*sin(el)        \
      + 0.009998*sin(g-el)      \
      + 0.003332*sin(g+el)      \
      + 0.000042*sin(g+g+el)    \
      - 0.000025*t*sin(g-el)    \
      - 0.000014*sin(g+g-el)    \
      - 0.000010*cos(g-el-j)

  # Precess_to new equator?
  if equinox is not None:
    X, Y, Z = precess_xyz(X, Y, Z, 1950, equinox)

  if not velocity:
    return [X, Y, Z, None, None, None]
  
  XVEL = -0.017200 * sin(el)          \
        -0.000288 * sin(g + el)       \
        -0.000005 * sin(2.00*g + el)  \
        -0.000004 * sin(c)            \
        +0.000003 * sin(c - 2.00*el)  \
        +0.000001 *t * sin(g+el)      \
        -0.000001 * sin(2.00*g-el)

  YVEL =  0.015780 * cos(el)           \
        +0.000264 * cos(g + el)        \
        +0.000005 * cos(2.00*g + el)   \
        +0.000004 * cos(c)             \
        +0.000003 * cos(c - 2.00*el)   \
        -0.000001 * t * cos(g + el)

  ZVEL = 0.006843 * cos(el)             \
        +0.000115 * cos(g  + el)        \
        +0.000002 * cos(2.00*g + el)    \
        +0.000002 * cos(c)              \
        +0.000001 * cos(c - 2.00*el)

  # Precess to new equator?

  if equinox != 1950.0:
    XVEL, YVEL, ZVEL = precess_xyz(XVEL, YVEL, ZVEL, 1950, equinox)
  
  return [X, Y, Z, XVEL, YVEL, ZVEL]












def helio_jd(date, ra, dec, B1950 = False, TIME_DIFF = False):
  """
    Convert geocentric (reduced) Julian date to heliocentric Julian date
  
    Parameters
    ----------
    date : float
        (Reduced) Julian date (2.4e6 subtracted)
    ra, dec : float
        Right ascension and declination in degrees
    B1950 : boolean
        If True, input coordinates are assumed to be given in equinox
        1950 coordinates.
    TIME_DIFF : boolean
        If True, this function returns the time difference
        (heliocentric JD - geocentric JD ) in seconds
    
    Returns
    -------
    HJD : float
        The heliocentric Julian date.
      
    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.
    
    :IDL - Documentation:
    
    NAME:
         HELIO_JD
    PURPOSE:
         Convert geocentric (reduced) Julian date to heliocentric Julian date
    EXPLANATION:
         This procedure correct for the extra light travel time between the Earth 
         and the Sun.
    
          An online calculator for this quantity is available at 
          http://www.physics.sfasu.edu/astro/javascript/hjd.html
    CALLING SEQUENCE:
          jdhelio = HELIO_JD( date, ra, dec, /B1950, /TIME_DIFF)
    
    INPUTS
          date - reduced Julian date (= JD - 2400000), scalar or vector, MUST
                  be double precision
          ra,dec - scalars giving right ascension and declination in DEGREES
                  Equinox is J2000 unless the /B1950 keyword is set
    
    OUTPUTS:
          jdhelio - heliocentric reduced Julian date.  If /TIME_DIFF is set, then
                    HELIO_JD() instead returns the time difference in seconds
                    between the geocentric and heliocentric Julian date.
    
    OPTIONAL INPUT KEYWORDS 
          /B1950 - if set, then input coordinates are assumed to be in equinox 
                   B1950 coordinates.
          /TIME_DIFF - if set, then HELIO_JD() returns the time difference
                   (heliocentric JD - geocentric JD ) in seconds 
    
    EXAMPLE:
          What is the heliocentric Julian date of an observation of V402 Cygni
          (J2000: RA = 20 9 7.8, Dec = 37 09 07) taken June 15, 1973 at 11:40 UT?
    
          IDL> juldate, [1973,6,15,11,40], jd      ;Get geocentric Julian date
          IDL> hjd = helio_jd( jd, ten(20,9,7.8)*15., ten(37,9,7) )  
    
          ==> hjd = 41848.9881
    
    Wayne Warren (Raytheon ITSS) has compared the results of HELIO_JD with the
    FORTRAN subroutines in the STARLINK SLALIB library (see 
    http://star-www.rl.ac.uk/).    
                                                     Time Diff (sec)
         Date               RA(2000)   Dec(2000)  STARLINK      IDL
    
    1999-10-29T00:00:00.0  21 08 25.  -67 22 00.  -59.0        -59.0
    1999-10-29T00:00:00.0  02 56 33.4 +00 26 55.  474.1        474.1
    1940-12-11T06:55:00.0  07 34 41.9 -00 30 42.  366.3        370.2
    1992-02-29T03:15:56.2  12 56 27.4 +42 10 17.  350.8        350.9
    2000-03-01T10:26:31.8  14 28 36.7 -20 42 11.  243.7        243.7
    2100-02-26T09:18:24.2  08 26 51.7 +85 47 28.  104.0        108.8
    PROCEDURES CALLED:
          bprecess, xyz, zparcheck
    
    REVISION HISTORY:
          Algorithm from the book Astronomical Photometry by Henden, p. 114
          Written,   W. Landsman       STX     June, 1989 
          Make J2000 default equinox, add B1950, /TIME_DIFF keywords, compute
          variation of the obliquity      W. Landsman   November 1999
  """

  # Because XYZ uses default B1950 coordinates, we'll convert everything to B1950

  if date > 2.4e6:
    PE.warn(PE.PyAValError("The given Julian Date ( " + str(date) + ") is exceedingly large far a reduced JD.",
                           solution="Did you forget to subtract 2.4e6?",
                           where="helio_jd"))

  if not B1950:
    bpresult = bprecess(ra,dec)
    ra1 = bpresult[0]
    dec1 = bpresult[1]
  else:
    ra1 = ra
    dec1 = dec

  radeg = 180.0/numpy.pi
  # I think, this is not needed in Python, even at this stage...
  # zparcheck,'HELIO_JD',date,1,[3,4,5],[0,1],'Reduced Julian Date'

  delta_t = (date - 33282.42345905)/36525.0
  epsilon_sec = 44.836 - 46.8495*delta_t - 0.00429*delta_t**2 + 0.00181*delta_t**3
  epsilon = (23.433333 + epsilon_sec/3600.0)/radeg
  ra1 = ra1/radeg
  dec1 = dec1/radeg

  x, y, z, tmp, tmp, tmp = xyz(date)

  # Find extra distance light must travel in AU, multiply by 1.49598e13 cm/AU,
  # and divide by the speed of light, and multiply by 86400 second/year

  time = -499.00522*( cos(dec1)*cos(ra1)*x + \
                  (tan(epsilon)*sin(dec1) + cos(dec1)*sin(ra1))*y)

  if TIME_DIFF:
    return time
  else:
    return (date + time/86400.0)


def jdcnv(dt):
  """
    Converts Gregorian dates to Julian days  
  
    Parameters
    ----------
    dt : DateTime object
        The date.
    
    Returns
    -------
    Julian day : float

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.
    
    :IDL - Documentation:

     NAME:
           JDCNV
     PURPOSE:
           Converts Gregorian dates to Julian days   

     EXPLANATION:
           For IDL versions V5.1 or greater, this procedure is superceded by
           JULDAY() function in the standard IDL distribution.   Note, however,
           that prior to V5.1 there wasa bug in JULDAY() that gave answers off
           by 0.5 days. 

     CALLING SEQUENCE:
           JDCNV, YR, MN, DAY, HR, JULIAN

     INPUTS:
           YR = Year, integer scalar or vector
           MN = Month  integer (1-12) scalar or vector
           DAY = Day   integer 1-31) scalar or vector 
           HR  = Hours and fractions of hours of universal time (U.T.), scalar
                  or vector

     OUTPUTS:
           JULIAN = Julian date (double precision) 

     EXAMPLE:
           To find the Julian Date at 1978 January 1, 0h (U.T.)

           IDL> JDCNV, 1978, 1, 1, 0., JULIAN

           will give JULIAN = 2443509.5
     NOTES:
           (1) JDCNV will accept vector arguments 
           (2) JULDATE is an alternate procedure to perform the same function

     REVISON HISTORY:
           Converted to IDL from Don Yeomans Comet Ephemeris Generator,
           B. Pfarr, STX, 6/15/88
           Converted to IDL V5.0   W. Landsman   September 1997
           Added checks on valid month, day ranges W. Landsman July 2008
  """

  # L = leap years, -1 for Jan, Feb, else 0
  L = int((dt.month-14.0)/12.0)
  julian = dt.day - 32075 + \
           int(1461*(dt.year+4800+L)/4.0) + \
           int(367*(dt.month - 2-L*12)/12.0) - int(int(3*((dt.year+4900+L)/100.0))/4.0)
  julian += ((dt.hour/24.0) + (dt.minute/(24.0*60.0)) + (dt.second/86400.) + (dt.microsecond/(86400.*1e6)) - 0.5)

  return julian





def get_juldate():
  """
    Return the current Julian Date
  
    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.
  
    :IDL - Documentation:

    NAME:
        GET_JULDATE
    PURPOSE:
        Return the current Julian Date

    EXPLANATION:
        In V5.4, GET_JULDATE became completely obsolete with the introduction
        of the /UTC keyword to SYSTIME().   So GET_JULDATE,jd is equivalent to
        jd = SYSTIME(/JULIAN,/UTC).

    CALLING SEQUENCE:
          GET_JULDATE,jd

    INPUTS:
          None

    OUTPUTS:
          jd = Current Julian Date, double precision scalar

    EXAMPLE:
          Return the current hour, day, month and year as integers

          IDL> GET_JULDATE, JD                  ;Get current Julian date
          IDL> DAYCNV, JD, YR, MON, DAY, HOURS  ;Convert to hour,day month & year

    METHOD:
          A call is made to SYSTIME(/JULIAN,/UTC).

    REVISION HISTORY:
          Written Wayne Landsman                March, 1991
          Converted to IDL V5.0   W. Landsman   September 1997
          Assume since V5.4 Use /UTC keyword to SYSTIME()  W. Landsman April 2006
  """
  now = datetime.datetime.utcnow()
  return jdcnv(now)





def juldate(date):
  """
    Convert from calendar to Reduced Julian Date
  
    This function returns the *reduced* Julian date, which
    is obtained by subtracting 2400000 from the Julian date.
    To convert the output into Modified Julian Date (MJD),
    another 0.5 days have to be subtracted.
  
    Parameters
    ----------
    date : DateTime object
        Calendar date
    
    Returns
    -------
    RJD : float
        The **reduced** Julian date.
  
    Notes
    -----

    .. note:: This function was ported from the IDL Astronomy User's Library.
  
    :IDL - Documentation:

    NAME:
      JULDATE
    PURPOSE:
      Convert from calendar to Reduced Julian Date

    EXPLANATION:
      Julian Day Number is a count of days elapsed since Greenwich mean noon 
      on 1 January 4713 B.C.  The Julian Date is the Julian day number
      followed by the fraction of the day elapsed since the preceding noon. 

      This procedure duplicates the functionality of the JULDAY() function in
      in the standard IDL distribution, but also allows interactive input and
      gives output as Reduced Julian date (=JD - 2400000.)  
      (Also note that prior to V5.1 there was a bug in JULDAY() that gave 
      answers offset by 0.5 days.)

    CALLING SEQUENCE:
      JULDATE, /PROMPT           ;Prompt for calendar Date, print Julian Date
                or
      JULDATE, date, jd      

    INPUT:
      DATE -  3 to 6-element vector containing year,month (1-12),day, and 
                optionally hour, minute, and second all specified as numbers
                (Universal Time).   Year should be supplied with all digits.
                Years B.C should be entered as negative numbers (and note that
                Year 0 did not exist).  If Hour, minute or seconds are not 
                supplied, they will default to 0. 

    OUTPUT:
        JD - Reduced Julian date, double precision scalar.  To convert to
                Julian Date, add 2400000.   JULDATE will print the value of
                JD at the terminal if less than 2 parameters are supplied, or 
                if the /PROMPT keyword is set

    OPTIONAL INPUT KEYWORD:
        /PROMPT - If this keyword is set and non-zero, then JULDATE will prompt
                for the calendar date at the terminal.

    RESTRICTIONS:
        The procedure HELIO_JD can be used after JULDATE, if a heliocentric
        Julian date is required.

    EXAMPLE:
        A date of 25-DEC-2006 06:25 UT may be expressed as either

        IDL> juldate, [2006, 12, 25, 6, 25], jd       
        IDL> juldate, [2006, 12, 25.2673611d], jd 

        In either case, one should obtain a Reduced Julian date of 
        JD = 54094.7673611

    PROCEDURE USED:
        GETOPT()
    REVISION HISTORY
        Adapted from IUE RDAF (S. Parsons)                      8-31-87
        Algorithm from Sky and Telescope April 1981   
        Added /PROMPT keyword, W. Landsman    September 1992
        Converted to IDL V5.0   W. Landsman   September 1997
        Make negative years correspond to B.C. (no year 0), work for year 1582
        Disallow 2 digit years.    W. Landsman    March 2000
  """
  jd = jdcnv(date)
  jd -= 2400000.0
  return jd
