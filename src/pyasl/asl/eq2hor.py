# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
from .sunpos import sunpos
from .observatory import observatory as observatory_function
from .idlMod import idlMod
from .idlCirrange import cirrange
from PyAstronomy.pyasl import _ic
from . import astroTimeLegacy as atleg
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo



def nutate(jd, radian=False, plot=False):
  """
    Computes the Earth's nutation in longitude and obliquity for a given (array) of Julian date.
    
    .. warning:: The output of the IDL routine is in units of arcseconds, whereas the default
                 if this routine returns degrees.
    
    Parameters
    ----------
    jd : float
         The Julian date
    radian : boolean, optional
         Results are returned in radian instead of in degrees.
         The default is False.
    plot : boolean, optional
         Results are plotted. The default is False.
         
    Returns
    -------
    Longitude : float
        The nutation in longitude (in deg by default).
    Obliquity : float
        The nutation in latitude (in deg by default).
        

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:  
    
    NAME:
          NUTATE
    PURPOSE:
          Return the nutation in longitude and obliquity for a given Julian date
    
    CALLING SEQUENCE:
          NUTATE, jd, Nut_long, Nut_obliq
    
    INPUT:
          jd - Julian ephemeris date, scalar or vector, double precision  
    OUTPUT:
          Nut_long - the nutation in longitude, same # of elements as jd
          Nut_obliq - nutation in latitude, same # of elements as jd
    
    EXAMPLE:
          (1) Find the nutation in longitude and obliquity 1987 on Apr 10 at Oh.
                 This is example 22.a from Meeus
           IDL> jdcnv,1987,4,10,0,jul
           IDL> nutate, jul, nut_long, nut_obliq
                ==> nut_long = -3.788    nut_obliq = 9.443
               
          (2) Plot the large-scale variation of the nutation in longitude 
                  during the 20th century
    
          IDL> yr = 1900 + indgen(100)     ;Compute once a year        
          IDL> jdcnv,yr,1,1,0,jul          ;Find Julian date of first day of year
          IDL> nutate,jul, nut_long        ;Nutation in longitude
          IDL> plot, yr, nut_long
    
          This plot will reveal the dominant (18.6 year) period, but a finer
          grid is needed to display the shorter periods in the nutation.
    METHOD:
          Uses the formula in Chapter 22 of ``Astronomical Algorithms'' by Jean 
          Meeus (1998, 2nd ed.) which is based on the 1980 IAU Theory of Nutation
          and includes all terms larger than 0.0003".
    
    PROCEDURES CALLED:
          POLY()                       (from IDL User's Library)
          CIRRANGE, ISARRAY()          (from IDL Astronomy Library)
    
    REVISION HISTORY:
          Written, W.Landsman (Goddard/HSTX)      June 1996       
          Converted to IDL V5.0   W. Landsman   September 1997
          Corrected minor typos in values of d_lng W. Landsman  December 2000
          Updated typo in cdelt term              December 2000
          Avoid overflow for more than 32767 input dates W. Landsman January 2005    
  """

  #form time in Julian centuries from 1900.0
  jdcen = (np.array(jd, ndmin=1) - 2451545.0)/36525.0
  
  #Mean elongation of the Moon
  coef_moon = [1.0/189474.0, -0.0019142, 445267.111480, 297.85036]
  d = np.polyval(coef_moon, jdcen)*np.pi/180.
  d = cirrange(d, radians=True)

  #Sun's mean anomaly
  coef_sun = [-1.0/3e5, -0.0001603, 35999.050340, 357.52772]
  sun = np.polyval(coef_sun, jdcen)*np.pi/180.
  sun = cirrange(sun, radians=True)

  # Moon's mean anomaly
  coef_mano = [1.0/5.625e4, 0.0086972, 477198.867398, 134.96298]
  mano = np.polyval(coef_mano, jdcen)*np.pi/180.
  mano = cirrange(mano, radians=True)

  # Moon's argument of latitude
  coef_mlat = [-1.0/3.27270e5, -0.0036825, 483202.017538, 93.27191]
  mlat = np.polyval(coef_mlat, jdcen)*np.pi/180.
  mlat = cirrange(mlat, radians=True)
  
  # Longitude of the ascending node of the Moon's mean orbit on the ecliptic,
  #  measured from the mean equinox of the date
  coef_moe = [1.0/4.5e5, 0.0020708, -1934.136261, 125.04452]
  omega = np.polyval(coef_moe, jdcen)*np.pi/180.
  omega = cirrange(omega, radians=True)

  d_lng = np.array([0.,-2.,0.,0.,0.,0.,-2.,0.,0.,-2,-2,-2,0,2,0,2,0,0,-2,0,2,0,0,-2,0,-2,0,0,2, \
           -2,0,-2,0,0,2,2,0,-2,0,2,2,-2,-2,2,2,0,-2,-2,0,-2,-2,0,-1,-2,1,0,0,-1,0,0, \
           2,0,2], float)

  m_lng = np.array([0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,1,0,-1,0,0,0,1,1,-1,0, \
           0,0,0,0,0,-1,-1,0,0,0,1,0,0,1,0,0,0,-1,1,-1,-1,0,-1], float)

  mp_lng = np.array([0,0,0,0,0,1,0,0,1,0,1,0,-1,0,1,-1,-1,1,2,-2,0,2,2,1,0,0,-1,0,-1, \
            0,0,1,0,2,-1,1,0,1,0,0,1,2,1,-2,0,1,0,0,2,2,0,1,1,0,0,1,-2,1,1,1,-1,3,0], float)

  f_lng = np.array([0,2,2,0,0,0,2,2,2,2,0,2,2,0,0,2,0,2,0,2,2,2,0,2,2,2,2,0,0,2,0,0, \
           0,-2,2,2,2,0,2,2,0,2,2,0,0,0,2,0,2,0,2,-2,0,0,0,2,2,0,0,2,2,2,2], float)

  om_lng = np.array([1,2,2,2,0,0,2,1,2,2,0,1,2,0,1,2,1,1,0,1,2,2,0,2,0,0,1,0,1,2,1, \
            1,1,0,1,2,2,0,2,1,0,2,1,1,1,0,1,1,1,1,1,0,0,0,0,0,2,0,0,2,2,2,2], float)

  sin_lng = np.array([-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217, \
             -158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38, -31, 29, 29, 26, -22, \
             21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7, -7, -7, -7, \
             6,6,6,-6,-6,5,-5,-5,-5,4,4,4,-4,-4,-4,3,-3,-3,-3,-3,-3,-3,-3 ], float)
 
  sdelt = np.array([-174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0., -0.5, 0., 0.1, \
           0.,0.,0.1, 0.,-0.1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., -0.1, 0., 0.1, \
           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.] , float)

  cos_lng = np.array([ 92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95,0,-70,-53,0, \
              -33, 26, 32, 27, 0, -24, 16,13,0,-12,0,0,-10,0,-8,7,9,7,6,0,5,3,-3,0,3,3, \
              0,-3,-3,3,3,0,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ], float)

  cdelt = np.array([8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0, -0.1, 0.3, \
           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., \
           0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.], float)

  # Sum the periodic terms 
  n = jdcen.size
  nut_lon = np.zeros(n)
  nut_obliq = np.zeros(n)
  arg = np.outer(d_lng,d) + np.outer(m_lng,sun) + np.outer(mp_lng,mano) + \
        np.outer(f_lng,mlat) + np.outer(om_lng,omega)
  arg = np.transpose(arg)
  sarg = np.sin(arg)
  carg = np.cos(arg)
  for i in smo.range(n):
    nut_lon[i] = 0.0001*np.sum( (sdelt*jdcen[i] + sin_lng)*sarg[i] )
    nut_obliq[i] = 0.0001*np.sum( (cdelt*jdcen[i] + cos_lng)*carg[i] )
  
  # Until here result are in arcseconds!
  # Convert to degrees
  nut_lon /= 3600.
  nut_obliq /= 3600.
  if radian:
    nut_lon *= (np.pi/180.)
    nut_obliq *= (np.pi/180.)
  
  if plot:
    if not _ic.check["matplotlib"]:
      raise(PE.PyARequiredImport("Could not import matplotlib.", \
            where="nutate", \
            solution=["Install matplotlib", "Switch `plot` flag to False."]))
    import matplotlib.pylab as plt
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    plt.plot(jd, nut_lon)
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel('Julian date')
    if radian==True:
      ax1.set_ylabel('Nutation longitude [red]')
      ax2.set_ylabel('Nutation obliquity [rad]')
    else:
      ax1.set_ylabel('Nutation longitude [deg]')
      ax2.set_ylabel('Nutation obliquity [deg]')
    plt.plot(jd, nut_obliq)
    plt.show()
  
  return nut_lon, nut_obliq


def co_nutate(jd, ra, dec, radian=False, plot=False, full_output=False):
  """
    Compute the changes in RA and DEC due to the Earth's nutation.
    
    Parameters
    ----------
    jd : float or array
         The Julian date. If given as array,
         its size must match that of `ra` and `dec`.
    ra : float or array
         The right ascension in degrees.
         If array, it must be same size as `dec`.
    dec : float or array
         The declination in degrees.
         If array, it must be same size as ra.
    radian : boolean, optional
         Results are returned in radian instead of in degrees.
         The default is False.
    plot : boolean, optional
         If True, the results are plotted.
         For single value `jd`, the change in `ra` and `dec` is plotted
         versus `ra` and `dec`.
         For an array of JDs, ra and dec is plotted versus JD.
         The default is False
    full_output : boolean, optional 
         If True, the result will also contain the obliquity of the ecliptic,
         the nutation in the longitude and the nutation in the
         obliquity of the ecliptic. The default is False.

    Returns
    -------
    dRa : float or array
        The change in right ascension [by default in deg].
    dDec : float or array
        The change in declination [by default in deg].
    True obliquity : float, optional
        The true obliquity of the ecliptic [by default in deg].
        Only if `full_output` is True.
    dLong : float or array, optional
        The nutation in longitude [by default in deg].
        Only if `full_output` is True.
    dObliquity : float or array, optional
        The nutation in the obliquity of the ecliptic [by default in deg].
        Only if `full_output` is True.
        
    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:  

     NAME:
        CO_NUTATE
     PURPOSE:
        Calculate changes in RA and Dec due to nutation of the Earth's rotation
    EXPLANATION:
        Calculates necessary changes to ra and dec due to
        the nutation of the Earth's rotation axis, as described in Meeus, Chap 23.
        Uses formulae from Astronomical Almanac, 1984, and does the calculations
        in equatorial rectangular coordinates to avoid singularities at the
        celestial poles.
    
    CALLING SEQUENCE:
        CO_NUTATE, jd, ra, dec, d_ra, d_dec, [EPS=, D_PSI =, D_EPS = ]
    INPUTS
       JD: Julian Date [scalar or vector]
       RA, DEC : Arrays (or scalars) of the ra and dec's of interest
    
      Note: if jd is a vector, ra and dec MUST be vectors of the same length.
    
    OUTPUTS:
       d_ra, d_dec: the corrections to ra and dec due to nutation (must then
                                   be added to ra and dec to get corrected values).
    OPTIONAL OUTPUT KEYWORDS:
       EPS: set this to a named variable that will contain the obliquity of the 
                ecliptic.
       D_PSI: set this to a named variable that will contain the nutation in the
              longitude of the ecliptic
       D_EPS: set this to a named variable that will contain the nutation in the
                          obliquity of the ecliptic
    EXAMPLE:
       (1) Example 23a in Meeus: On 2028 Nov 13.19 TD the mean position of Theta
           Persei is 2h 46m 11.331s 49d 20' 54.54".    Determine the shift in 
           position due to the Earth's nutation.
       
           IDL> jd = JULDAY(11,13,2028,.19*24)       ;Get Julian date
           IDL> CO_NUTATE, jd,ten(2,46,11.331)*15.,ten(49,20,54.54),d_ra,d_dec    
    
                 ====> d_ra = 15.843"   d_dec = 6.217"
    PROCEDURES USED:
       NUTATE 
    REVISION HISTORY:
       Written  Chris O'Dell, 2002
       Vector call to NUTATE   W. Landsman   June 2002
  """
  
  ra = np.array(ra)
  dec = np.array(dec)
  num = ra.size
  
  if ra.size != dec.size:
    raise(PE.PyAValError("`ra` and `dec` must have the same size. The current situation is " + \
          "size(ra) = " + str(ra.size) + ", size(dec) = " + str(dec.size) + ".", \
          where="co_nutate", \
          solution="Make arrays the same size or use float as input type."))
    
  # Julian centuries from J2000 of jd.
  jdcen = (np.array(jd) - 2451545.0)/36525.0
  if jdcen.size != ra.size:
    raise(PE.PyAValError("`jd` and `ra` must have the same size. The current situation is " + \
          "size(jd) = " + str(jdcen.size) + ", size(ra) = " + str(ra.size) + ".", \
          where="co_nutate", \
          solution="Make arrays the same size or use float as input type."))

  # Must calculate obliquity of ecliptic
  nut = nutate(jd)
  # Change degrees to seconds
  d_psi = nut[0]*3600.
  d_eps = nut[1]*3600.

  eps0 = 23.4392911*3600. - 46.8150*jdcen - 0.00059*jdcen**2 + 0.001813*jdcen**3
  # True obliquity of the ecliptic in radians
  eps = (eps0 + d_eps)/3600.*(np.pi/180.)

  ce = np.cos(eps)
  se = np.sin(eps)

  # convert ra-dec to equatorial rectangular coordinates
  x = np.cos(ra*np.pi/180.) * np.cos(dec*np.pi/180.)
  y = np.sin(ra*np.pi/180.) * np.cos(dec*np.pi/180.)
  z = np.sin(dec*np.pi/180.)

  # apply corrections to each rectangular coordinate
  x2 = x - (y*ce + z*se)*d_psi * (np.pi/(180.*3600.))
  y2 = y + (x*ce*d_psi - z*d_eps) * (np.pi/(180.*3600.))
  z2 = z + (x*se*d_psi + y*d_eps) * (np.pi/(180.*3600.))
    
  # convert back to equatorial spherical coordinates
  r = np.sqrt(x2**2 + y2**2 + z2**2)
  xyproj = np.sqrt(x2**2 + y2**2)

  ra2 = x2 * 0.
  dec2= x2 * 0.

  if num == 1:
    # Calculate Ra and Dec in RADIANS (later convert to DEGREES)
    if np.logical_and( xyproj == 0 , z != 0 ):
      # Places where xyproj==0 (point at NCP or SCP)
      dec2 = np.arcsin(z2/r)
      ra2 = 0.
    if xyproj != 0:
      # places other than NCP or SCP
      ra2 = np.arctan2( y2 , x2 )
      dec2 = np.arcsin( z2/r )      
  else:
    w1 = np.where( np.logical_and( xyproj == 0 , z != 0 ) )[0]
    w2 = np.where( xyproj != 0 )[0]
    # Calculate Ra and Dec in RADIANS (later convert to DEGREES)
    if len(w1) > 0:
      # Places where xyproj==0 (point at NCP or SCP)
      dec2[w1] = np.arcsin(z2[w1]/r[w1])
      ra2[w1] = 0.
    if len(w2) > 0:
      # Places other than NCP or SCP
      ra2[w2] = np.arctan2( y2[w2] , x2[w2] )
      dec2[w2] = np.arcsin( z2[w2]/r[w2] )

  # Convert into DEGREES
  ra2 = ra2/np.pi*180.
  dec2 = dec2/np.pi*180.
  d_psi /= 3600.
  d_eps /= 3600.

  if num == 1:
    if ra2 < 0.: ra2 += 360.
  else:
    w = np.where( ra2 < 0. )[0]
    if len(w) > 0: ra2[w] += 360.

  # Return changes in ra and dec
  d_ra = (ra2 - ra)
  d_dec = (dec2 - dec)

  if radian:
    # convert result to RADIAN
    d_ra *= (np.pi/180.)
    d_dec *= (np.pi/180.)
    d_psi *= (np.pi/180.)
    d_eps *= (np.pi/180.)
  else:
    eps = eps/np.pi*180. # eps in DEGREES

  if plot:
    if not _ic.check["matplotlib"]:
      raise(PE.PyARequiredImport("Could not import matplotlib.", \
            where="nutate", \
            solution=["Install matplotlib", "Switch `plot` flag to False."]))
    import matplotlib.pylab as plt
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    if jdcen.size > 1:
      plt.plot(jd, d_ra)
    else:
      plt.plot(ra, d_ra)
    ax2 = fig.add_subplot(212)
    if jdcen.size > 1:
      plt.plot(jd, d_dec)
    else:
      plt.plot(dec, d_dec)
    ax2.set_xlabel('Julian date')
    if radian==True:
      ax1.set_ylabel('RA [red]')
      ax2.set_ylabel('DEC [rad]')
    else:
      ax1.set_ylabel('RA [deg]')
      ax2.set_ylabel('DEC [deg]')
    plt.show()

  if full_output:
    return d_ra, d_dec, eps, d_psi, d_eps    
  else:
    return d_ra, d_dec


def co_aberration(jd, ra, dec, radian=False):
  """
    Computes the changes in RA and DEC due to annual aberration.
    
    Parameters
    ----------
    jd : float or array
         The Julian date(s). If array, must be the same size as
         `ra` and `dec`.
    ra : float or array
         The right ascension in degrees.
         If array, it must be the same size as `dec`.
    dec : float or array
         The declination in degrees.
         If array, it must be the same size as `ra`.
    radian : boolean, optional
         Results are returned in radian instead of degrees.
         The default is False.
             
    Returns
    -------
    dRa : float or array
        The change in right ascension [by default in deg].
    dDec : float or arrays
        The change in declination [by default in deg].

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:  

     NAME:
        CO_ABERRATION
    PURPOSE:
        Calculate changes to Ra and Dec due to the effect of annual aberration 
    EXPLANATION:
         as described in Meeus, Chap 23.
    CALLING SEQUENCE:
         co_aberration, jd, ra, dec, d_ra, d_dec, [EPS = ]
    INPUTS
          jd      : Julian Date [scalar or vector]
          ra, dec : Arrays (or scalars) of the ra  and dec's in degrees
      Note: if jd is a vector, then ra and dec must either be scalars, or 
                   vectors of the same length.
    
    OUTPUTS
          d_ra, d_dec: the corrections to ra and dec due to aberration in 
                       arcseconds.  (These values can be added to the true RA 
                       and dec to get the apparent position).   Note that d_ra
                        is *not* multiplied by cos(dec), so that 
                        apparent_ra = ra + d_ra/3600. 
    OPTIONAL INPUT KEYWORD:
          eps : set this to the true obliquity of the ecliptic (in radians), or
            it will be set for you if you don't know it (in that case, set it to
                    an empty variable).
    EXAMPLE:
      Compute the change in RA and Dec of Theta Persei (RA = 2h46m,11.331s, Dec =
      49d20',54.54") due to aberration on 2028 Nov 13.19 TD
    
         IDL> jdcnv,2028,11,13,.19*24,jd      ;Get Julian date
         IDL> co_aberration,jd,ten(2,46,11.331)*15,ten(49,20,54.54),d_ra,d_dec
    
         ==> d_ra = 30.045" (=2.003s)    d_dec = 6.697"
    NOTES:
     These formula are from Meeus, Chapters 23.  Accuracy is much better than 1 
      arcsecond.
    
      The maximum deviation due to annual aberration is 20.49" and occurs when the
      Earth velocity is perpendicular to the direction of the star.
    
    REVISION HISTORY:
      Written, June 2002,      Chris O'Dell, U. of Wisconsin
      Fix error with vector input   W. Landsman   June 2009
      June 2009 update fixed case where JD was scalar but RA,Dec were vectors, but 
      broke the case when both JD and RA,Dec were vectors Aug 2012 W. Landsman
  """

  ra = np.array(ra, ndmin=1)
  dec = np.array(dec, ndmin=1)
  jd = np.array(jd, ndmin=1)
  
  if ra.size != dec.size:
    raise(PE.PyAValError("`ra` and `dec` must have the same size. The current situation is " + \
          "size(ra) = " + str(ra.size) + ", size(dec) = " + str(dec.size) + ".", \
          where="co_aberration", \
          solution="Make arrays the same size or use float as input type."))
    
  if jd.size != ra.size:
    raise(PE.PyAValError("`jd` and `ra` must have the same size. The current situation is " + \
          "size(jd) = " + str(jd.size) + ", size(ra) = " + str(ra.size) + ".", \
          where="co_aberration", \
          solution="Make arrays the same size or use float as input type."))

  # Julian centuries from J2000 of jd.
  jdcen = ( jd - 2451545.0 ) / 36525.0
  # Must calculate obliquity of ecliptic
  res = nutate(jd)
  d_eps = res[1]

  eps0 = 23.4392911*3600. - 46.8150*jdcen - 0.00059*jdcen**2 + 0.001813*jdcen**3
  eps = (eps0/3600. + d_eps)*(np.pi/180.)  #true obliquity of the ecliptic in radians
  
  if jd.size == 1:
    sunlon = np.ravel(sunpos(jd, full_output=True)[3])
  else:
    sunlon = np.zeros( jd.size )
    for i in smo.range(jd.size):
      sunlon[i] = np.ravel(sunpos(jd[i], full_output=True)[3])

  # Earth's orbital eccentricity
  e = 0.016708634 - 0.000042037*jdcen - 0.0000001267*jdcen**2
  # Longitude of perihelion in degrees 
  pi = 102.93735 + 1.71946*jdcen + 0.00046*jdcen**2 
  # Constant of aberration in arcseconds
  k = 20.49552
  
  # Trigonometric Functions
  cd = np.cos(dec*np.pi/180.)
  sd = np.sin(dec*np.pi/180.)
  ce = np.cos(eps)
  te = np.tan(eps)
  cp = np.cos(pi*np.pi/180.)
  sp = np.sin(pi*np.pi/180.)
  cs = np.cos(sunlon*np.pi/180.)
  ss = np.sin(sunlon*np.pi/180.)
  ca = np.cos(ra*np.pi/180.)
  sa = np.sin(ra*np.pi/180.)

  term1 = (ca*cs*ce+sa*ss)/cd
  term2 = (ca*cp*ce+sa*sp)/cd
  term3 = (cs*ce*(te*cd-sa*sd)+ca*sd*ss)
  term4 = (cp*ce*(te*cd-sa*sd)+ca*sd*sp)

  # in ARCSECONDS
  d_ra = -k * term1 + e*k * term2
  d_dec = -k * term3 + e*k * term4
  d_ra /= 3600.
  d_dec /= 3600.

  if radian:
    # Convert result into radian
    d_ra *= (np.pi/180.)
    d_dec *= (np.pi/180.)

  return d_ra, d_dec


def ct2lst(jd, lon):
  """
    Converts the Local Civil Time (given as Julian Date) into Local Mean Sidereal Time.
    
    Parameters
    ----------
    jd : float or array
         The Local Civil Time as Julian date (UT).
         If given as array, its size must be equal to that of `lon`.
    lon : float or array
         The right ascension in DEGREES.
                 
    Returns
    -------
    Time : list
        A list holding the Local Mean Sidereal Time in hours (0 - 24)
        for the given Julian dates and right ascensions.        

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:  

    NAME:
        CT2LST
    PURPOSE:
        To convert from Local Civil Time to Local Mean Sidereal Time.
    
    CALLING SEQUENCE:
        CT2LST, Lst, Lng, Tz, Time, [Day, Mon, Year] 
                          or
        CT2LST, Lst, Lng, dummy, JD
    
    INPUTS:
        Lng  - The longitude in degrees (east of Greenwich) of the place for 
               which the local sidereal time is desired, scalar.   The Greenwich 
               mean sidereal time (GMST) can be found by setting Lng = 0.
        Tz  - The time zone of the site in hours, positive East  of the Greenwich
              meridian (ahead of GMT).  Use this parameter to easily account 
              for Daylight Savings time (e.g. -4=EDT, -5 = EST/CDT), scalar
              This parameter is not needed (and ignored) if Julian date is 
              supplied.    ***Note that the sign of TZ was changed in July 2008
              to match the standard definition.*** 
        Time or JD  - If more than four parameters are specified, then this is 
                  the time of day of the specified date in decimal hours.  If 
                  exactly four parameters are specified, then this is the 
                  Julian date of time in question, scalar or vector
    
    OPTIONAL INPUTS:
         Day -  The day of the month (1-31),integer scalar or vector
         Mon -  The month, in numerical format (1-12), integer scalar or vector
         Year - The 4 digit year (e.g. 2008), integer scalar or vector
    
    OUTPUTS:
          Lst   The Local Sidereal Time for the date/time specified in hours.
    
    RESTRICTIONS:
          If specified, the date should be in numerical form.  The year should
          appear as yyyy.
    
    PROCEDURE:
          The Julian date of the day and time is question is used to determine
          the number of days to have passed since 0 Jan 2000.  This is used
          in conjunction with the GST of that date to extrapolate to the current
          GST; this is then used to get the LST.    See Astronomical Algorithms
          by Jean Meeus, p. 84 (Eq. 11-4) for the constants used.
    
    EXAMPLE:
          Find the Greenwich mean sidereal time (GMST) on 2008 Jul 30 at 15:53 pm
          in Baltimore, Maryland (longitude=-76.72 degrees).   The timezone is 
          EDT or tz=-4
    
          IDL> CT2LST, lst, -76.72, -4,ten(15,53), 30, 07, 2008
    
                  ==> lst =  11.356505  hours  (= 11h 21m 23.418s)
    
          The Web site  http://tycho.usno.navy.mil/sidereal.html contains more
          info on sidereal time, as well as an interactive calculator.
    PROCEDURES USED:
          jdcnv - Convert from year, month, day, hour to julian date
    
    MODIFICATION HISTORY:
        Adapted from the FORTRAN program GETSD by Michael R. Greason, STX, 
                  27 October 1988.
        Use IAU 1984 constants Wayne Landsman, HSTX, April 1995, results 
                  differ by about 0.1 seconds  
        Longitudes measured *east* of Greenwich   W. Landsman    December 1998
        Time zone now measure positive East of Greenwich W. Landsman July 2008
        Remove debugging print statement  W. Landsman April 2009
  """

  jd_vals = np.array(jd)
  lon_vals = np.array(lon)
  
  if jd_vals.size != lon_vals.size:
    raise(PE.PyAValError("`jd` and `lon` need to be of the same length.", \
          where="ct2lst", \
          solution="Make arrays the same length."))
      
  
  # Useful constants, see Meeus, p.84
  c = np.array([280.46061837, 360.98564736629, 0.000387933, 38710000.0 ], float)
  jd2000 = 2451545.0
  t0 = jd_vals - jd2000
  time = t0/36525.0
  
  # Compute GST in seconds.
  theta = c[0] + (c[1] * t0) +  ( c[2] - time / c[3] ) * time**2
  
  # Compute LST in hours.
  lst = ( theta + lon_vals )/15.0
  neg = np.where(lst < 0.0)[0]
  if len(neg) > 0:
    lst[neg] = 24.0 + idlMod(lst[neg], 24.)
  
  # Local sidereal time in hours (0. to 24.)
  lst = idlMod(lst, 24.0)

  return lst


def hadec2altaz(ha, dec, lat, ws=False, radian=False):
  """
    Convert hour angle and declination into horizon (alt/az) coordinates.
    
    Parameters
    ----------
    ha : float or array
        Local apparent hour angle in DEGREES.
    dec : float or array
        Local apparent declination in DEGREES.
    lat : float or array
        Local latitude in DEGREES.
    radian : boolean, optional
        If True, the result is returned in radian
        instead of in degrees (default is False).
    ws : boolean, optional
        Set this to True, if the azimuth shall be measured West from South.
        Default is to measure azimuth East from North.

    Returns
    -------
    Altitude : list
        A list holding the Local Apparent Altitude [deg].
    Apparent Azimuth : list
        The Local Apparent Azimuth [deg].
        
    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:  

     NAME:
        HADEC2ALTAZ
     PURPOSE:
         Converts Hour Angle and Declination to Horizon (alt-az) coordinates.
     EXPLANATION:
         Can deal with NCP/SCP singularity.    Intended mainly to be used by
         program EQ2HOR
    
    CALLING SEQUENCE:
         HADEC2ALTAZ, ha, dec, lat ,alt ,az [ /WS ]
    
    INPUTS
        ha -  the local apparent hour angle, in DEGREES, scalar or vector
        dec -  the local apparent declination, in DEGREES, scalar or vector
        lat -  the local latitude, in DEGREES, scalar or vector
    
    OUTPUTS
        alt - the local apparent altitude, in DEGREES.
        az  - the local apparent azimuth, in DEGREES, all results in double
              precision
    OPTIONAL KEYWORD INPUT:
         /WS - Set this keyword for the output azimuth to be measured West from 
               South.    The default is to measure azimuth East from North.
    
    EXAMPLE:
        What were the apparent altitude and azimuth of the sun when it transited 
        the local meridian at Pine Bluff Observatory (Lat=+43.07833 degrees) on 
        April 21, 2002?   An object transits the local meridian at 0 hour angle.
        Assume this will happen at roughly 1 PM local time (18:00 UTC).
    
        IDL> jdcnv, 2002, 4, 21, 18., jd  ; get rough Julian date to determine 
                                          ;Sun ra, dec.
        IDL> sunpos, jd, ra, dec
        IDL> hadec2altaz, 0., dec, 43.078333, alt, az
    
          ===> Altitude alt = 58.90
               Azimuth  az = 180.0

    REVISION HISTORY:
         Written  Chris O'Dell Univ. of Wisconsin-Madison May 2002
  """
  
  ha = np.array(ha)
  dec = np.array(dec)
  lat = np.array(lat)
  
  if np.logical_or(ha.size != dec.size, dec.size != lat.size):
    raise(PE.PyAValError("`ha`, `dec`, and `lat` must be of the same size. " + \
          "Currently, size(ha) = " + str(ha.size) + ", size(dec) = " + str(dec.size) + ", " + \
          "size(lat) = " + str(lat.size), \
          where="hadec2altaz", \
          solution="Make the arrays the same size."))
  
  sh = np.sin(ha*np.pi/180.)
  ch = np.cos(ha*np.pi/180.)
  sd = np.sin(dec*np.pi/180.)
  cd = np.cos(dec*np.pi/180.)
  sl = np.sin(lat*np.pi/180.)
  cl = np.cos(lat*np.pi/180.)

  x = - ch * cd * sl + sd * cl
  y = - sh * cd
  z = ch * cd * cl + sd * sl
  r = np.sqrt(x**2 + y**2)

  # Now get Alt, Az
  az = np.arctan2(y,x) / (np.pi/180.)
  alt = np.arctan2(z,r) / (np.pi/180.)

  # Correct for negative AZ
  if ha.size==1:
    if az < 0: az += 360.
  else:
    w = np.where(az < 0)[0]
    if len(w) > 0: az[w] += 360.

  # Convert AZ into West from South, if desired
  if ws: az = idlMod( (az + 180.), 360.)
  
  if radian:
    alt *= np.pi/180.
    az *= np.pi/180.

  return alt, az


def co_refract_forward(alt, pressure=1010., temperature=10.0):
  """
    Converts the observed into the apparent (real) altitude.
    
    The *observed altitude* is the altitude that a star is seen
    to be with a telescope. This is where it appears in the sky.
    The observed altitude is always greater than the
    the *apparent altitude*, which is the altitude that a star would
    be at, if there were no atmosphere (sometimes called "true" altitude).
    
    Parameters
    ----------
    alt : float or array
        Observed altitude of an object in DEGREES.
    pressure : float or array, optional
        Atmospheric pressure in MILLIBAR.
        Default pressure is 1010 mbar. If a single value is
        given, it will be used for all given altitudes.
    temperature : float or array, optional
        Ground temperature in degrees Celsius.
        Default temperature is 10 Celsius. If a single value is
        given, it will be used for all given altitudes.
         
    Returns
    -------
    Altitude correction : array
        An array holding the altitude correction [deg]. To convert
        observed altitude into apparent (real) altitude, the
        correction needs to be subtracted from the observed
        altitude.

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:  

    function co_refract_forward, a, P=P, T=T

    INPUTS
       a = The observed (apparent) altitude, in DEGREES.
           May be scalar or vector.
    
    INPUT KEYWORDS
       P:  Pressure [in millibars]. Default is 1010 millibars. [scalar or vector]
       T:  Ground Temp [in Celsius].  Default is 0 Celsius. [scalar or vector]

    .. see also co_refract documentation.
  
  """
  
  alt = np.array(alt, ndmin=1)
  pres = np.array(pressure, ndmin=1)
  # Temperature in Kelvin
  temper = np.array(temperature, ndmin=1) + 273.15
  
  if pres.size != temper.size:
    raise(PE.PyAValError("Pressure and temperature parameters must have the same size. " + \
          "Currently, size(pressure) = " + str(pressure.size) + ", size(temperature) = " + \
          str(temperature.size), \
          where="co_refract_forward", \
          solution="Provide arrays of the same size."))

  if pres.size > 1:
    if pres.size != alt.size:
      raise(PE.PyAValError("Pressure and altitude parameters must have the same size. " + \
            "Currently, size(pressure) = " + str(pressure.size) + ", size(altitude) = " + \
            str(alt.size), \
            where="co_refract_forward", \
              solution="Provide arrays of the same size."))

  # You observed the altitude alt, and would like to know what the "apparent" 
  # altitude is (the one behind the atmosphere).
  R = 0.0166667 / np.tan( (alt + 7.31/(alt+4.4))*np.pi/180. )

  w = np.where(alt < 15.)[0]
  if len(w) > 0:
    R[w] = 3.569*(0.1594 + 0.0196*alt[w] + 0.00002*alt[w]**2)/(1.+0.505*alt[w]+0.0845*alt[w]**2)

  tpcor = pres/1010. * 283./temper
  R *= tpcor

  return R


def co_refract(alt, observer_alt=0.0, pressure=None, temperature=None, epsilon=0.25, \
               convert_to_observed=False, full_output=True):
  """
    Convert between apparent (real) altitude and observed altitude.
    
    This routine converts between the apparent (real) altitude 
    of the object, which does not include the
    influence of the atmosphere, and the observed
    altitude, which is the altitude as seen
    through the atmosphere.
    
    The `convert_to_observed` flag determines the direction of the
    conversion. By default, observed altitude is converted into
    apparent altitude.
    
    Parameters
    ----------
    alt : float or array
        Altitude of an object in DEGREES. Whether the value is
        interpreted as apparent or observed altitude depends on the
        `convert_to_observed` flag. By default, it refers to the
        apparent (real) altitude.
    observer_alt : float or array, optional
        Altitude of the observer in METER. Default is 0 (sea level).
    pressure : float or array, optional
        Atmospheric pressure in MILLIBAR.
        Default pressure is 1010 mbar.
        If `observer_alt` is given, an estimate for the real
        atmospheric pressure is calculated and used.
    temperature : float or array, optional
        Atmospheric temperature at the observing location in Celsius.
        If not specified, the temperature will be calculated
        assuming a ground temperature of 10 degrees Celsius.
    epsilon : float, optional
        If convert_to_observed is TRUE, it specifies the accuracy of
        the calculated altitude in ARCSECONDS that should be reached 
        by the iteration process.
    convert_to_observed : boolean, optional
        If set True, an iterative method is used to calculate
        the observed altitude of an object, which includes
        atmospheric refraction. If False (default), the given altitude
        will be interpreted as the observed altitude and the apparent (real)
        altitude will be calculated using :py:func:`co_refract_forward`.
    full_output : boolean, optional
        If True (default), pressure and temperature used in the calculation
        will be returned as well.

    Returns
    -------
    Altitude : array
        By default, this will be the observed altitude of the object in
        degrees. If `convert_to_observed` was set to False, the number
        refers to the apparent (real) altitude.
    Pressure : array
        The pressure [mbar] used in the calculations (only returned if
        `full_output` is True).
    Temperature : array
        The temperature used in the calculations [K] (only returned if
        `full_output` is True).
        
    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:  

    NAME:
      CO_REFRACT()      
    
    PURPOSE:
      Calculate correction to altitude due to atmospheric refraction.
    
    DESCRIPTION:
      CO_REFRACT can calculate both apparent altitude from observed altitude and 
      vice-versa.
    
    CALLING SEQUENCE:
      new_alt  = CO_REFRACT(old_alt, [ ALTITUDE= , PRESSURE= , $
                                     TEMPERATURE= , /TO_OBSERVED , EPSILON= ])
    
    INPUT:
      old_alt - Observed (apparent) altitude, in DEGREES.  (apparent if keyword 
                /TO_OBSERVED set).    May be scalar or vector.
    
    OUTPUT: 
        Function returns apparent (observed) altitude, in DEGREES. (observed if 
            keyword /TO_OBSERVED set).    Will be of same type as input 
            altitude(s).
    
    OPTIONAL KEYWORD INPUTS:
         ALTITUDE :  The height of the observing location, in meters.  This is 
                only used to determine an approximate temperature and pressure, 
                if these are not specified separately. [default=0, i.e. sea level]
         PRESSURE :  The pressure at the observing location, in millibars.
         TEMPERATURE:    The temperature at the observing location, in Kelvin.
         EPSILON:  When keyword /TO_OBSERVED has been set, this is the accuracy 
                  to  obtain via the iteration, in arcseconds [default = 0.25 
                   arcseconds].
         /TO_OBSERVED:  Set this keyword to go from Apparent->Observed altitude, 
                    using the iterative technique.
    
          Note, if altitude is set, but temperature or pressure are not, the 
          program will make an intelligent guess for the temperature and pressure.
    
    DESCRIPTION:
    
      Because the index of refraction of air is not precisely 1.0, the atmosphere
      bends all incoming light, making a star or other celestial object appear at
      a slightly different altitude (or elevation) than it really is.  It is 
      important to understand the following definitions:
    
      Observed Altitude:  The altitude that a star is SEEN to BE, with a telescope.
                          This is where it appears in the sky.  This is always 
                          GREATER than the apparent altitude.
    
      Apparent Altitude:  The altitude that a star would be at, if *there were no
                        atmosphere* (sometimes called "true" altitude). This is 
                        usually calculated from an object's celestial coordinates.
                        Apparent altitude is always LOWER than the observed 
                        altitude.
    
      Thus, for example, the Sun's apparent altitude when you see it right on the
      horizon is actually -34 arcminutes.
    
      This program uses couple simple formulae to estimate the effect for most 
      optical and radio wavelengths.  Typically, you know your observed altitude 
      (from an observation), and want the apparent altitude.  To go the other way,
      this program uses an iterative approach.
    
    EXAMPLE:
       The lower limb of the Sun is observed to have altitude of 0d 30'.   
       Calculate the the true (=apparent) altitude of the Sun's lower limb using 
       mean  conditions of air pressure and temperature
    
       IDL> print, co_refract(0.5)     ===>  0.025degrees (1.55')
    WAVELENGTH DEPENDENCE:
       This correction is 0 at zenith, about 1 arcminute at 45 degrees, and 34 
       arcminutes at the horizon FOR OPTICAL WAVELENGTHS.  The correction is 
       NON-NEGLIGIBLE at all wavelengths, but is not very easily calculable.  
       These formulae assume a wavelength of 550 nm, and will be accurate to 
       about 4 arcseconds for all visible wavelengths, for elevations of 10 
       degrees and higher.    Amazingly, they are also ACCURATE FOR RADIO 
       FREQUENCIES LESS THAN ~ 100 GHz.
    
       It is important to understand that these formulae really can't do better 
       than about 30 arcseconds of accuracy very close to the horizon, as 
       variable atmospheric effects become very important.
    
    REFERENCES:
       1.  Meeus, Astronomical Algorithms, Chapter 15.
       2.  Explanatory Supplement to the Astronomical Almanac, 1992.
       3.  Methods of Experimental Physics, Vol 12 Part B, Astrophysics, 
           Radio Telescopes, Chapter 2.5, "Refraction Effects in the Neutral 
           Atmosphere", by R.K. Crane.
    
    
    DEPENDENCIES:
       CO_REFRACT_FORWARD (contained in this file and automatically compiled).
    
    AUTHOR:
      Chris O'Dell
          Univ. of Wisconsin-Madison
      Observational Cosmology Laboratory
      Email: odell@cmb.physics.wisc.edu
    
    REVISION HISTORY:
       version 1 (May 31, 2002)
       Update iteration formula,   W. Landsman    June 2002
       Corrected slight bug associated with scalar vs. vector temperature and 
                  pressure inputs. 6/10/2002
       Fixed problem with vector input when /TO_OBSERVED set W. Landsman Dec 2005
       Allow arrays with more than 32767 elements W.Landsman/C.Dickinson Feb 2010
    
  """
  old_alt = np.array(alt, ndmin=1)
    
  if observer_alt is not None:
    # Observer altitude has been given
    obsalt = np.array(observer_alt, ndmin=1)
  else:
    # Assume that observer is on the ground (altitude 0).
    obsalt = np.zeros(old_alt.size)

  if old_alt.size != obsalt.size:
    raise(PE.PyAValError("The number of given observer altitudes (" + str(obsalt.size) + ") " + \
                         "and observed altitudes (" + str(old_alt.size) + ") differ.", \
          where="co_refract", \
          solution="Make the arrays the same size."))

  if temperature is not None:
    # Temperature is given
    temper = np.array(temperature, ndmin=1)
    # Convert temperature to Kelvin
    temper += 273.15
  else:
    # No temperature has been specified. Use approximation
    # based on altitude of observer.
    # 
    # Temperature lapse rate [deg C per meter]
    alpha = 0.0065
    temper = np.zeros(old_alt.size) + 283. - alpha*obsalt
    ind = np.where(obsalt > 11000.)[0]
    if len(ind) > 0:
      temper[ind] = 211.5

  # Estimate Pressure based on altitude, using U.S. Standard Atmosphere formula.
  if pressure is not None:
    # Pressure has been specified
    pres = np.array(pressure, ndmin=1)
  else:
    # Use default atmospheric pressure
    pres = 1010.*(1.0-6.5/288000.0*obsalt)**5.255

  if not convert_to_observed:
    altout = old_alt - co_refract_forward(old_alt, pressure=pres, temperature=temper-273.15)
  else:
    # Convert from real to observed altitude
    altout = np.zeros(old_alt.size)
    for i in smo.range(altout.size):
      dr = co_refract_forward(old_alt[i], pressure=pres[i], temperature=temper[i]-273.15)
      # Guess of observed location
      cur = old_alt[i] + dr
      while True:
        last = cur.copy()
        dr = co_refract_forward(cur, pressure=pres[i], temperature=temper[i]-273.15)
        cur = old_alt[i] + dr
        if np.abs(last - cur)*3600. < epsilon:
          break
      altout[i] = cur.copy()

  if full_output:
    return altout, pres, temper
  else:
    return altout


def eq2hor(jd, ra, dec, observatory=None, lon=None, lat=None, alt=None, B1950=False, precess=True, nutate=True, aberration=True, refract=True):
  """
    Convert celestial coordinates (RA/DEC) to local horizon coordinates (ALT/AZ).
    
    This routine is typically accurate to about 1 arcsec.
    It considers Earth's precession, nutation, aberration,
    and refraction (if keywords are True).
    
    Parameters
    ----------
    jd : float or array
        The Julian date(s)
    ra : float or array
        The right ascension in DEGREES.
    dec : float or array
        The declination in DEGREES.
    observatory : string, {HS}, optional
        A string identifying the observatory. If given,
        the observer's longitude, latitude, and altitude
        are set automatically (and must not be given
        separately then).
    lon : float
         East longitude of the observer in DEGREES.
         Specify West longitude with a negative sign.
         Default is the longitude of the Hamburger Sternwarte.
    lat : float
         Latitude of the observer in DEGREES.
         Default is the latitude of the Hamburger Sternwarte.
    alt : float
         Altitude of the observer in METER.
         Default is the altitude of the Hamburger Sternwarte.
    B1950 : boolean, optional
         If True, your RA and DEC coordinates are given for epoch B1950 FK4.
         If False, RA and DEC are given for epoch J2000 FK5.
         Default is FALSE.
    precess : boolean, optional
         If True (default), Earth's precess motion is considered
         in the calculations.
    nutate : boolean, optional
         If True (default), Earth's nutation is considered
         in the calculations.
    aberration : boolean, optional
         If True (default), the annual aberration is considered
         in the calculations.
    refraction : boolean, optional
         If True, the atmospheric refraction is considered
         in the calculations.
         
    Returns
    -------
    Altitude : float or array
        The altitude in degrees.
    Azimuth : float or array
        The azimuth in degrees (measured East from North).
    Hour angle : float or array
        The hour angle in degrees.

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:  

    NAME:
      EQ2HOR
    
    PURPOSE:
       Convert celestial  (ra-dec) coords to local horizon coords (alt-az).
    
    CALLING SEQUENCE:
    
       eq2hor, ra, dec, jd, alt, az, [ha, LAT= , LON= , /WS, OBSNAME= , $
                          /B1950 , PRECESS_= 0, NUTATE_= 0, REFRACT_= 0, $
                          ABERRATION_= 0, ALTITUDE= , /VERBOSE, _EXTRA= ]
    
    DESCRIPTION:
     This  code calculates horizon (alt,az) coordinates from equatorial
     (ra,dec) coords.   It is typically accurate to about 1 arcsecond or better (I
     have checked the output against the publicly available XEPHEM software). It
     performs precession, nutation, aberration, and refraction corrections.  The
     perhaps best thing about it is that it can take arrays as inputs, in all
     variables and keywords EXCEPT Lat, lon, and Altitude (the code assumes these
     aren't changing), and uses vector arithmetic in every calculation except
     when calculating the precession matrices.
    
    INPUT VARIABLES:
          RA   : Right Ascension of object  (J2000) in degrees (FK5); scalar or
                 vector.
          Dec  : Declination of object (J2000) in degrees (FK5), scalar or vector.
          JD   : Julian Date [scalar or vector]
    
          Note: if RA and DEC are arrays, then alt and az will also be arrays.
                If RA and DEC are arrays, JD may be a scalar OR an array of the
                same dimensionality.
    
    OPTIONAL INPUT KEYWORDS:
          lat   : north geodetic latitude of location in degrees
          lon   : EAST longitude of location in degrees (Specify west longitude
                  with a negative sign.)
          WS    : Set this to get the azimuth measured westward from south (not
                  East of North).
          obsname: Set this to a valid observatory name to be used by the
                 astrolib OBSERVATORY procedure, which will return the latitude
                 and longitude to be used by this program.
          B1950 : Set this if your ra and dec are specified in B1950, FK4
                 coordinates (instead of J2000, FK5)
          precess : Set this to 1 to force precession [default], 0 for no
                  precession correction
          nutate  : Set this to 1 to force nutation [default], 0 for no nutation.
          aberration : Set this to 1 to force aberration correction [default],
                        0 for no correction.
          refract : Set to 1 to force refraction correction [default], 0 for no
                        correction.
          altitude: The altitude of the observing location, in meters. [default=0].
          verbose: Set this for verbose output.  The default is verbose=0.
          extra: This is for setting TEMPERATURE or PRESSURE explicitly, which are
                  used by CO_REFRACT to calculate the refraction effect of the
                  atmosphere. If you don't set these, the program will make an
                  intelligent guess as to what they are (taking into account your
                  altitude).  See CO_REFRACT for more details.
    
    OUTPUT VARIABLES: (all double precision)
          alt    : altitude (in degrees)
          az     : azimuth angle (in degrees, measured EAST from NORTH, but see
                   keyword WS above.)
          ha     : hour angle (in degrees) (optional)
    
    DEPENDENCIES:
          NUTATE, PRECESS, OBSERVATORY, SUNPOS, ADSTRING()
          CO_NUTATE, CO_ABERRATION, CO_REFRACT, ALTAZ2HADEC, SETDEFAULTVALUE
    
    BASIC STEPS
      Apply refraction correction to find apparent Alt.
      Calculate Local Mean Sidereal Time
      Calculate Local Apparent Sidereal Time
      Do Spherical Trig to find apparent hour angle, declination.
      Calculate Right Ascension from hour angle and local sidereal time.
      Nutation Correction to Ra-Dec
      Aberration correction to Ra-Dec
          Precess Ra-Dec to current equinox.
    
    
    ORRECTIONS I DO NOT MAKE:
      *  Deflection of Light by the sun due to GR. (typically milliarcseconds,
           can be arseconds within one degree of the sun)
      *  The Effect of Annual Parallax (typically < 1 arcsecond)
      *  and more (see below)
    
    TO DO
       * Better Refraction Correction.  Need to put in wavelength dependence,
       and integrate through the atmosphere.
           * Topocentric Parallax Correction (will take into account elevation of
             the observatory)
       * Proper Motion (but this will require crazy lookup tables or something).
           * Difference between UTC and UT1 in determining LAST -- is this
             important?
           * Effect of Annual Parallax (is this the same as topocentric Parallax?)
       * Polar Motion
           * Better connection to Julian Date Calculator.
    
    EXAMPLE
    
     Find the position of the open cluster NGC 2264 at the Effelsburg Radio
     Telescope in Germany, on June 11, 2023, at local time 22:00 (METDST).
     The inputs will then be:
    
          Julian Date = 2460107.250
          Latitude = 50d 31m 36s
          Longitude = 06h 51m 18s
          Altitude = 369 meters
          RA (J2000) = 06h 40m 58.2s
          Dec(J2000) = 09d 53m 44.0s
    
     IDL> eq2hor, ten(6,40,58.2)*15., ten(9,53,44), 2460107.250d, alt, az, $
                  lat=ten(50,31,36), lon=ten(6,51,18), altitude=369.0, /verb, $
                   pres=980.0, temp=283.0
    
    The program produces this output (because the VERBOSE keyword was set)
    
    Latitude = +50 31 36.0   Longitude = +06 51 18.0
    Julian Date =  2460107.250000
    Ra, Dec:  06 40 58.2  +09 53 44.0   (J2000)
    Ra, Dec:  06 42 15.7  +09 52 19.2   (J2023.4422)
    Ra, Dec:  06 42 13.8  +09 52 26.9   (fully corrected)
    LMST = +11 46 42.0
    LAST = +11 46 41.4
    Hour Angle = +05 04 27.6  (hh:mm:ss)
    Az, El =  17 42 25.6  +16 25 10.3   (Apparent Coords)
    Az, El =  17 42 25.6  +16 28 22.8   (Observer Coords)
    
    Compare this with the result from XEPHEM:
    Az, El =  17h 42m 25.6s +16d 28m 21s
    
    This 1.8 arcsecond discrepancy in elevation arises primarily from slight
    differences in the way I calculate the refraction correction from XEPHEM, and
    is pretty typical.
    
    AUTHOR:
      Chris O'Dell
          Univ. of Wisconsin-Madison
      Observational Cosmology Laboratory
      Email: odell@cmb.physics.wisc.edu
     Revision History: 
       August 2012  Use Strict_Extra to flag spurious keywords W. Landsman
  """  
  
  if (observatory is not None) and ((lon is not None) or (lat is not None) or \
                                    (alt is not None)):
    raise(PE.PyAParameterConflict("You must specify either `observatory` OR lon, lat, and alt.", \
          solution="Use only one way to specify the observer's location."))
  
  if observatory is not None:
    obsdata = observatory_function(observatory)
    lon = obsdata["longitude"]
    lat = obsdata["latitude"]
    alt = obsdata["altitude"]
  
  if (lon is None) or (lat is None) or (alt is None):
    raise(PE.PyAParameterConflict("You need to specify all of lon, lat, and alt", \
          solution="Specify all parameters"))
  
  jd_vals = np.array(jd, ndmin=1)
  ra_vals = np.array(ra, ndmin=1)
  dec_vals = np.array(dec, ndmin=1)

  if ra_vals.size != dec_vals.size:
    raise(PE.PyAValError("The number of values given for RA and DEC must be equal", \
          where="eq2hor", \
          solution="Make the arrays the same length."))

  if jd_vals.size != ra_vals.size:
    raise(PE.PyAValError("The number of JDs must match the number of coordinates.", \
          where="eq2hor", \
          solution="Make the JD and Ra (or Dec) arrays have the same length."))

  # Precess coordinates to current date
  equinox_now = (jd_vals - 2451545.0)/365.25 + 2000.0
  
  if precess:
    if B1950:
      for i in smo.range(ra_vals.size):
        ra_vals[i], dec_vals[i] = atleg.precess(ra_vals[i], dec_vals[i], 1950.0, equinox_now[i])
    else:
      # Now B2000 is expected
      for i in smo.range(ra_vals.size):
        ra_vals[i], dec_vals[i] = atleg.precess(ra_vals[i], dec_vals[i], 2000.0, equinox_now[i])
        
  # Calculate NUTATION and ABERRATION Corrections to Ra-Dec
  dra1, ddec1, eps, d_psi, dump = co_nutate(jd_vals, ra_vals, dec_vals, full_output=True)
  dra2, ddec2 = co_aberration(jd_vals, ra_vals, dec_vals)

  # Make nutation and aberration corrections
  if nutate:
    ra_vals += dra1
    dec_vals += ddec1
  if aberration:
    ra_vals += dra2
    dec_vals += ddec2

  # Calculate LOCAL MEAN SIDEREAL TIME
  lmst = ct2lst(jd_vals, np.repeat(lon, jd_vals.size))
  # Convert LMST to degrees (btw, this is the RA of the zenith)
  lmst *= 15.0   
  
  # Calculate local APPARENT sidereal time
  # add correction in degrees
  last = lmst + d_psi * np.cos(eps*np.pi/180.)

  # Find hour angle (in DEGREES)
  ha = last - ra_vals
  if ha.size == 1:
    if ha < 0: ha += 360.
  else:
    w = np.where(ha < 0)[0]
    if len(w) > 0: ha[w] += 360.
  
  ha = idlMod(ha, 360.)

  # Now do the spherical trig to get APPARENT alt,az.
  altaz = hadec2altaz(ha, dec_vals, np.ones(ha.size)*lat)

  # Make Correction for ATMOSPHERIC REFRACTION
  # (use this for visible and radio wavelengths; author is unsure about other wavelengths.
  #  See the comments in CO_REFRACT.pro for more details.)
  if refract:
    alts = np.repeat(alt, jd_vals.size)
    altitude = co_refract(altaz[0], observer_alt=alts, full_output=False, convert_to_observed=True)
  else:
    altitude = altaz[0]

  az = altaz[1]
  
  return altitude, az, ha

