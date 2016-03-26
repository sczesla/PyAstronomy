# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
from .idlCirrange import cirrange
from .eq2hor import nutate
import six.moves as smo


def moonpos(jd, radian=False):
  """
    Computes RA and DEC of the Moon at given Julian date(s).
    
    Parameters
    ----------
    jd : float or array
        The Julian date.
    radian : boolean, optional
        If True, results are returned in RADIAN instead
        of DEGREES. Default is False.
    
    Returns
    -------
    RA : float or array
        Right ascension of the moon for given JD(s) in DEGREES.
    DEC : float or array
        Declination of the moon for given JD(s) in DEGREES.
    DISTANCE : float or array
        Distance of the moon from the Earth for given JD(s) in KILOMETERS.
    GEOLONGITUDE : float or array
        Apparent longitude of the moon for given JD(s) in DEGREES.
    GEOLATITUDE : float or array
        Apparent latitude of the moon for given JD(s) in DEGREES.

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:  

    PRO MOONPOS, jd, ra, dec, dis, geolong, geolat, RADIAN = radian
    
    NAME:                                     
          MOONPOS
    PURPOSE:
          To compute the RA and Dec of the Moon at specified Julian date(s).
    
    CALLING SEQUENCE:
          MOONPOS, jd, ra, dec, dis, geolong, geolat, [/RADIAN ]
    
    INPUTS:
          JD - Julian ephemeris date, scalar or vector, double precision suggested
    
    OUTPUTS:
          Ra  - Apparent right ascension of the moon in DEGREES, referred to the
                  true equator of the specified date(s) 
          Dec - The declination of the moon in DEGREES 
          Dis - The Earth-moon distance in kilometers (between the center of the
                Earth and the center of the Moon).
          Geolong - Apparent longitude of the moon in DEGREES, referred to the
                  ecliptic of the specified date(s)
          Geolat - Apparent longitude of the moon in DEGREES, referred to the
                  ecliptic of the specified date(s)
    
          The output variables will all have the same number of elements as the
          input Julian date vector, JD.   If JD is a scalar then the output 
          variables will be also.
    
    OPTIONAL INPUT KEYWORD:
          /RADIAN - If this keyword is set and non-zero, then all output variables 
                  are given in Radians rather than Degrees
    
    EXAMPLES:
          (1) Find the position of the moon on April 12, 1992
    
          IDL> jdcnv,1992,4,12,0,jd    ;Get Julian date
          IDL> moonpos, jd, ra ,dec     ;Get RA and Dec of moon
          IDL> print,adstring(ra,dec,1)
                  ==> 08 58 45.23  +13 46  6.1
    
          This is within 1" from the position given in the Astronomical Almanac
          
          (2) Plot the Earth-moon distance for every day at 0 TD in July, 1996
    
          IDL> jdcnv,1996,7,1,0,jd                   ;Get Julian date of July 1
          IDL> moonpos,jd+dindgen(31), ra, dec, dis  ;Position at all 31 days
          IDL> plot,indgen(31),dis, /YNOZ
    
    METHOD:
          Derived from the Chapront ELP2000/82 Lunar Theory (Chapront-Touze' and
          Chapront, 1983, 124, 50), as described by Jean Meeus in Chapter 47 of
          ``Astronomical Algorithms'' (Willmann-Bell, Richmond), 2nd edition, 
          1998.    Meeus quotes an approximate accuracy of 10" in longitude and
          4" in latitude, but he does not give the time range for this accuracy.
    
          Comparison of this IDL procedure with the example in ``Astronomical
          Algorithms'' reveals a very small discrepancy (~1 km) in the distance 
          computation, but no difference in the position calculation.
    
          This procedure underwent a major rewrite in June 1996, and the new
          calling sequence is *incompatible with the old* (e.g. angles now 
          returned in degrees instead of radians).
    
    PROCEDURES CALLED:
          CIRRANGE, ISARRAY(), NUTATE, TEN()  - from IDL Astronomy Library
          POLY() - from IDL User's Library
    MODIFICATION HISTORY:
          Written by Michael R. Greason, STX, 31 October 1988.
          Major rewrite, new (incompatible) calling sequence, much improved 
                  accuracy,       W. Landsman   Hughes STX      June 1996
          Added /RADIAN keyword  W. Landsman August 1997
          Converted to IDL V5.0   W. Landsman   September 1997
          Use improved expressions for L',D,M,M', and F given in 2nd edition of
               Meeus (very slight change),  W. Landsman    November 2000
          Avoid 32767 overflow   W. Landsman January 2005
  """

  jd = np.array(jd, ndmin=1)
  time = (jd - 2451545.0)/36525.0
  
  d_lng = np.array([0,2,2,0,0,0,2,2,2,2,0,1,0,2,0,0,4,0,4,2,2,1,1,2,2,4,2,0,2,2,1,2,0,0, \
           2,2,2,4,0,3,2,4,0,2,2,2,4,0,4,1,2,0,1,3,4,2,0,1,2,2])

  m_lng = np.array([0,0,0,0,1,0,0,-1,0,-1,1,0,1,0,0,0,0,0,0,1,1,0,1,-1,0,0,0,1,0,-1,0, \
           -2,1,2,-2,0,0,-1,0,0,1,-1,2,2,1,-1,0,0,-1,0,1,0,1,0,0,-1,2,1,0,0])

  mp_lng = np.array([1,-1,0,2,0,0,-2,-1,1,0,-1,0,1,0,1,1,-1,3,-2,-1,0,-1,0,1,2,0,-3,-2, \
            -1,-2,1,0,2,0,-1,1,0,-1,2,-1,1,-2,-1,-1,-2,0,1,4,0,-2,0,2,1,-2,-3,2,1,-1,3,-1])

  f_lng = np.array([0,0,0,0,0,2,0,0,0,0,0,0,0,-2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0, \
           0,0,0,-2,2,0,2,0,0,0,0,0,0,-2,0,0,0,0,-2,-2,0,0,0,0,0,0,0,-2])

  sin_lng = np.array([6288774,1274027,658314,213618,-185116,-114332,58793,57066,53322, \
             45758,-40923,-34720,-30383,15327,-12528,10980,10675,10034,8548,-7888,-6766, \
             -5163,4987,4036,3994,3861,3665,-2689,-2602,2390,-2348,2236,-2120,-2069,2048, \
             -1773,-1595,1215,-1110,-892,-810,759,-713,-700,691,596,549,537,520,-487, \
             -399,-381,351,-340,330,327,-323,299,294,0.0])

  cos_lng = np.array([-20905355,-3699111,-2955968,-569925,48888,-3149,246158,-152138, \
             -170733,-204586,-129620,108743,104755,10321,0,79661,-34782,-23210,-21636, \
             24208,30824,-8379,-16675,-12831,-10445,-11650,14403,-7003,0,10056,6322, \
             -9884,5751,0,-4950,4130,0,-3958,0,3258,2616,-1897,-2117,2354,0,0,-1423, \
             -1117,-1571,-1739,0,-4421,0,0,0,0,1165,0,0,8752.0])

  d_lat = np.array([0,0,0,2,2,2,2,0,2,0,2,2,2,2,2,2,2,0,4,0,0,0,1,0,0,0,1,0,4,4,0,4,2,2, \
                    2,2,0,2,2,2,2,4,2,2,0,2,1,1,0,2,1,2,0,4,4,1,4,1,4,2])

  m_lat = np.array([0,0,0,0,0,0,0,0,0,0,-1,0,0,1,-1,-1,-1,1,0,1,0,1,0,1,1,1,0,0,0,0,0,0, \
                    0,0,-1,0,0,0,0,1,1,0,-1,-2,0,1,1,1,1,1,0,-1,1,0,-1,0,0,0,-1,-2])

  mp_lat = np.array([0,1,1,0,-1,-1,0,2,1,2,0,-2,1,0,-1,0,-1,-1,-1,0,0,-1,0,1,1,0,0,3,0, \
                    -1,1, -2,0,2,1,-2,3,2,-3,-1,0,0,1,0,1,1,0,0,-2,-1,1,-2,2,-2,-1,1,1,-1,0,0])

  f_lat = np.array([ 1,1,-1,-1,1,-1,1,1,-1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,3,1,1,1,-1,-1,-1, \
                     1,-1,1,-3,1,-3,-1,-1,1,-1,1,-1,1,1,1,1,-1,3,-1,-1,1,-1,-1,1,-1,1,-1,-1, \
                    -1,-1,-1,-1,1])

  sin_lat = np.array([5128122,280602,277693,173237,55413,46271,32573,17198,9266,8822, \
                      8216,4324,4200,-3359,2463,2211,2065,-1870,1828,-1794,-1749,-1565,-1491, \
                      -1475,-1410,-1344,-1335,1107,1021,833,777,671,607,596,491,-451,439,422, \
                      421,-366,-351,331,315,302,-283,-229,223,223,-220,-220,-185,181,-177,176, \
                      166,-164,132,-119,115,107.0])

  # Mean longitude of the moon referred to mean equinox of the date
  coeff0 = [-1.0/6.5194e7, 1.0/538841.0, -0.0015786, 481267.88123421, 218.3164477]
  lprimed = np.polyval(coeff0, time)*np.pi/180.
  lprimed = cirrange(lprimed, radians=True)
  
  # Mean elongation of the Moon
  coeff1 = [-1.0/1.13065e8, 1.0/545868.0, -0.0018819, 445267.1114034, 297.8501921]
  d = np.polyval(coeff1, time)*np.pi/180.
  d = cirrange(d, radians=True)

  # Sun's mean anomaly
  coeff2 = [1.0/2.449e7, -0.0001536, 35999.0502909, 357.5291092]
  M = np.polyval(coeff2, time)*np.pi/180.
  M = cirrange(M, radians=True)
  
  # Moon's mean anomaly
  coeff3 = [-1.0/1.4712e7, 1.0/6.9699e4, 0.0087414, 477198.8675055, 134.9633964]
  Mprime = np.polyval(coeff3, time)*np.pi/180.
  Mprime = cirrange(Mprime, radians=True)

  # Moon's argument of latitude
  coeff4 = [1.0/8.6331e8, -1.0/3.526e7, -0.0036539, 483202.0175233, 93.2720950]
  F = np.polyval(coeff4, time)*np.pi/180.
  F = cirrange(F, radians=True)

  # Eccentricity of Earth's orbit around the Sun
  E = 1 - 0.002516*time - 7.4e-6*time**2
  E2 = E**2
  ecorr1 = np.where(np.abs(m_lng) == 1)[0]
  ecorr2 = np.where(np.abs(m_lat) == 1)[0]
  ecorr3 = np.where(np.abs(m_lng) == 2)[0]
  ecorr4 = np.where(np.abs(m_lat) == 2)[0]

  # Additional arguments
  A1 = (119.75 + 131.849*time) * np.pi/180.
  A2 = (53.09 + 479264.290*time) * np.pi/180.
  A3 = (313.45 + 481266.484*time) * np.pi/180.
  suml_add = 3958.*np.sin(A1) + 1962.*np.sin(lprimed - F) + 318.*np.sin(A2)
  sumb_add =  -2235.*np.sin(lprimed) + 382.*np.sin(A3) + 175.*np.sin(A1-F) + \
              175.*np.sin(A1 + F) + 127.*np.sin(lprimed - Mprime) - 115.*np.sin(lprimed + Mprime)

  # Sum the periodic terms 
  geolong = np.zeros(jd.size)
  geolat = np.zeros(jd.size)
  dis = np.zeros(jd.size)

  for i in smo.range(jd.size):
    sinlng = sin_lng.copy()
    coslng = cos_lng.copy()
    sinlat = sin_lat.copy()

    sinlng[ecorr1] = E[i]*sinlng[ecorr1]
    coslng[ecorr1] = E[i]*coslng[ecorr1]
    sinlat[ecorr2] = E[i]*sinlat[ecorr2]
    sinlng[ecorr3] = E2[i]*sinlng[ecorr3]
    coslng[ecorr3] = E2[i]*coslng[ecorr3]
    sinlat[ecorr4] = E2[i]*sinlat[ecorr4]

    arg = d_lng*d[i] + m_lng*M[i] +mp_lng*Mprime[i] + f_lng*F[i]
    geolong[i] = lprimed[i]/(np.pi/180.) + ( np.sum( sinlng*np.sin(arg) ) + suml_add[i] )/1.0e6

    dis[i] = 385000.56 + np.sum( coslng*np.cos(arg) )/1.0e3

    arg = d_lat*d[i] + m_lat*M[i] +mp_lat*Mprime[i] + f_lat*F[i]
    geolat[i] = ( np.sum( sinlat*np.sin(arg) ) + sumb_add[i] )/1.0e6

  # Find the nutation in longitude
  nut = nutate(jd)
  nlong = nut[0]
  elong = nut[1]
  geolong += nlong
  geolong = geolong*np.pi/180.    
  geolong = cirrange(geolong, radians=True)    
  
  lamb = geolong.copy()
  beta = geolat*np.pi/180.
  geolong = geolong/(np.pi/180.)
  
  # Find mean obliquity and convert lambda,beta to RA, Dec 
  c = np.array([21.448,-4680.93,-1.55,1999.25,-51.38,-249.67,-39.05,7.12,27.87,5.79,2.45])
  c = c[::-1]  
  epsilon = 23. + 26./60. + np.polyval(c, time/1.0e2)/3600.0
  eps = ( epsilon + elong )*np.pi/180. # True obliquity in radians

  ra = np.arctan2( np.sin(lamb)*np.cos(eps) - np.tan(beta)*np.sin(eps), np.cos(lamb) )
  ra = cirrange(ra, radians=True) 
  dec = np.arcsin( np.sin(beta)*np.cos(eps) + np.cos(beta)*np.sin(eps)*np.sin(lamb) )

  if not radian:
    ra = ra/(np.pi/180.)
    dec = dec/(np.pi/180.)
  else:
    geolong = lamb
    geolat = beta

  return ra, dec, dis, geolong, geolat
