from __future__ import print_function, division
import numpy as np
from .astroTimeLegacy import premat, daycnv, precess, helio_jd
from .idlMod import idlMod
from PyAstronomy.pyaC import pyaErrors as PE
import six
import six.moves as smo

def baryvel(dje, deq):
  """
    Calculate helio- and barycentric velocity.
    
    .. note:: The "JPL" option present in IDL is not provided here.
    
    Parameters
    ----------
    dje : float
        Julian ephemeris date
    deq : float
        Epoch of mean equinox of helio- and barycentric velocity output.
        If `deq` is zero, `deq` is assumed to be equal to `dje`.
    
    Returns
    -------
    dvelh : array
        Heliocentric velocity vector [km/s].
    dvelb : array
        Barycentric velocity vector [km/s].

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:
    
    pro baryvel, dje, deq, dvelh, dvelb, JPL = JPL

    NAME:
        BARYVEL  
    PURPOSE:    
        Calculates heliocentric and barycentric velocity components of Earth.

    EXPLANATION:
      BARYVEL takes into account the Earth-Moon motion, and is useful for 
      radial velocity work to an accuracy of  ~1 m/s.

    CALLING SEQUENCE:
      BARYVEL, dje, deq, dvelh, dvelb, [ JPL =  ] 

    INPUTS:
      DJE - (scalar) Julian ephemeris date.
      DEQ - (scalar) epoch of mean equinox of dvelh and dvelb. If deq=0
              then deq is assumed to be equal to dje.
    OUTPUTS: 
      DVELH: (vector(3)) heliocentric velocity component. in km/s 
      DVELB: (vector(3)) barycentric velocity component. in km/s

      The 3-vectors DVELH and DVELB are given in a right-handed coordinate 
      system with the +X axis toward the Vernal Equinox, and +Z axis 
      toward the celestial pole.      

    OPTIONAL KEYWORD SET:
      JPL - if /JPL set, then BARYVEL will call the procedure JPLEPHINTERP
            to compute the Earth velocity using the full JPL ephemeris.   
            The JPL ephemeris FITS file JPLEPH.405 must exist in either the 
            current directory, or in the directory specified by the 
            environment variable ASTRO_DATA.   Alternatively, the JPL keyword
            can be set to the full path and name of the ephemeris file.
            A copy of the JPL ephemeris FITS file is available in
                http://idlastro.gsfc.nasa.gov/ftp/data/         
    PROCEDURES CALLED:
      Function PREMAT() -- computes precession matrix
      JPLEPHREAD, JPLEPHINTERP, TDB2TDT - if /JPL keyword is set
    NOTES:
      Algorithm taken from FORTRAN program of Stumpff (1980, A&A Suppl, 41,1)
      Stumpf claimed an accuracy of 42 cm/s for the velocity.    A 
      comparison with the JPL FORTRAN planetary ephemeris program PLEPH
      found agreement to within about 65 cm/s between 1986 and 1994

      If /JPL is set (using JPLEPH.405 ephemeris file) then velocities are 
      given in the ICRS system; otherwise in the FK4 system.   
    EXAMPLE:
      Compute the radial velocity of the Earth toward Altair on 15-Feb-1994
         using both the original Stumpf algorithm and the JPL ephemeris

      IDL> jdcnv, 1994, 2, 15, 0, jd          ;==> JD = 2449398.5
      IDL> baryvel, jd, 2000, vh, vb          ;Original algorithm
              ==> vh = [-17.07243, -22.81121, -9.889315]  ;Heliocentric km/s
              ==> vb = [-17.08083, -22.80471, -9.886582]  ;Barycentric km/s
      IDL> baryvel, jd, 2000, vh, vb, /jpl   ;JPL ephemeris
              ==> vh = [-17.07236, -22.81126, -9.889419]  ;Heliocentric km/s
              ==> vb = [-17.08083, -22.80484, -9.886409]  ;Barycentric km/s

      IDL> ra = ten(19,50,46.77)*15/!RADEG    ;RA  in radians
      IDL> dec = ten(08,52,3.5)/!RADEG        ;Dec in radians
      IDL> v = vb[0]*cos(dec)*cos(ra) + $   ;Project velocity toward star
              vb[1]*cos(dec)*sin(ra) + vb[2]*sin(dec) 

    REVISION HISTORY:
      Jeff Valenti,  U.C. Berkeley    Translated BARVEL.FOR to IDL.
      W. Landsman, Cleaned up program sent by Chris McCarthy (SfSU) June 1994
      Converted to IDL V5.0   W. Landsman   September 1997
      Added /JPL keyword  W. Landsman   July 2001
      Documentation update W. Landsman Dec 2005
  """
  
  
  # Define constants
  dc2pi = 2 * np.pi 
  cc2pi = 2 * np.pi 
  dc1 = 1.0
  dcto = 2415020.0
  dcjul = 36525.0                     # days in Julian year
  dcbes = 0.313
  dctrop = 365.24219572               # days in tropical year (...572 insig)
  dc1900 = 1900.0
  AU = 1.4959787e8
  
  # Constants dcfel(i,k) of fast changing elements.
  dcfel = [1.7400353e00, 6.2833195099091e02,  5.2796e-6 \
          ,6.2565836e00, 6.2830194572674e02, -2.6180e-6 \
          ,4.7199666e00, 8.3997091449254e03, -1.9780e-5 \
          ,1.9636505e-1, 8.4334662911720e03, -5.6044e-5 \
          ,4.1547339e00, 5.2993466764997e01,  5.8845e-6 \
          ,4.6524223e00, 2.1354275911213e01,  5.6797e-6 \
          ,4.2620486e00, 7.5025342197656e00,  5.5317e-6 \
          ,1.4740694e00, 3.8377331909193e00,  5.6093e-6 ]
  
  dcfel = np.resize(dcfel, (8,3))
  
  # constants dceps and ccsel(i,k) of slowly changing elements.
  dceps = [4.093198e-1, -2.271110e-4, -2.860401e-8 ]
  ccsel = [1.675104e-2, -4.179579e-5, -1.260516e-7 \
          ,2.220221e-1,  2.809917e-2,  1.852532e-5 \
          ,1.589963e00,  3.418075e-2,  1.430200e-5 \
          ,2.994089e00,  2.590824e-2,  4.155840e-6 \
          ,8.155457e-1,  2.486352e-2,  6.836840e-6 \
          ,1.735614e00,  1.763719e-2,  6.370440e-6 \
          ,1.968564e00,  1.524020e-2, -2.517152e-6 \
          ,1.282417e00,  8.703393e-3,  2.289292e-5 \
          ,2.280820e00,  1.918010e-2,  4.484520e-6 \
          ,4.833473e-2,  1.641773e-4, -4.654200e-7 \
          ,5.589232e-2, -3.455092e-4, -7.388560e-7 \
          ,4.634443e-2, -2.658234e-5,  7.757000e-8 \
          ,8.997041e-3,  6.329728e-6, -1.939256e-9 \
          ,2.284178e-2, -9.941590e-5,  6.787400e-8 \
          ,4.350267e-2, -6.839749e-5, -2.714956e-7 \
          ,1.348204e-2,  1.091504e-5,  6.903760e-7 \
          ,3.106570e-2, -1.665665e-4, -1.590188e-7 ]
  ccsel = np.resize(ccsel, (17,3))
  
  # Constants of the arguments of the short-period perturbations.
  dcargs = [5.0974222e0, -7.8604195454652e2 \
           ,3.9584962e0, -5.7533848094674e2 \
           ,1.6338070e0, -1.1506769618935e3 \
           ,2.5487111e0, -3.9302097727326e2 \
           ,4.9255514e0, -5.8849265665348e2 \
           ,1.3363463e0, -5.5076098609303e2 \
           ,1.6072053e0, -5.2237501616674e2 \
           ,1.3629480e0, -1.1790629318198e3 \
           ,5.5657014e0, -1.0977134971135e3 \
           ,5.0708205e0, -1.5774000881978e2 \
           ,3.9318944e0,  5.2963464780000e1 \
           ,4.8989497e0,  3.9809289073258e1 \
           ,1.3097446e0,  7.7540959633708e1 \
           ,3.5147141e0,  7.9618578146517e1 \
           ,3.5413158e0, -5.4868336758022e2 ]
  dcargs = np.resize(dcargs, (15,2))
  
  # Amplitudes ccamps(n,k) of the short-period perturbations.
  ccamps = \
      [-2.279594e-5,  1.407414e-5,  8.273188e-6,  1.340565e-5, -2.490817e-7 \
      ,-3.494537e-5,  2.860401e-7,  1.289448e-7,  1.627237e-5, -1.823138e-7 \
      , 6.593466e-7,  1.322572e-5,  9.258695e-6, -4.674248e-7, -3.646275e-7 \
      , 1.140767e-5, -2.049792e-5, -4.747930e-6, -2.638763e-6, -1.245408e-7 \
      , 9.516893e-6, -2.748894e-6, -1.319381e-6, -4.549908e-6, -1.864821e-7 \
      , 7.310990e-6, -1.924710e-6, -8.772849e-7, -3.334143e-6, -1.745256e-7 \
      ,-2.603449e-6,  7.359472e-6,  3.168357e-6,  1.119056e-6, -1.655307e-7 \
      ,-3.228859e-6,  1.308997e-7,  1.013137e-7,  2.403899e-6, -3.736225e-7 \
      , 3.442177e-7,  2.671323e-6,  1.832858e-6, -2.394688e-7, -3.478444e-7 \
      , 8.702406e-6, -8.421214e-6, -1.372341e-6, -1.455234e-6, -4.998479e-8 \
      ,-1.488378e-6, -1.251789e-5,  5.226868e-7, -2.049301e-7,  0.e0 \
      ,-8.043059e-6, -2.991300e-6,  1.473654e-7, -3.154542e-7,  0.e0 \
      , 3.699128e-6, -3.316126e-6,  2.901257e-7,  3.407826e-7,  0.e0 \
      , 2.550120e-6, -1.241123e-6,  9.901116e-8,  2.210482e-7,  0.e0 \
      ,-6.351059e-7,  2.341650e-6,  1.061492e-6,  2.878231e-7,  0.e0 ]
  ccamps = np.resize(ccamps, (15,5))
  
  # Constants csec3 and ccsec(n,k) of the secular perturbations in longitude.
  ccsec3 = -7.757020e-8
  ccsec = [1.289600e-6, 5.550147e-1, 2.076942e00 \
          ,3.102810e-5, 4.035027e00, 3.525565e-1 \
          ,9.124190e-6, 9.990265e-1, 2.622706e00 \
          ,9.793240e-7, 5.508259e00, 1.559103e01 ]
  ccsec = np.resize(ccsec, (4,3))
  
  # Sidereal rates.
  dcsld = 1.990987e-7                   # sidereal rate in longitude
  ccsgd = 1.990969e-7                   # sidereal rate in mean anomaly
  
  # Constants used in the calculation of the lunar contribution.
  cckm = 3.122140e-5
  ccmld = 2.661699e-6
  ccfdi = 2.399485e-7
  
  # Constants dcargm(i,k) of the arguments of the perturbations of the motion
  # of the moon.
  dcargm = [5.1679830e0,  8.3286911095275e3 \
           ,5.4913150e0, -7.2140632838100e3 \
           ,5.9598530e0,  1.5542754389685e4 ]
  dcargm = np.resize(dcargm, (3,2))
  
  # Amplitudes ccampm(n,k) of the perturbations of the moon.
  ccampm = [ 1.097594e-1, 2.896773e-7, 5.450474e-2,  1.438491e-7 \
           ,-2.223581e-2, 5.083103e-8, 1.002548e-2, -2.291823e-8 \
           , 1.148966e-2, 5.658888e-8, 8.249439e-3,  4.063015e-8 ]
  ccampm = np.resize(ccampm, (3,4))
  
  # ccpamv(k)=a*m*dl,dt (planets), dc1mme=1-mass(earth+moon)
  ccpamv = [8.326827e-11, 1.843484e-11, 1.988712e-12, 1.881276e-12]
  dc1mme = 0.99999696e0
  
  # Time arguments.
  dt = (dje - dcto) / dcjul
  tvec = np.array([1e0, dt, dt*dt])
  
  # Values of all elements for the instant(aneous?) dje.
  temp = idlMod(np.dot(dcfel, tvec), dc2pi)

  dml = temp[0]
  forbel = temp[1:8]
  g = forbel[0]                         # old fortran equivalence
  
  deps = idlMod(np.sum(tvec*dceps), dc2pi)
  sorbel = idlMod(np.dot(ccsel, tvec), dc2pi)
  e = sorbel[0]                         # old fortran equivalence
  
  # Secular perturbations in longitude.
  dummy = np.cos(2.0)
  sn = np.sin(idlMod(np.dot(ccsec[::,1:3], tvec[0:2]), cc2pi))
  
  # Periodic perturbations of the emb (earth-moon barycenter).
  pertl = np.sum(ccsec[::,0] * sn) + (dt * ccsec3 * sn[2])
  pertld = 0.0
  pertr = 0.0
  pertrd = 0.0
  for k in smo.range(15):
    a = idlMod((dcargs[k,0] + dt*dcargs[k,1]), dc2pi)
    cosa = np.cos(a)
    sina = np.sin(a)
    pertl = pertl + ccamps[k,0]*cosa + ccamps[k,1]*sina
    pertr = pertr + ccamps[k,2]*cosa + ccamps[k,3]*sina
    if k < 11:
      pertld = pertld + (ccamps[k,1]*cosa-ccamps[k,0]*sina)*ccamps[k,4]
      pertrd = pertrd + (ccamps[k,3]*cosa-ccamps[k,2]*sina)*ccamps[k,4]
  
  # Elliptic part of the motion of the emb.
  phi = (e*e/4e0)*(((8e0/e)-e)*np.sin(g) +5*np.sin(2*g) +(13/3e0)*e*np.sin(3*g))
  f = g + phi
  sinf = np.sin(f)
  cosf = np.cos(f)
  dpsi = (dc1 - e*e) / (dc1 + e*cosf)
  phid = 2*e*ccsgd*((1 + 1.5*e*e)*cosf + e*(1.25 - 0.5*sinf*sinf))
  psid = ccsgd*e*sinf / np.sqrt(dc1 - e*e)
  
  # Perturbed heliocentric motion of the emb.
  d1pdro = dc1+pertr
  drd = d1pdro * (psid + dpsi*pertrd)
  drld = d1pdro*dpsi * (dcsld+phid+pertld)
  dtl = idlMod((dml + phi + pertl), dc2pi)
  dsinls = np.sin(dtl)
  dcosls = np.cos(dtl)
  dxhd = drd*dcosls - drld*dsinls
  dyhd = drd*dsinls + drld*dcosls
  
  # Influence of eccentricity, evection and variation on the geocentric
  # motion of the moon.
  pertl = 0.0
  pertld = 0.0
  pertp = 0.0
  pertpd = 0.0
  for k in smo.range(3):
    a = idlMod((dcargm[k,0] + dt*dcargm[k,1]), dc2pi)
    sina = np.sin(a)
    cosa = np.cos(a)
    pertl = pertl + ccampm[k,0]*sina
    pertld = pertld + ccampm[k,1]*cosa
    pertp = pertp + ccampm[k,2]*cosa
    pertpd = pertpd - ccampm[k,3]*sina
  
  # Heliocentric motion of the earth.
  tl = forbel[1] + pertl
  sinlm = np.sin(tl)
  coslm = np.cos(tl)
  sigma = cckm / (1.0 + pertp)
  a = sigma*(ccmld + pertld)
  b = sigma*pertpd
  dxhd = dxhd + a*sinlm + b*coslm
  dyhd = dyhd - a*coslm + b*sinlm
  dzhd= -sigma*ccfdi*np.cos(forbel[2])
  
  # Barycentric motion of the earth.
  dxbd = dxhd*dc1mme
  dybd = dyhd*dc1mme
  dzbd = dzhd*dc1mme
  for k in smo.range(4):
    plon = forbel[k+3]
    pomg = sorbel[k+1]
    pecc = sorbel[k+9]
    tl = idlMod((plon + 2.0*pecc*np.sin(plon-pomg)), cc2pi)
    dxbd = dxbd + ccpamv[k]*(np.sin(tl) + pecc*np.sin(pomg))
    dybd = dybd - ccpamv[k]*(np.cos(tl) + pecc*np.cos(pomg))
    dzbd = dzbd - ccpamv[k]*sorbel[k+13]*np.cos(plon - sorbel[k+5])
  
  # Transition to mean equator of date.
  dcosep = np.cos(deps)
  dsinep = np.sin(deps)
  dyahd = dcosep*dyhd - dsinep*dzhd
  dzahd = dsinep*dyhd + dcosep*dzhd
  dyabd = dcosep*dybd - dsinep*dzbd
  dzabd = dsinep*dybd + dcosep*dzbd
  
  # Epoch of mean equinox (deq) of zero implies that we should use
  # Julian ephemeris date (dje) as epoch of mean equinox.
  if deq == 0:
    dvelh = AU * np.array([dxhd, dyahd, dzahd])
    dvelb = AU * np.array([dxbd, dyabd, dzabd])
    return dvelh, dvelb
  
  # General precession from epoch dje to deq.
  deqdat = (dje-dcto-dcbes) / dctrop + dc1900
  prema = np.transpose(premat(deqdat, deq, FK4=True))
  
  dvelh = AU * np.dot( [dxhd, dyahd, dzahd], prema )
  dvelb = AU * np.dot( [dxbd, dyabd, dzabd], prema )
  return dvelh, dvelb


def baryCorr(jd, ra, dec, deq=0.0):
  """
    Calculate barycentric correction.
    
    This function uses the :py:func:`baryvel` function
    to calculate the helio- and barycentric motion of
    the Earth and projects it onto the direction to
    the star.
    
    .. note:: Positive return values indicate that the Earth
              moves toward the star.
    
    Parameters
    ----------
    jd : float
        The time at which to calculate the correction.
    ra : float
        Right ascension in degrees.
    dec : float
        Declination in degrees.
    deq : float, optional
        The mean equinox of barycentric velocity
        calculation (see :py:func:`bryvel`). If zero,
        it is assumed to be the same as `jd`.
    
    Returns
    -------
    Projected heliocentric velocity : float
        Heliocentric velocity toward star [km/s]
    Projected barycentric velocity : float
        Barycentric velocity toward star [km/s]
  """
  vh, vb = baryvel(jd, deq)
  rarad = ra/180. * np.pi
  derad = dec/180. * np.pi
 
  # Projected velocity toward star
  vbp = vb[0]*np.cos(derad)*np.cos(rarad) + \
        vb[1]*np.cos(derad)*np.sin(rarad) + vb[2]*np.sin(derad)
  vhp = vh[0]*np.cos(derad)*np.cos(rarad) + \
        vh[1]*np.cos(derad)*np.sin(rarad) + vh[2]*np.sin(derad) 
  
  return vhp, vbp


def helcorr(obs_long, obs_lat, obs_alt, ra2000, dec2000, jd, debug=False):
  """
    Calculate barycentric velocity correction.
    
    This function calculates the motion of an observer in
    the direction of a star. In contract to :py:func:`baryvel`
    and :py:func:`baryCorr`, the rotation of the Earth is
    taken into account.
    
    .. note:: This function was ported from the REDUCE IDL package.
              See Piskunov & Valenti 2002, A&A 385, 1095 for a detailed
              description of the package and/or visit
              http://www.astro.uu.se/~piskunov/RESEARCH/REDUCE/
    
    .. warning:: Contrary to the original implementation the longitude
                 increases toward the East and the right ascension is
                 given in degrees instead of hours. The JD is given as is,
                 in particular, nothing needs to be subtracted.
    
    Parameters
    ----------
    obs_long : float
        Longitude of observatory (degrees, **eastern** direction is positive)
    obs_lat : float
        Latitude of observatory [deg]
    obs_alt : float
        Altitude of observatory [m]
    ra2000 : float
        Right ascension of object for epoch 2000.0 [deg]
    dec2000 : float
        Declination of object for epoch 2000.0 [deg]
    jd : float
        Julian date for the middle of exposure.
    
    Returns
    -------
    Barycentric correction : float
        The barycentric correction accounting for the rotation
        of the Earth, the rotation of the Earth's center around
        the Earth-Moon barycenter, and the motion of the Earth-Moon 
        barycenter around the center of the Sun [km/s].
    HJD : float
        Heliocentric Julian date for middle of exposure.

    Notes
    -----

    :IDL REDUCE - Documentation:


    Calculates heliocentric Julian date, barycentric and heliocentric radial
    velocity corrections from:
    
    INPUT:
    <OBSLON> Longitude of observatory (degrees, western direction is positive)
    <OBSLAT> Latitude of observatory (degrees)
    <OBSALT> Altitude of observatory (meters)
    <RA2000> Right ascension of object for epoch 2000.0 (hours)
    <DE2000> Declination of object for epoch 2000.0 (degrees)
    <JD> Julian date for the middle of exposure
    [DEBUG=] set keyword to get additional results for debugging
    
    OUTPUT:
    <CORRECTION> barycentric correction - correction for rotation of earth,
       rotation of earth center about the earth-moon barycenter, earth-moon 
       barycenter about the center of the Sun.
    <HJD> Heliocentric Julian date for middle of exposure
    
    Algorithms used are taken from the IRAF task noao.astutils.rvcorrect
    and some procedures of the IDL Astrolib are used as well.
    Accuracy is about 0.5 seconds in time and about 1 m/s in velocity.
    
    History:
    written by Peter Mittermayer, Nov 8,2003
    2005-January-13   Kudryavtsev   Made more accurate calculation of the sidereal time.
                                    Conformity with MIDAS compute/barycorr is checked.
    2005-June-20      Kochukhov Included precession of RA2000 and DEC2000 to current epoch

"""
  from PyAstronomy.pyaC import degtorad

  # This reverts the original longitude convention. After this,
  # East longitudes are positive
  obs_long = -obs_long

  if jd < 2.4e6:
    PE.warn(PE.PyAValError("The given Julian Date (" + str(jd) + ") is exceedingly small. Did you subtract 2.4e6?"))

  # Covert JD to Gregorian calendar date
  xjd = jd
  
  year, month, day, ut = tuple(daycnv(xjd))

  # Current epoch
  epoch = year + month/12. + day/365.

  # Precess ra2000 and dec2000 to current epoch, resulting ra is in degrees
  ra = ra2000
  dec = dec2000
  ra, dec = precess(ra, dec, 2000.0, epoch)  

  # Calculate heliocentric julian date
  rjd = jd-2.4e6
  hjd = helio_jd(rjd, ra, dec) + 2.4e6

  # DIURNAL VELOCITY (see IRAF task noao.astutil.rvcorrect)
  # convert geodetic latitude into geocentric latitude to correct
  # for rotation of earth
  dlat = -(11.*60.+32.743)*np.sin(2.0*degtorad(obs_lat)) \
         +1.1633*np.sin(4.0*degtorad(obs_lat)) - 0.0026*np.sin(6.0*degtorad(obs_lat))
  lat = obs_lat + dlat/3600.0

  # Calculate distance of observer from earth center
  r = 6378160.0 * (0.998327073+0.001676438*np.cos(2.0*degtorad(lat)) \
     -0.00000351 * np.cos(4.0*degtorad(lat)) + 0.000000008*np.cos(6.0*degtorad(lat))) \
     + obs_alt

  # Calculate rotational velocity (perpendicular to the radius vector) in km/s
  # 23.934469591229 is the sidereal day in hours for 1986
  v = 2.*np.pi * (r/1000.) / (23.934469591229*3600.)

  # Calculating local mean sidereal time (see astronomical almanach)
  tu = (rjd-51545.0)/36525.0
  gmst = 6.697374558 + ut + \
        (236.555367908*(rjd-51545.0) + 0.093104*tu**2 - 6.2e-6*tu**3)/3600.0
  lmst = idlMod(gmst-obs_long/15., 24)

  # Projection of rotational velocity along the line of sight
  vdiurnal = v*np.cos(degtorad(lat))*np.cos(degtorad(dec))*np.sin(degtorad(ra-lmst*15))

  # BARICENTRIC and HELIOCENTRIC VELOCITIES
  vh, vb = baryvel(xjd,0)

  # Project to line of sight
  vbar = vb[0]*np.cos(degtorad(dec))*np.cos(degtorad(ra)) + vb[1]*np.cos(degtorad(dec))*np.sin(degtorad(ra)) + \
         vb[2]*np.sin(degtorad(dec))
  vhel = vh[0]*np.cos(degtorad(dec))*np.cos(degtorad(ra)) + vh[1]*np.cos(degtorad(dec))*np.sin(degtorad(ra)) + \
         vh[2]*np.sin(degtorad(dec))
  
  # Use barycentric velocity for correction
  corr = (vdiurnal + vbar) 

  if debug:
    print('')
    print('----- HELCORR.PRO - DEBUG INFO - START ----')
    print('(obs_long (East positive),obs_lat,obs_alt) Observatory coordinates [deg,m]: ', -obs_long, obs_lat, obs_alt)
    print('(ra,dec) Object coordinates (for epoch 2000.0) [deg]: ', ra,dec)
    print('(ut) Universal time (middle of exposure) [hrs]: ', ut)
    print('(jd) Julian date (middle of exposure) (JD): ', jd)
    print('(hjd) Heliocentric Julian date (middle of exposure) (HJD): ', hjd)
    print('(gmst) Greenwich mean sidereal time [hrs]: ', idlMod(gmst, 24))
    print('(lmst) Local mean sidereal time [hrs]: ', lmst)
    print('(dlat) Latitude correction [deg]: ', dlat)
    print('(lat) Geocentric latitude of observer [deg]: ', lat)
    print('(r) Distance of observer from center of earth [m]: ', r)
    print('(v) Rotational velocity of earth at the position of the observer [km/s]: ', v)
    print('(vdiurnal) Projected earth rotation and earth-moon revolution [km/s]: ', vdiurnal)
    print('(vbar) Barycentric velocity [km/s]: ', vbar)
    print('(vhel) Heliocentric velocity [km/s]: ', vhel)
    print('(corr) Vdiurnal+vbar [km/s]: ', corr)
    print('----- HELCORR.PRO - DEBUG INFO - END -----')
    print('')
  
  return corr, hjd
