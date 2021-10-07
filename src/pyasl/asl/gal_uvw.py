# -*- coding: utf-8 -*-
import numpy as np

def gal_uvw(ra, dec, pmra, pmdec, distance, vrad, lsr=None):
    """
    Calculate the Galactic space velocity (U,V,W)
    
    .. note::  In contrast to the IDL implementation, U is taken positive toward the
               galactic center.
    
    Parameters
    ----------
    ra, dec : float or array
        Position in degrees
    pmra, pmdec : float or array
        Proper motion in mas/yr
    distance : float or array
        Distance in pc
    vrad : float or array
        Radial velocity in km/s
    lsr : three-float tuple, optional
        If given, the output UVW velocities are corrected for the
        solar motion (U,V,W)_Sun to the local standard of rest (LSR).
        Use, e.g., lsr = (8.5, 13.38, 6.49) from Coșkunoǧlu et al. 2011,
        MNRAS 412, 1237-1245. Be aware that the value of the solar motion
        through the LSR remains poorly determined.

    Returns
    -------
    U, V, W : float or arrays
        U: Velocity (km/s) positive toward the Galactic center
        V: Velocity (km/s) positive in the direction of Galactic rotation
        W: Velocity (km/s) positive toward the North Galactic Pole
    
    
    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.
    :IDL - Documentation:
    
    NAME:
        GAL_UVW
    PURPOSE:
        Calculate the Galactic space velocity (U,V,W) of star
    EXPLANATION:
        Calculates the Galactic space velocity U, V, W of star given its
        (1) coordinates, (2) proper motion, (3) distance (or parallax), and
        (4) radial velocity.
    CALLING SEQUENCE:
        GAL_UVW [/LSR, RA=, DEC=, PMRA= ,PMDEC=, VRAD= , DISTANCE=
                    PLX= ]
    OUTPUT PARAMETERS:
            U - Velocity (km/s) positive toward the Galactic *anti*center
            V - Velocity (km/s) positive in the direction of Galactic rotation
            W - Velocity (km/s) positive toward the North Galactic Pole
    REQUIRED INPUT KEYWORDS:
            User must supply a position, proper motion,radial velocity and distance
            (or parallax).    Either scalars or vectors can be supplied.
        (1) Position:
            RA - Right Ascension in *Degrees*
            Dec - Declination in *Degrees*
        (2) Proper Motion
            PMRA = Proper motion in RA in arc units (typically milli-arcseconds/yr)
            PMDEC = Proper motion in Declination (typically mas/yr)
        (3) Radial Velocity
            VRAD = radial velocity in km/s
        (4) Distance or Parallax
            DISTANCE - distance in parsecs
                    or
            PLX - parallax with same distance units as proper motion measurements
                typically milliarcseconds (mas)

    OPTIONAL INPUT KEYWORD:
            /LSR - If this keyword is set, then the output velocities will be
                corrected for the solar motion (U,V,W)_Sun = (-8.5, 13.38, 6.49)
                (Coskunoglu et al. 2011 MNRAS) to the local standard of rest.
                Note that the value of the solar motion through the LSR remains
                poorly determined.
        EXAMPLE:
            (1) Compute the U,V,W coordinates for the halo star HD 6755.
                Use values from Hipparcos catalog, and correct to the LSR
            ra = ten(1,9,42.3)*15.    & dec = ten(61,32,49.5)
            pmra = 627.89  &  pmdec = 77.84         ;mas/yr
            dis = 144    &  vrad = -321.4
            gal_uvw,u,v,w,ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,vrad=vrad,dis=dis,/lsr
                ===>  u=154  v = -493  w = 97        ;km/s

            (2) Use the Hipparcos Input and Output Catalog IDL databases (see
            http://idlastro.gsfc.nasa.gov/ftp/zdbase/) to obtain space velocities
            for all stars within 10 pc with radial velocities > 10 km/s

            dbopen,'hipparcos,hic'      ;Need Hipparcos output and input catalogs
            list = dbfind('plx>100,vrad>10')      ;Plx > 100 mas, Vrad > 10 km/s
            dbext,list,'pmra,pmdec,vrad,ra,dec,plx',pmra,pmdec,vrad,ra,dec,plx
            ra = ra*15.                 ;Need right ascension in degrees
            GAL_UVW,u,v,w,ra=ra,dec=dec,pmra=pmra,pmdec=pmdec,vrad=vrad,plx = plx
            forprint,u,v,w              ;Display results
    METHOD:
            Follows the general outline of Johnson & Soderblom (1987, AJ, 93,864)
            except that U is positive outward toward the Galactic *anti*center, and
            the J2000 transformation matrix to Galactic coordinates is taken from
            the introduction to the Hipparcos catalog.
    REVISION HISTORY:
            Written, W. Landsman                       December   2000
            fix the bug occuring if the input arrays are longer than 32767
            and update the Sun velocity           Sergey Koposov June 2008
            vectorization of the loop -- performance on large arrays
            is now 10 times higher                Sergey Koposov December 2008
    """

    cosd = np.cos(np.deg2rad(dec))
    sind = np.sin(np.deg2rad(dec))
    cosa = np.cos(np.deg2rad(ra))
    sina = np.sin(np.deg2rad(ra))

    # Equivalent of 1 AU/yr in km/s   
    k = 4.74047     
    a_g = np.array([[0.0548755604, +0.4941094279, -0.8676661490], \
                    [0.8734370902, -0.4448296300, -0.1980763734], \
                    [0.4838350155, 0.74698224450, +0.4559837762]])

    # Parallax [mas]
    plx = 1e3 / distance
    
    vec1 = vrad
    vec2 = k * pmra / plx
    vec3 = k * pmdec / plx

    u = (a_g[0,0] * cosa * cosd + a_g[1,0] * sina * cosd + a_g[2,0] * sind) * vec1 + (-a_g[0,0] * sina + a_g[1,0] * cosa) * vec2 + (-a_g[0,0] * cosa * sind - a_g[1,0] * sina * sind + a_g[2,0] * cosd) * vec3
    v = (a_g[0,1] * cosa * cosd + a_g[1,1] * sina * cosd + a_g[2,1] * sind) * vec1 + (-a_g[0,1] * sina + a_g[1,1] * cosa) * vec2 + (-a_g[0,1] * cosa * sind - a_g[1,1] * sina * sind + a_g[2,1] * cosd) * vec3
    w = (a_g[0,2] * cosa * cosd + a_g[1,2] * sina * cosd + a_g[2,2] * sind) * vec1 + (-a_g[0,2] * sina + a_g[1,2] * cosa) * vec2 + (-a_g[0,2] * cosa * sind - a_g[1,2] * sina * sind + a_g[2,2] * cosd) * vec3
    # Change sign convention for U (now positive toward center of Galaxy)
    u = -u

    if (lsr is not None):   
        u = u + lsr[0]
        v = v + lsr[1]
        w = w + lsr[2]

    return (u,v,w)
