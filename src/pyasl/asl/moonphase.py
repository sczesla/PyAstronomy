from __future__ import print_function, division
import numpy as np
from .moonpos import moonpos
from .sunpos import sunpos


def moonphase(jd):
  """
    Computes the illuminated fraction of the Moon at given Julian date(s).
    
    Parameters
    ----------
    jd : float or array
         The Julian date.
         
    Returns
    -------
    Fraction : float or array
        The illuminated fraction [0 - 1] of the Moon.
        Has the same size as `jd`.

    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.

    :IDL - Documentation:  
    
    NAME:
          MPHASE
    PURPOSE:
          Return the illuminated fraction of the Moon at given Julian date(s) 
    
    CALLING SEQUENCE:
          MPHASE, jd, k
    INPUT:
          JD - Julian date, scalar or vector, double precision recommended
    OUTPUT:
          k - illuminated fraction of Moon's disk (0.0 < k < 1.0), same number
              of elements as jd.   k = 0 indicates a new moon, while k = 1 for
              a full moon.
    EXAMPLE:
          Plot the illuminated fraction of the moon for every day in July 
          1996 at 0 TD (~Greenwich noon).
    
          IDL> jdcnv, 1996, 7, 1, 0, jd         ;Get Julian date of July 1
          IDL> mphase, jd+dindgen(31), k        ;Moon phase for all 31 days
          IDL> plot, indgen(31),k               ;Plot phase vs. July day number
  """

  jd = np.array(jd, ndmin=1)
  
  # Earth-Sun distance (1 AU)
  edist = 1.49598e8         
  
  mpos = moonpos(jd)
  ram = mpos[0]*np.pi/180.
  decm = mpos[1]*np.pi/180.
  dism = mpos[2]
  
  spos = sunpos(jd)
  ras = spos[1]*np.pi/180.
  decs = spos[2]*np.pi/180.

  phi = np.arccos( np.sin(decs)*np.sin(decm) + np.cos(decs)*np.cos(decm)*np.cos(ras-ram) )
  inc = np.arctan2( edist * np.sin(phi), dism - edist*np.cos(phi) )
  k = (1 + np.cos(inc))/2.
  
  return np.ravel(k)
