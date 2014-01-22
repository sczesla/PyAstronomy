from PyAstronomy import pyaC
from numpy import arctan2, sin, cos, tan

def positionAngle(ra1, dec1, ra2, dec2, positive=True):
  """
    Compute the position angle.
  
    The position angle is measured from the first position
    from North through East. If the `positive` flag is set
    True (default) the result will be given as an angle from
    0 to 360 degrees. If the flag is set False, the scale is
    -180 to 180 degrees with negative number increasing from
    North through West; note that this is the behavior of
    the posAng IDL routine.
    
    Parameters
    ----------
    ra1 : float
        Right ascension of first object [deg].
    dec1 : float
        Declination of first object [deg].
    ra2 : float
        Right ascension of second object [deg].
    dec2 : float
        Declination of second object [deg].
    positive : boolean, optional
        If True (default), the output will be
        given as an angle between 0 and 360
        degrees. Otherwise, the angle ranges
        from -180 to +180 degrees.
    
    Returns
    -------
    Position angle : float
        The position angle in degrees.
    
  
    Notes
    -----
    
    .. note:: This function was ported from the IDL Astronomy User's Library.
    
    :IDL - Documentation:  
    
     NAME:
           POSANG
     PURPOSE:
           Computes rigorous position angle of source 2 relative to source 1
           
     EXPLANATION:
           Computes the rigorous position angle of source 2 (with given RA, Dec) 
           using source 1 (with given RA, Dec) as the center.
     
     CALLING SEQUENCE:
           POSANG, U, RA1, DC1, RA2, DC2, ANGLE
    
     INPUTS:
           U    -- Describes units of inputs and output:
                   0:  everything radians
                   1:  RAx in decimal hours, DCx in decimal
                           degrees, ANGLE in degrees
           RA1  -- Right ascension of point 1
           DC1  -- Declination of point 1
           RA2  -- Right ascension of point 2
           DC2  -- Declination of point 2
    
       OUTPUTS:
           ANGLE-- Angle of the great circle containing [ra2, dc2] from
                   the meridian containing [ra1, dc1], in the sense north
                   through east rotating about [ra1, dc1].  See U above 
                   for units.
    
       PROCEDURE:
           The "four-parts formula" from spherical trig (p. 12 of Smart's
           Spherical Astronomy or p. 12 of Green' Spherical Astronomy).
    
       EXAMPLE:
           For the star 56 Per, the Hipparcos catalog gives a position of 
           RA = 66.15593384, Dec = 33.94988843 for component A, and 
           RA = 66.15646079, Dec =  33.96100069 for component B.   What is the
           position angle of B relative to A?
    
           IDL> RA1 = 66.15593384/15.d   & DC1 = 33.95988843
           IDL> RA2 = 66.15646079/15.d   & DC2 = 33.96100069
           IDL> posang,1,ra1,dc1,ra2,dc2, ang
                will give the answer of ang = 21.4 degrees
       NOTES:
           (1) If RA1,DC1 are scalars, and RA2,DC2 are vectors, then ANGLE is a
           vector giving the position angle between each element of RA2,DC2 and 
           RA1,DC1.   Similarly, if RA1,DC1 are vectors, and RA2, DC2 are scalars,
           then DIS is a vector giving the position angle of each element of RA1, 
           DC1 and RA2, DC2.    If both RA1,DC1 and RA2,DC2 are vectors then ANGLE 
           is a vector giving the position angle between each element of RA1,DC1 
           and the corresponding element of RA2,DC2.    If then vectors are not the
           same length, then excess elements of the longer one will be ignored.
    
           (2) Note that POSANG is not commutative -- the position angle between
            A and B is theta, then the position angle between B and A is 180+theta 
       PROCEDURE CALLS:
            ISARRAY()
       HISTORY:
           Modified from GCIRC, R. S. Hill, RSTX, 1 Apr. 1998
           Use V6.0 notation W.L. Mar 2011
  """

  # Convert into rad               
  rarad1 = pyaC.degtorad(ra1)
  rarad2 = pyaC.degtorad(ra2)
  dcrad1 = pyaC.degtorad(dec1)
  dcrad2 = pyaC.degtorad(dec2)

  radif  = rarad2 - rarad1

  angle  = arctan2(sin(radif), cos(dcrad1)*tan(dcrad2)-sin(dcrad1)*cos(radif))

  result = pyaC.radtodeg(angle)
  
  if positive and (result < 0.0):
    result+= 360.0

  return result  
