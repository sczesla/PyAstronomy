from __future__ import print_function, division
import numpy
from PyAstronomy.pyaC import pyaErrors as PE

def aitoff(l, b):
  """
    Carry out Aitoff projection. 

    Parameters
    ----------
    l, b : float or array
        The longitude and latitude [deg]
    
    Returns
    -------
    x, y : float or array
        Aitoff-projected coordinates (`x`, `y`)
    
    Notes
    -----

    .. note:: This function was ported from the IDL Astronomy User's Library.
    
    :IDL - Documentation:
  
    pro aitoff,l,b,x,y
    +
     NAME:
           AITOFF
     PURPOSE:
           Convert longitude, latitude to X,Y using an AITOFF projection.
     EXPLANATION:
           This procedure can be used to create an all-sky map in Galactic 
           coordinates with an equal-area Aitoff projection.  Output map 
           coordinates are zero longitude centered.
    
     CALLING SEQUENCE:
           AITOFF, L, B, X, Y 
    
     INPUTS:
           L - longitude - scalar or vector, in degrees
           B - latitude - same number of elements as L, in degrees
    
     OUTPUTS:
           X - X coordinate, same number of elements as L.   X is normalized to
                   be between -180 and 180
           Y - Y coordinate, same number of elements as L.  Y is normalized to
                   be between -90 and 90.
    
     NOTES:
           See AIPS memo No. 46, page 4, for details of the algorithm.  This
           version of AITOFF assumes the projection is centered at b=0 degrees.
    
     REVISION HISTORY:
           Written  W.B. Landsman  STX          December 1989
           Modified for Unix:
                   J. Bloch        LANL SST-9      5/16/91 1.1
           Converted to IDL V5.0   W. Landsman   September 1997
  """
  wasFloat = False
  if isinstance(l, float):
    l = numpy.array([l]); b = numpy.array([b])
    wasFloat = True
  
  sa = l.copy()
  x180 = numpy.where(sa > 180.0)[0]
  if len(x180) > 0: sa[x180] -= 360.
  alpha2 = sa/(2*(180.0/numpy.pi))
  delta = b/(180.0/numpy.pi)   
  r2 = numpy.sqrt(2.)    
  f = 2*r2/numpy.pi   
  cdec = numpy.cos(delta)    
  denom = numpy.sqrt(1. + cdec*numpy.cos(alpha2))
  x = cdec*numpy.sin(alpha2)*2.*r2/denom
  y = numpy.sin(delta)*r2/denom
  x = x*(180.0/numpy.pi)/f
  y = y*(180.0/numpy.pi)/f
  
  if wasFloat:
    return float(x), float(y)
  else:
    return x, y




def inverseAitoff(x, y):
  """
    Carry out an inverse Aitoff projection.
    
    This function reverts to aitoff projection made by the
    function *aitoff*. The result is either two floats or
    arrays (depending on whether float or array was used
    as input) representing longitude and latitude. Both are
    given in degrees with -180 < longitude < +180 and
    -90 < latitude < 90.
    
    Parameters
    ----------
    x : float or array
        A value between -180. and +180.
        (see convention in *aitoff* function).
    y : float or array
        A value between -90. and +90
        (see convention in *aitoff* function).
    
    Returns
    -------
    Deprojected coordinates : float or array
        If arrays are used for input, the function returns an
        array for the longitude, for the latitude, and
        an index array containing those array indices for which
        the reprojection could be carried out.
  """
  wasFloat = False
  if isinstance(x, float):
    x = numpy.array([x]); y = numpy.array([y])
    wasFloat = True
  
  # First, rescale x and y
  x = x/180.0 * 2.0 * numpy.sqrt(2.0)
  y = y/90.0 * numpy.sqrt(2.0)
  
  zsqr = 1.0 - (x/4.0)**2 - (y/2.0)**2
  # Check whether x,y coordinates are within the ellipse of invertible values.
  indi = numpy.where((x**2/8. + y**2/2. - 1.0) <= 0.0)[0]
  if len(indi) == 0:
    raise(PE.PyAValError("Deprojection is not possible.", where="inverseAitoff", why="No values inside valid space."))
  
  z = numpy.sqrt(zsqr[indi])
  l = 2.0 * numpy.arctan( z*x[indi]/(2.0 * (2.0*z**2 - 1.0)) )
  b = numpy.arcsin(z * y[indi])
  
  l = l*180.0 / numpy.pi
  b = b*180.0 / numpy.pi
  
  if wasFloat:
    return float(l), float(b)
  else:
    return l, b, indi
  
