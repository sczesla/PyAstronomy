from __future__ import division
import numpy
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import pyasl

class _ZList:
  """
    Calculate the projected distance between the centers of star and planet.
    
    The resulting values will be stored in the `_zlist` attribute. Additionally,
    the indices belonging to in-transit points will be saved in the `_intrans`
    attribute. These can be used as input to calculate the transit light-curve.
    
    Parameters
    ----------
    orbit : string
        Either "circular" or "keplerian".
  """
  
  def _zlistCirc(self, time):
    """ 
      Calculate projected distance in case of circular orbit.
      
      Required parameters: per, i, a (assumes transit center at 0)
      
      Parameters
      ----------
      time : array
          The time points (unit has to be consistent throughout all parameters
          used for calculation).
    """
    # Initialize the z-list (NaN values will be neglected)
    self._zlist = numpy.empty(len(time))
    self._zlist[:] = numpy.NAN
    # w = circular frequency
    w=2.0*numpy.pi/self["per"]
    # The angle needed to realize the orbit inclination
    alpha=(90.0-self["i"])/180.0*numpy.pi
    # Determine the phases (to ensure that the secondary eclipse is not considered)
    phase = time/self["per"]; phase -= numpy.floor(phase)
    # In these cases the 'planet' will be in front of the primary
    self._intrans = numpy.where( numpy.logical_or( phase > 0.75, phase < 0.25 ) )[0]

    # Calculate the 'orbit' and store z-values in __zlist
    self._zlist[self._intrans] = self["a"] * numpy.sqrt( numpy.sin(w*time[self._intrans])**2 + \
                                 numpy.sin(alpha)**2*numpy.cos(w*time[self._intrans])**2 )
  

  def _zlistKep(self, time):
    """ 
      Calculate the projected, normalized distance of star and planet centers.
      
      Required parameters: per, i, a, e, w, T0pa, Omega, p
      Additionally requires: collisionCheck
      
      Note that in this method the line of sight is the -z direction.
      
      Parameters
      ----------
      time : array
          The time points (unit has to be consistent throughout all
          parameters used for calculation).
    """
    # Assign values to KeplerEllipse
    self._ke.a = self["a"]
    self._ke.e = self["e"]
    self._ke.i = self["i"]
    self._ke.Omega = self["Omega"]
    self._ke.per = self["per"]
    self._ke.tau = self["tau"]
    self._ke.w = self["w"]
    
    # Initialize the z-list (NaN values will be neglected)
    self._zlist = numpy.empty(len(time))
    self._zlist[:] = numpy.NAN
    # Get position of body from Keplerian orbit
    pos = self._ke.xyzPos(time)
    
    # Default line of sight is in +z direction (i.e., observer
    # located at -z)
    
    # Calculate projected distance between centers of star and
    # planet sqrt(x**2 + y**2)
    z = numpy.sqrt(pos[::,0]**2 + pos[::,1]**2)
    
    # If enabled, check for body collisions
    if self._collisionCheck:
      r = numpy.sqrt(pos[::,0]**2+pos[::,1]**2+pos[::,2]**2)
      indi = numpy.where(r < 1.+self["p"])[0]
      if len(indi) > 0:
        raise(PE.PyAValError("There is a body collision on this orbit.", \
                             where="_ZList"))
    

    # Determine which values are relevant. In particular, these are those
    # for which the planet is in front of the star (z > 0) and the distance
    # between stellar and planetary distance is lower than 1+p.
    self._intrans = numpy.where(numpy.logical_and(z <= (1.+self["p"]), pos[::,2] < 0.0))[0]
    
    # Calculate the 'orbit' and store z-values in _zlist
    self._zlist[self._intrans] = z[self._intrans]
  
  
  def __init__(self, orbit, cc=True):
    if orbit == "circular":
      self._calcZList = self._zlistCirc
    elif orbit == "keplerian":
      self._collisionCheck = cc
      self._ke = pyasl.KeplerEllipse(1.,1.)
      self._calcZList = self._zlistKep