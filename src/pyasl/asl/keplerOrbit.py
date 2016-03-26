from __future__ import print_function, division
import numpy
from numpy import pi, abs, sqrt, cos, sin, arccos, arctan, tan
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo

class MarkleyKESolver:
  """
    Solve Kepler's Equation.
    
    This class implements a solver for Kepler's Equation:
    .. math::
        M = E - sin(E),
    where M is the "mean anomaly" and E is the "eccentric anomaly".
    The implementation follows the prescription given by
    Markley (Markley 1995, CeMDA, 63, 101).
  """
  
  pi2 = pi**2

  def _alpha(self, e, M):
    """
      Solve Eq. 20
    """
    return ( 3.*self.pi2 + 1.6*pi*(pi-abs(M))/(1.+e) )/(self.pi2 - 6.)
  
  def _d(self, alpha, e):
    """
      Solve Eq. 5
    """
    return 3.*(1. - e) + alpha*e
  
  def _r(self, alpha, d, M, e):
    """
      Solve Eq. 10
    """
    return 3.*alpha*d * (d-1.+e)*M + M**3
  
  def _q(self, alpha, d, e, M):
    """
      Solve Eq. 9
    """
    return 2.*alpha*d*(1.-e) - M**2
  
  def _w(self, r, q):
    """
      Solve Eq. 14
    """
    return (abs(r) + sqrt(q**3 + r**2))**(2./3.)
  
  def _E1(self, d, r, w, q, M):
    """
      Solve Eq. 15
    """
    return ( 2.*r*w/(w**2 + w*q + q**2) + M ) / d
  
  def _f01234(self, e, E, M):
    """
      Solve Eq. 21, 25, 26, 27, and 28 (f, f', f'', f''', and f'''')
    """
    f0 = E - e*sin(E) - M
    f1 = 1. - e*cos(E)
    f2 = e*sin(E)
    return f0, f1, f2, 1.-f1, -f2
  
  def _d3(self, E, f):
    """
      Solve Eq. 22 
    """
    return -f[0]/(f[1] - 0.5*f[0]*f[2]/f[1])
  
  def _d4(self, E, f, d3):
    """
      Solve Eq. 23
    """
    return -f[0]/( f[1] + 0.5*d3*f[2] + (d3**2)*f[3]/6. )
  
  def _d5(self, E, f, d4):
    """
      Solve Eq. 24
    """
    return -f[0]/( f[1] + 0.5*d4*f[2] + d4**2*f[3]/6. + d4**3*f[4]/24.)
  
  def precisionTest(self, trials=10000):
    """
      Carry out a test of the achieved precision.
      
      Generate random numbers for the mean anomaly
      and the eccentricity and calculate the eccentric
      anomaly. Use the number in Kepler's Equation
      and compare the resulting mean anomaly with
      the input.
      
      Parameters
      ----------
      trials : int, optional
          The number of trial evaluations.
      
      Returns
      -------
      Deviation : float
          The largest determined deviation between
          input and output mean anomaly.
    """
    maxdev = 0.0
    for i in smo.range(trials):
      M = numpy.random.random() * 100. - 50.
      e = numpy.random.random()
      E = self.getE(M, e)
      M1 = E - e*sin(E)
      M = M - numpy.floor(M/(2.*pi))*2*pi
      maxdev = max(maxdev, abs(M-M1))
    print("Maximal deviation in resulting mean anomaly (M): ", maxdev)
    return maxdev
  
  def getE(self, M, e):
    """
      Solve Kepler's Equation for the "eccentric anomaly", E.
      
      Parameters
      ----------
      M : float
          Mean anomaly.
      e : float
          Eccentricity
      
      Returns
      -------
      Eccentric anomaly: float
          The solution of Kepler's Equation
    """
    # For the mean anomaly, use values between
    # -pi and pi.
    flip = False
    self.M = M - (numpy.floor(M/(2.*pi)) * 2.*pi)
    if self.M > pi:
      self.M = 2.*pi - self.M
      # Flip the sign of result
      # if this happened
      flip = True
    self.e = e
    if M == 0.0: return 0.0
    self.alpha = self._alpha(e, self.M)
    self.d = self._d(self.alpha, e)
    self.r = self._r(self.alpha, self.d, self.M, e)
    self.q = self._q(self.alpha, self.d, e, self.M)
    self.w = self._w(self.r, self.q)
    self.E1 = self._E1(self.d, self.r, self.w, self.q, self.M)
    self.f = self._f01234(e, self.E1, self.M)
    self.d3 = self._d3(self.E1, self.f)
    self.d4 = self._d4(self.E1, self.f, self.d3)
    self.d5 = self._d5(self.E1, self.f, self.d4)
    # Eq. 29
    self.E5 = self.E1 + self.d5
    if flip:
      self.E5 = 2.*pi - self.E5
    self.E = self.E5
    return self.E5


class KeplerEllipse(object):
  """
    Calculate a Kepler orbit.
    
    The definitions and most of the formulae used in this class
    derive from the book "Orbital Motion" by A.E. Roy.
    
    :Orientation of the ellipse in the coordinate system:
        For zero inclination the ellipse is located in the x-y plane.
        If the eccentricity is increased, the periastron will lie
        in +x direction. If the inclination is increased, the ellipse
        will be rotating around the x-axis, so that +y is rotated
        toward +z. An increase in Omega corresponds to a rotation
        around the z-axis so that +x is rotated toward +y.
        Changing `w`, i.e., the argument of the periastron, will
        not change the plane of the orbit, but rather represent a
        rotation of the orbit in the plane. In particular, the
        periapsis is shifted in the direction of motion.
        
    :Orbital angular momentum:
        For all parameters but semi-major axis and orbital period set to zero,
        the (orbital) angular momentum points into the +z direction. For an
        inclination of 90 deg (the remaining parameters remaining zero),
        it points in the -y direction.
    
    :Orientation of the ellipse in the sky:
        To project the ellipse onto the sky, the coordinate system
        should be oriented so that the +x direction points North and
        the +y direction points East (direction of increasing right
        ascension). The +z axis must be chosen so that the coordinate
        system becomes right handed. If the line of sight (LOS) points
        in the +z direction, i.e., the observer is located on the
        negative z axis, the parameters assume their conventional
        meaning.
    
    :The ascending and descending nodes:
        For systems outside the Solar System, the ascending node is the
        point where the body "crosses" the plane of the sky away from the
        observer. Likewise, the descending node is the point where the
        plane is crossed with the body approaching the observer. For the
        coordinate system described above and a value of zero for the longitude
        of the ascending node, the latter is in the North and rotates
        toward East (i.e., +y) when the longitude of the ascending node
        is increased.  
        
    :The argument and longitude of periapsis:
        The argument of periapsis is the angle between the ascending node
        and the periapsis of the body measured in the direction of motion.
        For exoplanets with circular orbits, for which no well-defined periapsis
        exists, the argument of periapsis is often chosen so that time
        of periapsis and central transit time coincide. For the planet, this
        is the case if the argument of periapsis is -90 deg. However, in the exoplanet
        literature, the argument of periapsis often refers to the *stellar* orbit
        (see, e.g., Pollacco et al. 2008, MNRAS 385, 1576-1584, Sect. 3.2.1). In
        this case, the corresponding value is +90 deg.
        The so-called longitude of the periapsis is given by the sum of the
        longitude of the ascending node and the argument of periapsis.
        
    Parameters
    ----------
    a : float
        Semi-major axis
    per : float
        Orbital period
    e : float, optional
        Orbital eccentricity (0-1).
    tau : float, optional
        Time of periapsis passage.
    Omega : float, optional
        Longitude of the ascending node [deg].
    w : float, optional
        Argument of periapsis [deg]. Note that the longitude
        if periapsis is given by Omega+w.
    i : float, optional
        Orbit inclination [deg].
    ks : Class object, optional
        The solver for Kepler's Equation. Default is the
        `MarkleyKESolver`. Each solver must have a `getE`
        method, which takes either a float or array of float
        of mean anomalies and the eccentricity, and returns
        the associated eccentric anomalies.
    
    Attributes
    ----------
    i : float
        Orbit inclination [deg].
    w : float
        Argument of periapsis [deg]
    Omega : float
        Longitude of the ascending node [deg]
    e : float
        Eccentricity
    a : float
        Semi-major axis
    per : float
        Orbital period
    tau : float
        Time of periapsis passage
    ks : Class object
        Solver for Kepler's equation
    _n : float
        Circular frequency.
  """
  
  def _getEccentricAnomaly(self, t):
    """
      Calculate eccentric anomaly.
      
      Parameters
      ----------
      t : array of float
          The times at which to calculate the eccentric anomaly, E.
      
      Returns
      -------
      E : Array of float
          The eccentric anomaly.
    """
    M = self.meanAnomaly(t)
    if not hasattr(t, "__iter__"):
      return self.ks.getE(M, self.e)
    E = numpy.zeros(len(t))
    for i in smo.range(len(t)):
      E[i] = self.ks.getE(M[i], self.e)
    return E
       
  def meanAnomaly(self, t):
    """
      Calculate the mean anomaly.
      
      Parameters
      ----------
      t : float or array
          The time axis.
      
      Returns
      -------
      Mean anomaly : float or array
          The mean anomaly (whether float or array
          depends on input).
    """
    return self._n*(t-self.tau)
  
  def radius(self, t, E=None):
    """
      Calculate the orbit radius.
      
      Parameters
      ----------
      t : float or array
          The time axis.
      E : float or array, optional
          If known, the eccentric anomaly corresponding
          to the time points. If not given, the numbers
          will be calculated.
      
      Returns
      -------
      Radius : float or array
          The orbit radius at the given points in time.
          Type depends on input type.
    """
    if E is None:
      E = self._getEccentricAnomaly(t)
    return self.a * (1. - self.e*cos(E))
  
  def xyzPos(self, t, getTA=False):
    """
      Calculate orbit position.
      
      Parameters
      ----------
      t : float or array
          The time axis.
      getTA : boolean, optional
          If True, returns the "true anomaly" as a function
          of time (default is False).
      
      Returns
      -------
      Position : array
          The x, y, and z coordinates of the body at the given time.
          If the input was an array, the output will be an array of
          shape (input-length, 3), holding the positions at the given
          times.
      True anomaly : float or array
          Is returned only if `getTA` is set to True. The true anomaly
          at the specified times.
    """
    E = self._getEccentricAnomaly(t)
    r = self.radius(t, E=E)
    f = arctan( sqrt((1.+self.e)/(1.-self.e)) * tan(E/2.) ) * 2.0
    wf = self._w + f
    cos_Omega = cos(self._Omega)
    sin_Omega = sin(self._Omega)
    cos_i = cos(self._i)
    sin_i = sin(self._i)
    if not hasattr(wf, "__iter__"):
      cos_wf = cos(wf)
      sin_wf = sin(wf)
      xyz = numpy.array([cos_Omega*cos_wf - sin_Omega*sin_wf*cos_i,
                         sin_Omega*cos_wf + cos_Omega*sin_wf*cos_i,
                         sin_wf*sin_i
                         ]) * r
    else:
      # Assume it is an array
      xyz = numpy.zeros( (len(t), 3) )
      for i in smo.range(len(t)):
        cos_wf = cos(wf[i])
        sin_wf = sin(wf[i])
        xyz[i,::] = numpy.array([cos_Omega*cos_wf - sin_Omega*sin_wf*cos_i,
                         sin_Omega*cos_wf + cos_Omega*sin_wf*cos_i,
                         sin_wf*sin_i
                         ]) * r[i]
    if not getTA:
      return xyz
    else:
      return xyz, f
  
  def xyzVel(self, t):
    """
      Calculate orbit velocity.
      
      Parameters
      ----------
      t : float or array
          The time axis.
      
      Returns
      -------
      Velocity : array
          The x, y, and z components of the body's velocity at the
          given time. If the input was an array, the output will be
          an array of shape (input-length, 3), holding the velocity
          components at the given times. The unit is that of the
          semi-major axis divided by that of the period.
    """
    # From AE ROY "Orbital motion" p. 102
    cos_Omega = cos(self._Omega)
    sin_Omega = sin(self._Omega)
    cos_i = cos(self._i)
    sin_i = sin(self._i)
    cos_w = cos(self._w)
    sin_w = sin(self._w)

    E = self._getEccentricAnomaly(t)
    l1 = cos_Omega*cos_w - sin_Omega*sin_w*cos_i
    l2 = -cos_Omega*sin_w - sin_Omega*cos_w*cos_i
    m1 = sin_Omega*cos_w + cos_Omega*sin_w*cos_i
    m2 = -sin_Omega*sin_w + cos_Omega*cos_w*cos_i
    n1 = sin_w*sin_i
    n2 = cos_w*sin_i
    b = self.a * sqrt(1. - self.e**2)
    r = self.radius(t, E)
    nar = self._n * self.a / r
    if not hasattr(t, "__iter__"):
      bcos_E = b*cos(E)
      asin_E = self.a*sin(E)
      vel = nar * numpy.array([l2*bcos_E - l1*asin_E,
                               m2*bcos_E - m1*asin_E,
                               n2*bcos_E - n1*asin_E])
    else:
      # Assume it is an array
      vel = numpy.zeros( (len(t), 3) )
      for i in smo.range(len(t)):
        bcos_E = b*cos(E[i])
        asin_E = self.a*sin(E[i])
        vel[i,::] = nar[i] * numpy.array([l2*bcos_E - l1*asin_E,
                                          m2*bcos_E - m1*asin_E,
                                          n2*bcos_E - n1*asin_E])
    return vel
  
  def xyzPeriastron(self):
    """
      The position of the periastron.
      
      Returns
      -------
      Periastron : array of float
          The x, y, and z coordinates of the periastron. 
    """
    return self.xyzPos(self.tau)
  
  def xyzApastron(self):
    """
      The position of the apastron.
      
      The apastron is the point of greatest distance.
      
      Returns
      -------
      Apastron : array of float
          The x, y, and z coordinates of the apastron
    """
    return self.xyzPos(self.tau + 0.5*self.per)

  def xyzCenter(self):
    """
      Center of the ellipse
      
      Returns
      -------
      Center : array of float
          x, y, and z coordinates of the center of the Ellipse.
    """
    return (self.xyzPeriastron() + self.xyzApastron())/2.0

  def xyzFoci(self):
    """
      Calculate the foci of the ellipse
      
      Returns
      -------
      Foci : Tuple of array of float
          A tuple containing two arrays, which hold the x, y, and z
          coordinates of the foci.
    """
    peri = self.xyzPeriastron()
    apas = self.xyzApastron()
    center = (peri + apas) / 2.0
    direc = (peri - apas) / numpy.sqrt( ((peri - apas)**2).sum() )
    ae = self.a * self.e
    return (center + ae*direc, center - ae*direc)
    
  def xyzNodes_LOSZ(self, los="+z"):
    """
      Calculate the nodes of the orbit for LOS in +/-z direction.
      
      The nodes of the orbit are the points at which
      the orbit cuts the plane of the sky. In this case,
      these are the points at which the z-coordinate
      vanishes, i.e., the x-y plane is regarded the plane
      of the sky.
      
      Parameters
      ----------
      los : string, {+z,-z}, optional
          Line of sight points either in +z direction (observer
          at -z) or vice versa. Changing the direction
          interchanges the ascending and descending node.
      
      Returns
      -------
      Nodes : Tuple of two coordinate arrays
          Returns the xyz coordinates of both nodes. The first is the
          ascending node and the second is the descending node.
    """
    # f = -w (z-component vanishes there)
    E = arctan( tan(-self._w/2.0) * sqrt((1.-self.e)/(1.+self.e)) ) * 2.
    M = E - self.e * sin(E)
    t = M/self._n + self.tau
    node1 = self.xyzPos(t)
    # Velocity is used to distinguish nodes
    v1 = self.xyzVel(t)
    # f = -w + pi
    E = arctan( tan((-self._w + pi)/2.0) * sqrt((1.-self.e)/(1.+self.e)) ) * 2.
    M = E - self.e * sin(E)
    t = M/self._n + self.tau
    node2 = self.xyzPos(t)
    # Find the ascending and descending node
    from PyAstronomy.pyasl import LineOfSight
    l = LineOfSight(los).los
    if abs(l[2]) != 1.0:
      raise(PE.PyAValError("Invalid line of sight (LOS): " + str(los), \
                           where="xyzNodes_LOSZ", \
                           solution="Use '-z' or '+z'."))

    if l[2] == 1.0:
      # Looking in +z direction
      if v1[2] > 0.0:
        # First node is ascending
        return (node1, node2)
      else:
        return (node2, node1)
    else:
      # Looking in -z direction
      if v1[2] < 0.0:
        # First node is ascending
        return (node1, node2)
      else:
        return (node2, node1)
    
  def xyzNodes(self):
    """
      Calculate the nodes of the orbit.
      
      The nodes of the orbit are the points at which
      the orbit cuts the observing plane. In this case,
      these are the points at which the z-coordinate
      vanishes, i.e., the x-y plane is regarded the plane
      of observation.
      
      Returns
      -------
      Nodes : Tuple of two coordinate arrays
          Returns the xyz coordinates of both nodes. 
    """
    raise(PE.PyADeprecationError("xyzNodes is deprecated.", \
                                 solution="Please use 'xyzNodes_LOSZ' instead."))

  def orbAngMomentum(self, t=0.0):
    """
      The specific orbital angular momentum of a body in the current orbit.
      
      Parameters
      ----------
      t : float, optional
          The time used to calculate the angular momentum. As
          is it constant for a Kepler orbit, this parameter is
          of no relevance.
      
      Results
      -------
      Specific angular momentum : array with three elements
          The specific angular momentum (i.e., per unit mass)
          of a body on the current orbit.
    """
    r = self.xyzPos(t)
    v = self.xyzVel(t)
    return numpy.cross(r, v)
  
  def projPlaStDist(self, t):
    """
      Calculate the sky-projected planet-star separation.
      
      Parameters
      ----------
      t : float or array
          The time axis.
      
      Returns
      -------
      Position : array
          The sky-projected planet-star separation at the given time.
          If the input was an array, the output will be an array, 
          holding the separations at the given times.
    """
    p = self.a * (1.- self.e**2)
    E = self._getEccentricAnomaly(t)
    f = arctan( sqrt((1.+self.e)/(1.-self.e)) * tan(E/2.) ) * 2.0
    wf = self._w + f
    if not hasattr(wf, "__iter__"):
      psdist = p/(1.+self.e*cos(f))*sqrt(1.-sin(self._i)**2*sin(wf)**2)

    else:
      # Assume it is an array
      psdist = numpy.zeros(len(t))
      for i in smo.range(len(t)):
        psdist[i] = p/(1.+self.e*cos(f[i]))*sqrt(1.-sin(self._i)**2*sin(wf[i])**2)

    return psdist
  
  def yzCrossingTime(self):
    """
      Calculate times of crossing the yz-plane.
      
      This method calculates the times at which
      the yz-plane is crossed by the orbit. This
      is equivalent to finding the times where
      x=0.
      
      Returns
      -------
      Time 1 : float
          First crossing time defined as having POSITIVE
          y position.
      Time 2 : float
          Second crossing time defined as having NEGATIVE
          y position.
    """
    if abs(self._Omega) < 1e-16:
      f = -self._w + pi/2.0
    else:
      f = -self._w + arctan(1.0/(tan(self._Omega)*cos(self._i)))
    E = 2.*arctan(sqrt((1-self.e)/(1.+self.e)) * tan(f/2.))
    t1 = (E - self.e*sin(E))/self._n + self.tau
    p1 = self.xyzPos(t1)
    f += pi
    E = 2.*arctan(sqrt((1-self.e)/(1.+self.e)) * tan(f/2.))
    t2 = (E - self.e*sin(E))/self._n + self.tau
    
    t1 -= self._per * numpy.floor(t1/self._per)
    t2 -= self._per * numpy.floor(t2/self._per)
    
    if p1[1] >= 0.0:
      # y position of p1 is > 0
      return (t1, t2)
    else:
      return (t2, t1)
 
  def xzCrossingTime(self):
    """
      Calculate times of crossing the xz-plane.
      
      This method calculates the times at which
      the xz-plane is crossed by the orbit. This
      is equivalent to finding the times where
      y=0.
      
      Returns
      -------
      Time 1 : float
          First crossing time defined as having POSITIVE
          x position.
      Time 2 : float
          Second crossing time defined as having NEGATIVE
          x position.
    """
    f = -self._w - arctan(tan(self._Omega)/cos(self._i))
    E = 2.*arctan(sqrt((1-self.e)/(1.+self.e)) * tan(f/2.))
    t1 = (E - self.e*sin(E))/self._n + self.tau
    p1 = self.xyzPos(t1)
    f += pi
    E = 2.*arctan(sqrt((1-self.e)/(1.+self.e)) * tan(f/2.))
    t2 = (E - self.e*sin(E))/self._n + self.tau
    
    t1 -= self._per * numpy.floor(t1/self._per)
    t2 -= self._per * numpy.floor(t2/self._per)
    
    if p1[0] >= 0.0:
      # y position of p1 is > 0
      return (t1, t2)
    else:
      return (t2, t1)
        
  def _setPer(self, per):
    self._per = per
    self._n = 2.0*pi/self._per
  
  def _seti(self, i):
    self._i = i/180.*pi
  
  def _setw(self, w):
    self._w = w/180.*pi
  
  def _setOmega(self, omega):
    self._Omega = omega/180.*pi
  
  per = property(lambda self: self._per, _setPer, doc="The orbital period.")
  i = property(lambda self: self._i/pi*180, _seti)
  w = property(lambda self: self._w/pi*180, _setw)
  Omega = property(lambda self: self._Omega/pi*180, _setOmega)
  
  def __init__(self, a, per, e=0, tau=0, Omega=0, w=0, i=0, ks=MarkleyKESolver):
    # i, w, Omega are properties so that the numbers can be given in
    # deg always. The underscored attributes are in rad.
    self.i = i
    self.w = w
    self.Omega = Omega
    self.e = e
    self.a = a
    self.per = per
    self.tau = tau
    self.ks = ks()



def phaseAngle(pos, los='-z'):
  """
    Calculate the phase angle.
    
    The phase angle is the angle between the star and the Earth (or Sun)
    as seen from the planet (e.g., Seager et al. 1999, ApJ, 504).
    The range of the phase angle is 0 - 180 degrees. In the calculations,
    it is assumed to be located at the center of the coordinate system.
    
    Parameters
    ----------
    pos : array
        Either a one-dimensional array with xyz coordinate or a
        [3,N] array containing N xyz positions.
    los : LineOfSight object
        A `LineOfSight` object from the pyasl giving the line of
        sight.
    
    Returns
    -------
    Phase angle : The phase angle in degrees. Depending on the input,
        it returns a single value or an array. 
  """
  from PyAstronomy.pyasl import LineOfSight
  l = LineOfSight(los)
  if pos.shape == (3,):
    # It is a single value
    return numpy.arccos( numpy.sum(-pos * (-l.los)) / \
            numpy.sqrt(numpy.sum(pos**2)))/numpy.pi*180. 
  else:
    # It is an array of positions
    N = len(pos[::,0])
    result = numpy.zeros(N)
    for i in smo.range(N):
      print(i, numpy.sum((-pos[i,::]) * (-l.los)), pos[i,::])
      result[i] = numpy.arccos( numpy.sum((-pos[i,::]) * (-l.los)) / \
                  numpy.sqrt(numpy.sum(pos[i,::]**2)) )
    return result/numpy.pi*180.
