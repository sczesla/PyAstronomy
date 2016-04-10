.. _keplerorbitpyasl:

Calculate a Keplerian (two body) orbit
========================================

.. p23ready

Although the two-body problem has long been solved,
calculation the orbit position of a body in an eccentric
orbit --- maybe a planet --- as a function of time
is not trivial. The major
complication is solving Kepler's Equation. The
classes defined here do this job.

The definitions and most of the formulae used in this class
derive from the book "Orbital Motion" by A.E. Roy.

Orbital elements and orientation of the orbit
----------------------------------------------------

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


Example: Invoking the solver for Kepler's Equation
---------------------------------------------------

This example demonstrates how to use the solver
for Kepler's Equation.

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Instantiate the solver
    ks = pyasl.MarkleyKESolver()
    
    # Solves Kepler's Equation for a set
    # of mean anomaly and eccentricity.
    # Uses the algorithm presented by
    # Markley 1995.
    M = 0.75
    e = 0.3
    print("Eccentric anomaly: ", ks.getE(M, e))


Example: Calculating the orbit
------------------------------- 

Here we show how the orbit can be calculated.

::
    
    from __future__ import print_function, division
    import numpy as np
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    
    # Instantiate a Keplerian elliptical orbit with
    # semi-major axis of 1.3 length units,
    # a period of 2 time units, eccentricity of 0.5,
    # longitude of ascending node of 70 degrees, an inclination
    # of 10 deg, and a periapsis argument of 110 deg.
    ke = pyasl.KeplerEllipse(1.3, 2., e=0.5, Omega=70., i=10.0, w=110.0)
    
    # Get a time axis
    t = np.linspace(0, 1.9, 200)
    
    # Calculate the orbit position at the given points
    # in a Cartesian coordinate system.
    pos = ke.xyzPos(t)
    print("Shape of output array: ", pos.shape)
    
    # x, y, and z coordinates for 50th time point
    print("x, y, z for 50th point: ", pos[50, ::])
    
    # Calculate orbit radius as a function of the
    radius = ke.radius(t)
    
    # Calculate velocity on orbit
    vel = ke.xyzVel(t)
    
    # Find the nodes of the orbit (Observer at -z)
    ascn, descn = ke.xyzNodes_LOSZ()
    
    # Plot x and y coordinates of the orbit
    plt.subplot(2,1,1)
    plt.title("Periapsis (red diamond), Asc. node (green circle), desc. node (red circle)")
    plt.xlabel("East ->")
    plt.ylabel("North ->")
    plt.plot([0], [0], 'k+', markersize=9)
    plt.plot(pos[::,1], pos[::,0], 'bp')
    # Point of periapsis
    plt.plot([pos[0,1]], [pos[0,0]], 'rd')
    # Nodes of the orbit
    plt.plot([ascn[1]], [ascn[0]], 'go', markersize=10)
    plt.plot([descn[1]], [descn[0]], 'ro', markersize=10)
    # Plot RV
    plt.subplot(2,1,2)
    plt.xlabel("Time")
    plt.ylabel("Radial velocity [length/time]")
    plt.plot(t, vel[::,2], 'r.-')
    plt.show()


Module API
---------------

.. currentModule:: PyAstronomy.pyasl

The module defines the following classes:

  - :py:class:`KeplerEllipse`
  - :py:class:`MarkleyKESolver`

The `KeplerEllipse` class calculates the orbit and provides
some convenience functions. For instance, the foci of the ellipse,
and the peri- and apastron positions can be calculated.

The `MarkleyKESolver` class implements a solver for Kepler's
equation, which is needed to calculate the orbit as a function
of time.
