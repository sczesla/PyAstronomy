Calculate a Keplerian (two body) orbit
========================================

Although the two-body problem has long been solved,
calculation the orbit position of a body in an eccentric
orbit --- maybe a planet --- as a function of time
is not trivial. The major
complication is solving Kepler's Equation. The
classes defined here do this job.

We first give examples demonstrating the
usage of the classes. The API described below.

Example: Invoking the solver for Kepler's Equation
---------------------------------------------------

This example demonstrates how to use the solver
for Kepler's Equation.

::

  from PyAstronomy import pyasl
  
  # Instantiate the solver
  ks = pyasl.MarkleyKESolver()
  
  # Solves Kepler's Equation for a set
  # of mean anomaly and eccentricity.
  # Uses the algorithm presented by
  # Markley 1995.
  M = 0.75
  e = 0.3
  print "Eccentric anomaly: ", ks.getE(M, e)

Example: Calculating the orbit
------------------------------- 

Here we show how the orbit can be calculated.

::

  import numpy
  from PyAstronomy import pyasl
  import matplotlib.pylab as plt
  
  # Instantiate a Keplerian elliptical orbit with
  # semi-major axis of 1.3 length units,
  # period of 2 time units, eccentricity of 0.5, and
  # longitude of ascending node of 70 degrees.
  ke = pyasl.KeplerEllipse(1.3, 2., e=0.5, Omega=70.)
  
  # Get a time axis
  t = numpy.linspace(0, 6.5, 200)
  
  # Calculate the orbit position at the given points
  # in a Cartesian coordinate system.
  pos = ke.xyzPos(t)
  print "Shape of output array: ", pos.shape
  
  # x, y, and z coordinates for 50th time point
  print "x, y, z for 50th point: ", pos[50, ::]
  
  # Calculate orbit radius as a function of the
  radius = ke.radius(t)
  
  # Plot x and y coordinates of the orbit
  plt.subplot(2,1,1)
  plt.plot(pos[::,0], pos[::,1], 'bp')
  # Plot orbit radius as a function of time
  plt.subplot(2,1,2)
  plt.plot(t, radius, 'bp')
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
