The n-dimensional fitting tutorial
====================================

.. currentmodule:: PyAstronomy.funcFit

The funcFit package supports n-dimensional fitting, which means
that both the domain and the range may be multidimensional.

An example for a multidimensional range would be the fitting of
a circle or orbit, when time (1d) is mapped to image coordinates (2d).
Fitting, for example, the structure of a point spread function (PSF)
would be an example for a multidimensional domain, viz., the image
coordinates, which are mapped onto a one-dimensional range (the flux
or intensity). 

Fitting a circular orbit
---------------------------

This example demonstrates how to fit a 2d model (location in plane) depending
on a single variable (time), so that there is a mapping of the form

.. math:: f: (R \times R^m) \rightarrow R^n ,

where the :math:`R^m` denotes the parameter vector.

In particular, we assume
x,y-position measurements at a number of times.
Furthermore, all x and y measurements have an error.

::

  import numpy as np
  import matplotlib.pylab as plt
  from PyAstronomy import funcFit as fuf
  
  # Get the circular model and assign
  # parameter values
  c = fuf.Circle2d()
  c["r"] = 1.0
  c["t0"] = 0.0
  c["per"] = 3.0
  
  # Evaluate the model at a number of
  # time stamps
  t = np.linspace(0.0, 10.0, 20)
  pos = c.evaluate(t)
  
  # Add some error to the "measurement"
  pos += np.reshape(np.random.normal(0.0, 0.2, pos.size), pos.shape)
  err = np.reshape(np.ones(pos.size), pos.shape) * 0.2
  
  # Define free parameters and fit the model
  c.thaw(["r", "t0", "per"])
  c.fit(t, pos, yerr=err)
  c.parameterSummary()
  
  # Evaluate the model at a larger number of
  # points for plotting
  tt = np.linspace(0.0, 10.0, 200)
  model = c.evaluate(tt)
  
  # Plot the result
  plt.errorbar(pos[::,0], pos[::,1], yerr=err[::,1], \
               xerr=err[::,0], fmt='bp')
  plt.plot(model[::,0], model[::,1], 'r--')
  plt.show()

.. note:: The only difference between the "normal" one-dimensional case and that
          exemplified here, is the dimension of the `y` and error arrays.

.. _creatingCoordinateArrays:

Creating "coordinate arrays"
-----------------------------

If one has multiple coordinate axes, such as in the case of an image,
it is often convenient to represent the coordinates using an appropriate
array. `funcFit` provides the `coordinateGrid` function to help
constructing an appropriate array. In particular, for n coordinate axes,
the function will construct an array with n+1 dimensions, which
gives the physical coordinates for every array index.

.. note:: The user is free to use any format to give the coordinates. For
          example, the `meshgrid` in numpy could be applied to obtain a
          similar result. The :ref:`gauss2dExample` example demonstrates
          both use cases.

::

    from __future__ import print_function, division
    from PyAstronomy import funcFit as fuf
    import numpy as np
    
    # Constructing the two individual coordinate axes
    x = np.linspace(-2.,2.,50)
    y = np.linspace(-2.,2.,50)
    
    # Applying funcFit's "coordinateGrid" helper function
    # to built appropriate array-index -> coordinate mapping
    # needed for nD fitting.
    g = fuf.coordinateGrid(x, y)
    
    print("(x, y) coordinates at index (11, 28): ", g[11,28])



.. _gauss2dExample:

Fitting image data with a 2d-Gaussian model
---------------------------------------------

The following two examples demonstrate how to fit 2d image data
with a 2d Gaussian model using funcFit, i.e., we have the mapping

.. math:: f: (R^2 \times R^m) \rightarrow R ,

where, again, :math:`R^m` denotes the parameter vector, :math:`R^2` is the image
coordinate, and the result is the "level", "flux", or whatever the image shows.

In this case, we need to take care of the parameter representation.
In principle, we may choose the coordinate format ad libitum, however,
the `evaluate` method of the fitting objects has to understand the format.
Below, we show two conceivable approaches to specify the parameters.

Let us say, we are given both an x- and a y-axis
specifying the coordinates in the image. As demonstrated in :ref:`creatingCoordinateArrays`,
the function :py:func:`coordinateGrid` can be used to
construct an appropriate coordinate array mapping array index to
physical coordinates. The class :py:class:`GaussFit2d` expects such a
coordinate array, and the example below demonstrates how to use it.

.. note:: Below the example, the `evaluate` method of `GaussFit2d` is
          shown to illustrate the use of the coordinate array.

::

  from PyAstronomy import funcFit as fuf
  import numpy as np
  import matplotlib.pylab as plt
  
  # Constructing the individual coordinate axes
  x = np.linspace(-2.,2.,50)
  y = np.linspace(-2.,2.,50)
  # Applying funcFit's "coordinateGrid" helper function
  # to built appropriate array-index -> coordinate mapping
  # needed for nD fitting.
  g = fuf.coordinateGrid(x, y)
  
  # Create the 2d-Gaussian model and assign
  # some model parameters.
  gf = fuf.GaussFit2d()
  gf["sigx"] = 0.75
  gf["sigy"] = 0.4
  gf["A"] = 1.0
  gf["rho"] = 0.4
  
  # Get the "data" by evaluating the model
  # and adding some noise. Note that the coordinate
  # mapping (array g) is passed to evaluate here.
  im = gf.evaluate(g)
  im += np.reshape(np.random.normal(0.0, 0.1, 2500), (50,50))
  err = np.ones((50,50))*0.1
  
  # Thaw parameters and fit
  gf.thaw(["A", "rho"])
  gf.fit(g, im, yerr=err)
  
  # Show the resulting parameter values ...
  gf.parameterSummary()
  
  # ... and plot the result.
  plt.title("Image data")
  plt.imshow(np.transpose(im), origin="lower")
  plt.show()
  plt.title("Residuals")
  plt.imshow(np.transpose(im - gf.evaluate(g)), origin="lower")
  plt.show()


The `evaluate` method of the `GaussFit2d` class has the following from:
  
::

  # The "evaluate" method of class GaussFit2d
  # Note how the coordinate grid (co argument) is used in the evaluation
  # of the model. Basically co[::,::,0] gives the x-coordinate for every
  # image pixel. Likewise, co[::,::,1] defines the y-coordinate.

  def evaluate(self, co):
    """
      Evaluates the model for current parameter values.
      
      Parameters
      ----------
      co : array
           Specifies the points at which to evaluate the model.
    """
    if (self["sigx"] <= 0.0) or (self["sigy"] <= 0.0):
      raise(PE.PyAValError("Width(s) of Gaussian must be larger than zero.", \
                           solution="Change width ('sigx/y')."))
    if self["rho"] > 1.0:
      raise(PE.PyAValError("The correlation coefficient must be 0 <= rho <= 1.", \
                           solution="Change width ('sigx/y')."))
    result = 1.0/(2.*pi*self["sigx"]*self["sigy"]*sqrt(1.-self["rho"]**2)) * \
        exp( ((co[::,::,0]-self["mux"])**2/self["sigx"]**2 + (co[::,::,1]-self["muy"])**2/self["sigy"]**2 - \
            2.*self["rho"]*(co[::,::,0]-self["mux"])*(co[::,::,1]-self["muy"])/(self["sigx"]*self["sigy"])) / \
            (-2.*(1.-self["rho"]**2)) )
    return result


Alternatively, you may find it more hand the coordinate axes directly to the
`evaluate` function. This is possible and demonstrated in the following
example, which relies on the :py:class:`GaussFit2dTuple` class, whose
`evaluate` method expects a tuple of coordinate axes.
Note the difference in the `evaluate` method. In this case, we
use numpy's `meshgrid` function to create coordinate mapping similar
to the one used before.

::

  from PyAstronomy import funcFit as fuf
  import numpy as np
  import matplotlib.pylab as plt
  
  # Constructing the individual coordinate axes
  x = np.linspace(-2.,2.,50)
  y = np.linspace(-2.,2.,50)
  
  # Create the 2d-Gaussian model and assign
  # some model parameters.
  gf = fuf.GaussFit2dTuple()
  gf["sigx"] = 0.75
  gf["sigy"] = 0.4
  gf["A"] = 1.0
  gf["rho"] = 0.4
  
  # Get the "data" by evaluating the model
  # and adding some noise. Note that the coordinate
  # mapping (array g) is passed to evaluate here.
  im = gf.evaluate((x,y))
  im += np.reshape(np.random.normal(0.0, 0.1, 2500), (50,50))
  err = np.ones((50,50))*0.1
  
  # Thaw parameters and fit
  gf.thaw(["A", "rho"])
  gf.fit((x,y), im, yerr=err)
  
  # Show the resulting parameter values ...
  gf.parameterSummary()
  
  # ... and plot the result.
  plt.title("Image data")
  plt.imshow(np.transpose(im), origin="lower")
  plt.show()
  plt.title("Residuals")
  plt.imshow(np.transpose(im - gf.evaluate((x,y))), origin="lower")
  plt.show() 

The associated `evaluate` method looks like this:

::

  # The "evaluate" method of class GaussFit2dTuple
  # A coordinate grid is constructed using numpy's meshgrid
  # function. However, you may do whatever you prefer to
  # map the input to the model coordinates. 

  def evaluate(self, co):
    """
      Evaluates the model for current parameter values.
      
      Parameters
      ----------
      co : array
           Specifies the points at which to evaluate the model.
    """
    if (self["sigx"] <= 0.0) or (self["sigy"] <= 0.0):
      raise(PE.PyAValError("Width(s) of Gaussian must be larger than zero.", \
                           solution="Change width ('sigx/y')."))
    if self["rho"] > 1.0:
      raise(PE.PyAValError("The correlation coefficient must be 0 <= rho <= 1.", \
                           solution="Change width ('sigx/y')."))
    xx, yy = meshgrid(co[0], co[1])
    result = 1.0/(2.*pi*self["sigx"]*self["sigy"]*sqrt(1.-self["rho"]**2)) * \
        exp( ((xx-self["mux"])**2/self["sigx"]**2 + (yy-self["muy"])**2/self["sigy"]**2 - \
            2.*self["rho"]*(xx-self["mux"])*(yy-self["muy"])/(self["sigx"]*self["sigy"])) / \
            (-2.*(1.-self["rho"]**2)) )
    return result
  
  