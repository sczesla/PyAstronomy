from PyAstronomy.pyasl import KeplerEllipse
from PyAstronomy import funcFit as fuf
from numpy import pi

class KeplerEllipseModel(fuf.OneDFit):
  """
    A model of a Keplerian orbit.
    
    This class uses the *KeplerEllipse* from the PyA's pyasl
    to calculate a Keplerian orbit. It may be used to fit
    complete 3d information on the orbit. Projections using
    only one or two dimensions are also possible.
    
    The constructor allows to specify *relevant axes*, which
    are those axes considered in the calculation. The actual
    model is, however, only one-dimensional. The values
    returned by evaluate have the order
    a1, b1, c1, a2, b2, c3, ... . Where a, b, and c stand for
    the first, second, and third axis and the number specifies
    the data point. Note that in this case, the resulting
    model has not the same number of points as the time
    axis.
    
    *Fit parameters*
      -  `a`    - The semi-major axis (same units as the data)
      - `per`   - The period (same time units as data)
      - `e`     - The eccentricity
      - `tau`   - Time of periastron passage (same time units as data)
      - `Omega` - Longitude of the ascending node [deg]
      - `w`     - Longitude of periastron [deg]
      - `i`     - Inclination angle [deg]
    
    Parameters
    ----------
    relevantAxes : string
        A string containing any combination of x, y, and z.
        The string specifies the axes (and their order) to
        be considered in the calculations.
  """

  def __init__(self, relevantAxes="xyz"):
    self.ke = KeplerEllipse(1.0, 1.0)
    fuf.OneDFit.__init__(self, ["a", "per", "e", "tau", "Omega", "w", "i"])
    self["a"] = 1.0
    self["per"] = 1.0
    # Which axes to consider?
    # x=0, y=1, z=2
    self.axes = ()
    for i, axis in enumerate("xyz"):
      if axis in relevantAxes:
        self.axes += (i,)

  def evaluate(self, t):
    """
      Calculates and returns model according to the
      current parameter values.
      
      Although more than one axis may be relevant
      the output will be one dimensional. If, e.g.,
      the relevant axes are x and y, the order of
      the output will be x0, y0, x1, y1, ... .

      Parameters
      ----------
      t : array
          Times at which to evaluate the model.
    """
    self.ke.i = self["i"]
    self.ke.w = self["w"]
    self.ke.Omega = self["Omega"]
    self.ke.e = abs(self["e"])
    self.ke.a = self["a"]
    self.ke.per = self["per"]
    self.ke.tau = self["tau"]
    self.ke.n = 2.0*pi/self["per"]
    # Get the data pertaining to the relevant axes
    pos = self.ke.xyzPos(t)[::,self.axes]
    # Reshape to 1d form
    # If the relevant axes are x and y, the order
    # of the output will be x0, y0, x1, y1, x2, y2, ...
    pos = pos.reshape(pos.size)
    return pos

