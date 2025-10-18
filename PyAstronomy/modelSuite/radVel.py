import numpy
from PyAstronomy.funcFit import OneDFit

class SinRadVel(OneDFit):
  """
    Spectroscopic radial-velocity (RV) shift due to circular orbital motion.

    *Fit parameters*:
     - P       - float, Orbital period [d]
     - T0      - float, Central transit time [d]
     - K       - float, radial velocity semi-amplitude [km/s]
     - rv0     - float, constant offset in radial velocity [km/s]

   By default all parameters remain frozen.

   The convention is that at phase zero, also the orbital
   radial velocity is zero. With increasing phase the object
   becomes bluer first.
  """

  def __init__(self):
    OneDFit.__init__(self,["P", "T0", "K", "rv0"])
    self.setRootName("SinRV")



  def evaluate(self, x):
    """
      Calculates and returns radial velocity shift according to
      current model parameters.

      Parameters
      ----------
      x : array
          The time stamps at which to calculate the model RV curve.
    """
    y = -numpy.sin(2.0*numpy.pi*(x-self["T0"])/self["P"]) * self["K"] + self["rv0"]
    return y


