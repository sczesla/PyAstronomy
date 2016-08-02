from __future__ import print_function, division
from PyAstronomy import funcFit as fuf
import numpy as np

class AtanProfile(fuf.OneDFit):
  """
    A profile based on the arc tangent function.
    
    This class implements the following profile:
    
    .. math::
      f(x) = \\frac{A}{2\\arctan(\\sigma)} \\times \\left(\\arctan\\left(\\frac{x-\mu}{scale} + \sigma\\right) +
             \\arctan\\left(-\\frac{x-\mu}{scale} + \sigma\\right)\\right) +
             \mu \\times x + off
    
    which can provide a relatively flat top and steep edges.
    
    *Fit parameters*
      -  `A`    - The amplitude. In this case, the height (not the area under)
                  the profile reached for :math:`x=0`. Note that for
                  :math:`\mu \\not = 0` the highest point may be elsewhere,
                  which is neglected here.
      - `scale` - A scale parameter affecting the width of the profile. Note,
                  however, that also :math:`\sigma` affects the width.
      - `mu`    - The center of the profile.
      - `off`   - An offset
      - `lin`   - A gradient in the offset.
    
    The width of the profile may be approximated by the inflection points, which
    are given by
    
    .. math::
      \\frac{\\partial^2 f(x)}{\partial x^2} = 0 \\rightarrow
      x_{1,2} = \mu \\pm\\frac{scale}{3}\\left(-3+3\sigma^2+6\\sqrt{\sigma^4+\sigma^2+1}\\right)^{1/2}
    
  """

  def __init__(self):
    fuf.OneDFit.__init__(self, ["scale", "sig", "mu", "A", "off", "lin"])

  def evaluate(self, x):
    """
      Calculates and returns model according to the
      current parameter values.

      Parameters
      ----------
      x : Array
          The positions at which to evaluate the model.
    """
    # Shift by mu
    x = x - self["mu"]
    # The heart of the profile
    y = np.arctan(x/self["scale"] + self["sig"]) + np.arctan(-x/self["scale"] + self["sig"])
    # Make the highest point (actually most extreme point)
    # equal to A
    y *= (self["A"] / (2.*np.arctan(self["sig"])))
    
    # Add offset and gradient
    y += self["off"]
    y += self["lin"] * (x + self["mu"])
    return y
 
  def inflectionPoints(self):
    """
      Calculate the inflection points.
      
      The inflection points of the profile depend on
      both :math:`\sigma` and :math:`\mu`.
      
      Returns
      -------
      Inflection points : tuple
          Locations of the inflection points. Smaller one first.
    """
    d = abs(self["scale"])/3.0 * \
        np.sqrt(-3. + 3.*self["sig"]**2 + 6.*np.sqrt(self["sig"]**4 + self["sig"]**2 + 1.0))
    return self["mu"]-d, self["mu"]+d



class AtanProfileDamped(fuf.OneDFit):
  """
    A profile based on the arc tangent function.
    
    This class implements the following profile:
    
    .. math::
      d(x) = f(x) \\times H(|x-\mu| - |ifp-\mu|) \\times
             \\exp\\left(\\frac{|x-\mu| - |ifp-\mu|}{\\tau}\\right) +
             \mu \\times x + off
    
    Here :math:`f(x)` is the profile described in :py:class:`AtanProfile`,
    H denotes the Heaviside function, and ifp is the location of the
    inflection point. The parameter :math:`\\tau` can be used to provide
    an additional drop at the edges of the profile.  
    
    
    *Fit parameters*
      -  `A`    - The amplitude. In this case, the height (not the area under)
                  the profile reached for :math:`x=0`. Note that for
                  :math:`\mu \\not = 0` the highest point may be elsewhere,
                  which is neglected here.
      - `scale` - A scale parameter affecting the width of the profile. Note,
                  however, that also :math:`\sigma` affects the width.
      - `tau`   - This parameter controls an additional drop at the edges
                  of the profile.
      - `mu`    - The center of the profile.
      - `off`   - An offset
      - `lin`   - A gradient in the offset.
    
    The width of the profile may be approximated by the inflection points, which
    are given by
    
    .. math::
      \\frac{\\partial^2 f(x)}{\partial x^2} = 0 \\rightarrow
      x_{1,2} = \mu \\pm\\frac{scale}{3}\\left(-3+3\sigma^2+6\\sqrt{\sigma^4+\sigma^2+1}\\right)^{1/2}
    
  """

  def __init__(self):
    fuf.OneDFit.__init__(self, ["scale", "sig", "mu", "A", "off", "lin", "tau"])

  def evaluate(self, x):
    """
      Calculates and returns model according to the
      current parameter values.

      Parameters
      ----------
      x : Array
          The positions at which to evaluate the model.
    """
    # Shift by mu
    x = x - self["mu"]
    # The heart of the profile
    y = np.arctan(x/self["scale"] + self["sig"]) + np.arctan(-x/self["scale"] + self["sig"])
    # Make the highest point (actually most extreme point)
    # equal to A
    y *= (self["A"] / (2.*np.arctan(self["sig"])))
    
    # Produce additional drop
    difp = abs(self.inflectionPoints()[0] - self["mu"])
    indi = np.where(np.abs(x) > difp)[0]
    y[indi] *= np.exp(-np.abs(np.abs(x[indi])-difp)**2/self["tau"])
    
    # Add offset and gradient
    y += self["off"]
    y += self["lin"] * (x + self["mu"])
    return y
 
  def inflectionPoints(self):
    """
      Calculate the inflection points.
      
      The inflection points of the profile depend on
      both :math:`\sigma` and :math:`\mu`.
      
      Returns
      -------
      Inflection points : tuple
          Locations of the inflection points. Smaller one first.
    """
    d = abs(self["scale"])/3.0 * \
        np.sqrt(-3. + 3.*self["sig"]**2 + 6.*np.sqrt(self["sig"]**4 + self["sig"]**2 + 1.0))
    return self["mu"]-d, self["mu"]+d