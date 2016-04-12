# -*- coding: utf-8 -*-
from __future__ import division
import numpy
from numpy import sqrt, sin, cos, pi, exp
from numpy import arcsin as asin
from numpy import arccos as acos
from scipy.optimize import fsolve
from scipy.special import i0,i1,erf
from PyAstronomy.funcFit import OneDFit, turnIntoRebin
from PyAstronomy.modelSuite.XTran import palTrans


class RmcL(OneDFit):
  """
    Analytical Rossiter-McLaughlin effect.
  
    This class implements the analytical model radial velocity (RV) curves
    for the Rossiter-McLaughlin effect given by *Ohta et. al 2005*.
  
    *Fit parameters*:
     - epsilon - linear limb dark
     - gamma   - Rp/Rs (ratio of planetary and stellar radius)
     - P       - Orbital period [d]
     - T0      - Central transit time
     - i       - Inclination of orbit [rad]
     - Is      - Inclination of stellar rotation axis [rad]
     - Omega   - Angular rotation velocity (star) [rad/s]
     - lambda  - Sky-projected angle between stellar
       rotation axis and normal of orbit plane [rad]
     - a       - Semi major axis [stellar radii]
   
   By default all parameters remain frozen.
   
   .. note::
     According to the input parameter units, the units of the model RV
     curve are **stellar-radii per second**.
  """

  def __init__(self):
    OneDFit.__init__(self,["epsilon", "gamma", "P", "T0", "i", "Is", "Omega", "lambda", "a"])
    self.setRootName("Ohta05")

  def planetDistance(self, f):
    return self["a"]

  def W1(self, rho):
    result = sqrt(1.0 - rho**2) - self["gamma"]**2 * (2.0-rho**2)/(8.0*(1.0-rho**2)**(3.0/2.0))
    return result
  
  def W2(self, rho):
    result =  sqrt(1.0 - rho**2) - self["gamma"]**2 * (4.0-3.0*rho**2)/(8.0*(1.0-rho**2)**(3.0/2.0))
    return result

  def XpVec(self, f):
    result = numpy.zeros(3)
    result[0] = -cos(self["lambda"])*sin(f) - sin(self["lambda"])*cos(self["i"])*cos(f)
    result[1] = sin(self["i"])*cos(f)
    result[2] = sin(self["lambda"]*sin(f)) - cos(self["lambda"])*cos(self["i"])*cos(f)
    result *= self.rp(f)
    return result
  
  def Xp(self, f):
    result = -cos(self["lambda"])*sin(f) - sin(self["lambda"])*cos(self["i"])*cos(f)
    result *= self.planetDistance(f)
    return result

  def Zp(self, f):
    result = sin(self["lambda"])*sin(f) - cos(self["lambda"])*cos(self["i"])*cos(f)
    result *= self.planetDistance(f)
    return result

  def etap(self, Xp, Zp):
    return sqrt(Xp**2 + Zp**2) - 1.0
  
  def zeta(self, etap):
    result = (2.0*etap + self["gamma"]**2 + etap**2) / (2.0*(1.0+etap))
    return result

  def rhoFromVec(self, XpVec):
    return sqrt(XpVec[0]**2 + XpVec[2]**2)

  def rho(self, Xp, Zp):
    return sqrt(Xp**2 + Zp**2)

  def trueAnomaly(self, time):
    result = ((time-self["T0"])/self["P"] - numpy.floor((time-self["T0"])/self["P"])) * 2.0*numpy.pi
    return result

  def z0(self, etap, indi):
    result = numpy.zeros(etap.size)
    result[indi] = sqrt((self["gamma"]**2 - etap[indi]**2)*( (etap[indi]+2.0)**2-self["gamma"]**2) ) / \
             (2.0*(1.0+etap[indi]))
    return result
  
  def x0(self, etap):
    return 1.0-(self["gamma"]**2 - etap**2)/(2.0*(1.0+etap))

  def g(self, x, e, g, x0):
    result = (1.0-x**2) * asin(sqrt((g**2-(x-1.0-e)**2)/(1.0-x**2))) + \
       sqrt(2.0*(1.0+e)*(x0-x)*(g**2-(x-1.0-e)**2))
    return result
  
  def xc(self, zeta, x0):
    return x0+(zeta-self["gamma"])/2.0

  def W3(self, x0, zeta, xc, etap):
    result = pi/6.0*(1.0-x0)**2*(2.0+x0) + \
      pi/2.0*self["gamma"]*(self["gamma"]-zeta) * \
      self.g(xc, etap, self["gamma"], x0) / self.g((1.0-self["gamma"]), -self["gamma"], self["gamma"], x0) * \
      self.W1(1.0-self["gamma"])
    return result
  
  def W4(self, x0, zeta, xc, etap):
    result = pi/8.*(1.0-x0)**2*(1.0+x0)**2 + \
      pi/2.*self["gamma"]*(self["gamma"]-zeta)*xc * \
      self.g(xc, etap, self["gamma"], x0) / self.g((1.0-self["gamma"]), -self["gamma"], self["gamma"], x0) * \
      self.W2(1.0-self["gamma"])
    return result

  def evaluate(self, xOrig):
    """
      Calculates and returns RV curve according to current model parameters.
      
      .. note:: The units of the model RV curve are **stellar-radii per second**.
      
      Parameters
      ----------
      xOrig : array
          The time stamps at which to calculate the model RV curve.
          Note that the orbit period and central transit time are used
          to convert time into "true anomaly".
    """
    x = self.trueAnomaly(xOrig)
    Xp = self.Xp(x)
    Zp = self.Zp(x)
    rho = self.rho(Xp, Zp)
    etap = self.etap(Xp, Zp)
    zeta = self.zeta(etap)
    x0 = self.x0(etap)
    xc = self.xc(zeta, x0)

    # dphase is the phase difference between the primary transit and the time points
    # It is used to exclude the secondary transit from the calculations 
    dphase = numpy.abs((xOrig-self["T0"])/self["P"])
    dphase = numpy.minimum( dphase-numpy.floor(dphase), numpy.abs(dphase-numpy.floor(dphase)-1))

    y = numpy.zeros(len(x))
    indi = numpy.where(numpy.logical_and(rho < (1.0-self["gamma"]), dphase < 0.25))[0]

    y[indi] = Xp[indi]*self["Omega"]*sin(self["Is"])* self["gamma"]**2 * \
        (1.0 - self["epsilon"]*(1.0 - self.W2(rho[indi]))) / \
        (1.0 - self["gamma"]**2 - self["epsilon"]*(1./3. - self["gamma"]**2*(1.0-self.W1(rho[indi]))))
    
    indi = numpy.where(numpy.logical_and( \
                numpy.logical_and(rho >= 1.-self["gamma"], rho < 1.0+self["gamma"]), dphase < 0.25))[0]
    z0 = self.z0(etap, indi)
    
    y[indi] = (Xp[indi]*self["Omega"]*sin(self["Is"])*( \
        (1.0-self["epsilon"]) * (-z0[indi]*zeta[indi] + self["gamma"]**2*acos(zeta[indi]/self["gamma"])) + \
        (self["epsilon"]/(1.0+etap[indi]))*self.W4(x0[indi], zeta[indi], xc[indi], etap[indi]))) / \
        (pi*(1.-1.0/3.0*self["epsilon"]) - (1.0-self["epsilon"]) * (asin(z0[indi])-(1.+etap[indi])*z0[indi] + \
        self["gamma"]**2*acos(zeta[indi]/self["gamma"])) - self["epsilon"]*self.W3(x0[indi], zeta[indi], xc[indi], etap[indi]))

    return y


RmcL_Rebin = turnIntoRebin(RmcL)


class RmcL_Hirano(OneDFit):

  def __init__(self):
    """
      This class implements analytical expressions for the Rossiter-McLaughlin \
      effect for the case when the anomalous radial velocity is obtained by \
      cross-correlation with a stellar spectrum, according to *Hirano et. al 2010*.

      .. note::
        The intrinsic line profile and the rotation kernel are both \
        approximated by Gaussians, while the planet is assumed to be \
        'sufficiently small enough'.

      *Fit parameters*:
       - linLimb - linear limb-darkening parameter
       - quadLimb - quadratic limb-darkening parameter
       - gamma   - Rp/Rs (ratio of planetary and stellar radius)
       - P       - Orbital period [d]
       - T0      - Central transit time
       - i       - Inclination of orbit [rad]
       - Is      - Inclination of stellar rotation axis [rad]
       - Omega   - Angular rotation velocity (star) [rad/s]
       - lambda  - Sky-projected angle between stellar rotation axis and normal of orbit plane [rad]
       - a       - Semi major axis [stellar radii]
       - vbeta   - natural line broadening [km/s]
       - vSurf   - maximum surface velocity [km/s]

     By default all parameters remain frozen.

     .. note::
       According to the input parameter units, the units of the model RV curve are **stellar-radii per second**.
       The original formulation by Hirano et al. uses two dispersion parameters, \
       `beta` and `sigma`. `beta` describes the width of an intrinsic line profile \
       modeled by a Gaussian; it corresponds to a 'broadening velocity' \
       vbeta = beta * c / lam0, where lam0 denotes a typical wavelength scale. \
       The dispersion parameter `sigma` describes line broadening due to stellar \
       rotation; it corresponds to a velocity of vSurf*sin(Is) = alpha*sigma*c/lam0, \
       where lam0 denotes the typical wavelength and alpha is a scaling \
       parameters depending on the limb-darkening coefficients.

    """
    OneDFit.__init__(self,["linLimb", "quadLimb", "gamma", "P", "T0", "i", "Is", "Omega", "lambda", "a", "vbeta", "vSurf"])
    self.flux = None
    self.setRootName("Hira10")

  def planetDistance(self, f):
    return self["a"]

  def Xp(self, f):
    result = -cos(self["lambda"])*sin(f) - sin(self["lambda"])*cos(self["i"])*cos(f)
    result *= self.planetDistance(f)
    return result

  def trueAnomaly(self, time):
    result = ((time-self["T0"])/self["P"] - numpy.floor((time-self["T0"])/self["P"])) * 2.0*numpy.pi
    return result

  def supplyFlux(self, flux):
    """
      This method can be used to supply own simultaneously obtained photometric data.
      The continuum flux should be normalized to one.
    """
    self.flux = flux
    return 0

  def EqF6(self,x):
    """
      Relation between rotation and Gaussian broadening kernels by least-squares fitting
    """
    result = 1./sqrt(2.)-(12.*(self["linLimb"]+2.*self["quadLimb"])*exp(-x**2))/(x**2*(-6.+2.*self["linLimb"]+self["quadLimb"])) \
           + (6.*exp(-x**2/2.))/(x**3*(-6.+2.*self["linLimb"]+self["quadLimb"])) * ( -2.*x**3*(-1.+self["linLimb"]+self["quadLimb"])*i0(x**2/2.) \
           + 2.*x*(-2.*self["quadLimb"]+x**2*(-1.+self["linLimb"]+self["quadLimb"])) * i1(x**2/2.)
           + sqrt(pi)*exp(x**2/2.)*(self["linLimb"]+2.*self["quadLimb"])*erf(x) )
    return result

  def evaluate(self, time):
    """
      Calculates and returns RV curve according to current model parameters.

      .. note:: The units of the model RV curve are **stellar-radii per second**.

      Parameters
      ----------
      time : array
          The time stamps at which to calculate the model RV curve.
          Note that the orbit period and central transit time are
          used to convert time into "true anomaly".
    """
    x = self.trueAnomaly(time)
    Xp = self.Xp(x)

    if self.flux != None:
      f=1.-self.flux
    else:
      plc = palTrans.PalLC()
      plc.assignValue({"p": self["gamma"], \
                       "per": self["P"], \
                       "a": self["a"], \
                       "i": self["i"]*180./numpy.pi, \
                       "linLimb": self["linLimb"], \
                       "quadLimb": self["quadLimb"], \
                       "T0": self["T0"], \
                       "b": 0.0})
      f=1.-plc.evaluate(time)

    deltaV = Xp*self["Omega"]*sin(self["Is"])* f

    # Calculate scaling factor alpha, Eq. F6 in Hirano et al. 2010
    alpha = fsolve(self.EqF6, x0=1.0, xtol=1e-8)[0]

    # Determine effective rotational velocity
    vsiniEff = self["vSurf"]*sin(self["Is"]) / alpha

    # Apply the Hirano et al. 2010 correction factor, Eq. 37 in Hirano et al. 2010
    deltaV = deltaV * ((2.*self["vbeta"]**2+2.*vsiniEff**2)/(2.*self["vbeta"]**2+vsiniEff**2))**(3./2.) * (1.-( Xp*self["Omega"]*sin(self["Is"]))**2/(2.*self["vbeta"]**2+vsiniEff**2))
    return deltaV

