# -*- coding: utf-8 -*-
from __future__ import division
from PyAstronomy.funcFit import OneDFit
import numpy
from PyAstronomy.modelSuite.XTran import _ZList

class LimBrightTrans(_ZList, OneDFit):
  """
    Planetary transit light-curves for spherical shell model.
  
    This class implements a model calculating the light curve of a planet
    transiting an optically thin spherical shell of negligible thickness
    (e.g., a stellar chromosphere).
        
    The model provided by Schlawin et al. 2010 assumes that the thickness
    of the shell is much smaller than the size of the planet.
    The shell is optically thin and thus provides natural limb-brightening.
    The obscured part of the stellar surface is calculated based on computing the
    volume of the intersection of a sphere with a cylinder and then
    taking a partial derivative with respect to the radius of the
    sphere to find its surface area.    
    
    The code closely follows the IDL procedure located at \
    http://www.astro.washington.edu/agol/.
  
    *Fit parameters*:
     - `p`     - Rp/Rs (ratio of planetary and stellar radius)
     - `a`     - Semi-major axis of planetary orbit [stellar radii].
  
     - `per`   - Orbital period [d]
     - `T0`    - Central transit time
     - `i`     - Inclination of orbit [rad]
   
    By default all parameters remain frozen.
  """

  def __init__(self):
    _ZList.__init__(self, "circular")
    OneDFit.__init__(self, ["p", "a", "i", "T0", "per"])
    self.freeze(["p", "a", "i", "T0", "per"])

    self._zlist=None

  def __ell1(self,k):
    """
      Computes polynomial approximation for the complete elliptic
      integral of the first kind (Hasting's approximation)
    """
    m1=1.-k**2
    a0=1.38629436112
    a1=0.09666344259
    a2=0.03590092383
    a3=0.03742563713
    a4=0.01451196212
    b0=0.5
    b1=0.12498593597
    b2=0.06880248576
    b3=0.03328355346
    b4=0.00441787012
    ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
    ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*numpy.log(m1)
    return ek1-ek2

  def __ell2(self,k):
    """
      Computes polynomial approximation for the complete elliptic
      integral of the second kind (Hasting's approximation)
    """
    m1=1.-k**2
    a1=0.44325141463
    a2=0.06260601220
    a3=0.04757383546
    a4=0.01736506451
    b1=0.24998368310
    b2=0.09200180037
    b3=0.04069697526
    b4=0.00526449639
    ee1=1.+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
    ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*numpy.log(1./m1)
    return ee1+ee2

  def __ell3(self,n,k):
    """
      Computes the complete elliptical integral of the third kind using
      the algorithm of Bulirsch (1965)
    """
    kc=numpy.sqrt(1.-k**2.); p=n+1.
    if numpy.min(p) < 0.:
      print('Negative p')
    m0=1.;  c=1.;  p=numpy.sqrt(p);  d=1./p;  e=kc
    loop=True
    while loop:
      f = c;  c = d/p+f;  g = e/p;  d = (f*g+d)*2.
      p = g + p; g = m0; m0 = kc + m0
      if numpy.max(numpy.abs(1.-kc/g)) > 1e-13:
        kc = 2.*numpy.sqrt(e); e=kc*m0
      else:
        loop=False
    return 0.5*numpy.pi*(c*m0+d)/(m0*(m0+p))

  def evaluate(self, time):
    """
      Calculate a light curve according to the analytical model
      by Schlawin et al. 2010.

      Parameters
      ----------
      time : array
          An array of time points at which the light curve shall be calculated

     Returns
     -------
     Model : array
         The analytical light curve is stored in the property `lightcurve`.
    
     Notes
     -----
     
     .. note:: time = 0 -> Planet is exactly in the line of sight (phase = 0).
    """
    self._calcZList(time - self["T0"])

    a=numpy.zeros(len(self._zlist))
    indi=numpy.where(self._zlist+self["p"] < 1.)[0]
    if len(indi)>0:
      k=numpy.sqrt(4.*self._zlist[indi]*self["p"]/(1.-(self._zlist[indi]-self["p"])**2))
      a[indi]=4./numpy.sqrt(1.-(self._zlist[indi]-self["p"])**2)*(((self._zlist[indi]-self["p"])**2-1.)*self.__ell2(k) \
                -(self._zlist[indi]**2-self["p"]**2)*self.__ell1(k)+(self._zlist[indi]+self["p"])/(self._zlist[indi]-self["p"]) \
                *self.__ell3(4.*self._zlist[indi]*self["p"]/(self._zlist[indi]-self["p"])**2,k))

    indi=numpy.where(numpy.logical_and(self._zlist+self["p"] > 1.,self._zlist-self["p"] < 1.))[0]
    if len(indi)>0:
      k=numpy.sqrt((1.-(self._zlist[indi]-self["p"])**2)/4./self._zlist[indi]/self["p"])
      a[indi]=2./(self._zlist[indi]-self["p"])/numpy.sqrt(self._zlist[indi]*self["p"])*(4.*self._zlist[indi]*self["p"]*(self["p"]-self._zlist[indi])*self.__ell2(k) \
              +(-self._zlist[indi]+2.*self._zlist[indi]**2*self["p"]+self["p"]-2.*self["p"]**3)*self.__ell1(k) \
              +(self._zlist[indi]+self["p"])*self.__ell3(-1.+1./(self._zlist[indi]-self["p"])**2,k))

    self.lightcurve=1.-(4.*numpy.pi*(self["p"] > self._zlist)*1.+a)/4./numpy.pi
    return self.lightcurve
