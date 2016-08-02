# -*- coding: utf-8 -*-
from __future__ import division
import numpy
import PyAstronomy.funcFit as fuf
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.modelSuite.XTran import _ZList


# Check whether mpmath and/or the elliptical integrals from the boost library
# can be imported.

mpmath_imported = False
boostEll_imported = False

try:
  import mpmath
  mpmath_imported = True
except ImportError:
  mpmath_imported = False

try:
  import ell
  boostEll_imported = True
except ImportError:
  boostEll_imported = False


class _Case:
  """ 
    *Case-designation* according to Pal 2008 Table 1.
    
    Parameters
    ----------
    s : int
        Step identifier in Pal 2008 Table 1.
    c : string
        Case identifier in Pal 2008 Table 1.
  """

  def __init__(self, s, c):
    self.step = s; self.case = c
  
  def _str_(self):
    return "Step: "+str(self.step)+"  case: "+str(self.case)
    



class PalLC(_ZList, fuf.OneDFit):
  """
    Calculate and fit analytical transit light-curves using the formulae
    provided by Pal 2008.
    
    This class uses a **circular planetary orbit**.
    
    .. note :: The **evaluation of elliptical integrals** is essential
               in calculating the transit model. While both
               the *mpmath* module and the *Boost* libraries implement
               those integrals, it is the Boost library, which
               evaluates them by far more quickly. Yet, the support
               for Boost has to be added manually. 
    
    *Model parameters*:
      - `p` - Radius ratio between planet and star.
      - `a` - Semi-major axis of planetary orbit [stellar radii].
      - `i` - Inclination of orbit in degrees (90 deg is *edge on* view).
      - `linLib` - Linear limb-darkening coefficient.
      - `quadLimb` - Quadratic limb-darkening coefficient.
      - `T0` - Time offset of transit center.
      - `per` - Period of planetary orbit.
      - `b` - Describes the flux ratio between a stellar companion and the main star (default is 0).

    This class inherits the functionality of funcFit's OneDFit object.
    You can, therefore, set and get the parameter using the brackets:
      e.g., pallc["p"] = 0.12345

    .. warning::
        Time units have to be consistent.
  """
  
  def __init__(self):
    _ZList.__init__(self, "circular")
    fuf.OneDFit.__init__(self, ["p", "a", "i", "linLimb", "quadLimb", "T0", "per", "b"])
    self.freeze(["p", "a", "i", "linLimb", "quadLimb", "T0", "per", "b"])
    self.setRootName("PalCirc08")
    
    self._zlist=None
    
    if (not mpmath_imported) and (not boostEll_imported):
      raise(PE.PyARequiredImport("Neither mpmath nor the elliptical integrals from the boost library could be found!", \
                                 where="PalLC::__init__", solution="Install mpmath or make Boost library available (more complicated - see documentation)"))
    self.useBoost = False
    if boostEll_imported:
      self.useBoost = True
      self.ell = ell.ell()

  def _selectCases(self):
    """ 
      Determines which 'case' (according to Pal '08, Table 1) must be
      attributed to an individual value of 'z' and returns a
      list of 'Case' class objects.
    """
        
    result = []
    p = self["p"]
    for z in self._zlist:
      if z == 0.0 and p < 1.0:
        result.append(_Case(1,'A')); continue
      if z <= p-1.0:
        result.append(_Case(2,'Ag')); continue
      if z < p and z < (1.0-p):
        result.append(_Case(3,'B')); continue
      if z < p and z == 1.0-p:
        result.append(_Case(4,'Bt')); continue
      if z < p:
        result.append(_Case(5,'Bg')); continue
      if z == p and z < 1.0-p:
        result.append(_Case(6,'C')); continue
      if z == p and z == 0.5:
        result.append(_Case(7,'Ct')); continue
      if z == p:
        result.append(_Case(8,'Cg')); continue
      if z < 1.0-p:
        result.append(_Case(9,'D')); continue
      if z == (1.0-p):
        result.append(_Case(10,'E')); continue
      if z < 1.0+p:
        result.append(_Case(11,'F')); continue
      # This means that there is no occultation at all
      result.append(_Case(12,'G'))
    return result

  def _ci(self,z,a):
    return 2/(9.0*numpy.pi*numpy.sqrt(1.0-a))

  def _cik(self,z,a,b):
    return 1.0-5.0*z**2.0+self["p"]**2.0+a*b

  def _cie(self,z,a):
    return (z**2.0+7*self["p"]**2.0-4.0)*(1.0-a)

  def _cipi(self,z):
    return -3.0*(self["p"]+z)/(self["p"]-z)

  def _cg(self,z):
    return 1.0/(9.0*numpy.pi*numpy.sqrt(self["p"]*z))

  def _cgk(self,z):
    return 3.0-6.0*(1.0-self["p"]*self["p"])**2.0-2.0*self["p"]*z*(z*z+7*self["p"]*self["p"]-4.0+5.0*self["p"]*z)

  def _cge(self,z):
    return 4.0*self["p"]*z*(z*z+7.0*self["p"]*self["p"]-4.0)

  def _cgpi(self,z):
    return -3.0*(self["p"]+z)/(self["p"]-z)

  def _ti(self,z):
    phat=numpy.sqrt(self["p"]*(1.0-self["p"]))
    return 2.0/(3.0*numpy.pi)*numpy.arccos(1.0-2.0*self["p"])-4.0/(9.0*numpy.pi)*(3.0+2.0*self["p"]-8.0*self["p"]**2.0)*phat

  def _g0(self,z):
    k0=numpy.arccos((self["p"]**2.0+z**2.0-1.0)/(2.0*self["p"]*z))
    k1=numpy.arccos((z**2.0+1.0-self["p"]**2.0)/(2.0*z))
    return (self["p"]*self["p"]*k0+k1-numpy.sqrt(z*z-0.25*(1+z*z-self["p"]*self["p"])**2.0))/numpy.pi

  def _g2(self,z,a,b):
    k0=numpy.arccos((self["p"]**2.0+z**2.0-1.0)/(2.0*self["p"]*z))
    k1=numpy.arccos((z**2.0+1.0-self["p"]**2.0)/(2.0*z))
    return (k1+self["p"]*self["p"]*(self["p"]*self["p"]+2.0*z*z)*k0-0.25*(1.0+5.0*self["p"]*self["p"]+z*z)*numpy.sqrt((1.0-a)*(b-1.0)))/(2.0*numpy.pi)

  def _returnCoeff(self,step,z):
    a=(self["p"]-z)*(self["p"]-z)
    b=(self["p"]+z)*(self["p"]+z)
    if (step==1):
      pprime=numpy.sqrt(1.0-self["p"]**2.0)
      f0=self["p"]**2.0
      f1=2.0/3.0*(1.0-pprime**3.0)
      fk=0.0
      fe=0.0
      fpi=0.0
      f2=0.5*self["p"]**4.0
      k=0.0
      n=0.0
    elif (step==2):
      f0=1.0
      f1=2.0/3.0
      fk=0.0
      fe=0.0
      fpi=0.0
      f2=0.5
      k=0.0
      n=0.0
    elif (step==3):
      f0=self["p"]**2.0
      f1=2.0/3.0
      fk=self._ci(z,a)*self._cik(z,a,b)
      fe=self._ci(z,a)*self._cie(z,a)
      fpi=self._ci(z,a)*self._cipi(z)
      f2=0.5*self["p"]**2.0*(self["p"]**2.0+2.0*z**2.0)
      k=numpy.sqrt((4.0*self["p"]*z)/(1.0-a))
      n=-4.0*self["p"]*z/a
    elif (step==4):
      f0=self["p"]**2.0
      f1=self._ti(z)
      fk=0.0
      fe=0.0
      fpi=0.0
      f2=0.5*self["p"]**2.0*(self["p"]**2.0+2.0*z**2.0)
      k=0.0
      n=0.0
    elif (step==5):
      f0=self._g0(z)
      f1=2.0/3.0
      fk=self._cg(z)*self._cgk(z)
      fe=self._cg(z)*self._cge(z)
      fpi=self._cg(z)*self._cgpi(z)
      f2=self._g2(z,a,b)
      k=numpy.sqrt((1.0-a)/(4.0*self["p"]*z))
      n=(a-1.0)/a
    elif (step==6):
      f0=self["p"]**2.0
      f1=1.0/3.0
      fk=2.0/(9.0*numpy.pi)*(1.0-4.0*self["p"]**2.0)
      fe=8.0/(9.0*numpy.pi)*(2.0*self["p"]**2.0-1.0)
      fpi=0.0
      f2=3.0/2.0*self["p"]**4.0
      k=2*self["p"]
      n=0.0
    elif (step==7):
      f0=0.25
      f1=1.0/3.0-4.0/(9.0*numpy.pi)
      fk=0.0
      fe=0.0
      fpi=0.0
      f2=3.0/32.0
      k=0.0
      n=0.0
    elif (step==8):
      f0=self._g0(z)
      f1=1.0/3.0
      fk=-1.0/(9.0*numpy.pi*self["p"])*(1.0-4.0*self["p"]**2.0)*(3.0-8.0*self["p"]**2.0)
      fe=1.0/(9.0*numpy.pi)*16.0*self["p"]*(2.0*self["p"]**2.0-1.0)
      fpi=0.0
      f2=self._g2(z,a,b)
      k=1.0/(2.0*numpy.pi)
      n=0.0
    elif (step==9):
      f0=self["p"]**2.0
      f1=0.0
      fk=self._ci(z,a)*self._cik(z,a,b)
      fe=self._ci(z,a)*self._cie(z,a)
      fpi=self._ci(z,a)*self._cipi(z)
      f2=0.5*self["p"]**2.0*(self["p"]**2.0+2.0*z**2.0)
      k=numpy.sqrt((4.0*self["p"]*z)/(1.0-a))
      n=-4.0*self["p"]*z/a
    elif (step==10):
      f0=self["p"]**2.0
      f1=self._ti(z)
      fk=0.0
      fe=0.0
      fpi=0.0
      f2=0.5*self["p"]**2.0*(self["p"]**2.0+2.0*z**2.0)
      k=0.0
      n=0.0
    elif (step==11):
      f0=self._g0(z)
      f1=0.0
      fk=self._cg(z)*self._cgk(z)
      fe=self._cg(z)*self._cge(z)
      fpi=self._cg(z)*self._cgpi(z)
      f2=self._g2(z,a,b)
      k=numpy.sqrt((1.0-a)/(4.0*self["p"]*z))
      n=(a-1.0)/a
    elif (step==12):
      f0=0.0
      f1=0.0
      fk=0.0
      fe=0.0
      fpi=0.0
      f2=0.0
      k=0.0
      n=0.0
    return [f0,f1,fk,fe,fpi,f2,float(k),float(n)]

  def whichEllInts(self):
    """
      Says which module is used for the evaluation of the elliptical functions.
      
      .. note::
         The support for BOOST has to be added manually!
      
      Returns
      -------
      Identifier : string
          Either 'boost' or 'mpmath' depending on whether the *BOOST*
          libraries are used or the *mpmath* python module.
    """
    if self.useBoost:
      return "boost"
    return "mpmath"
    
  def evaluate(self, time):
    """ 
     Calculate a light curve according to the analytical models
     given by Pal 2008.
        
     Parameters
     ----------
     time : array
         An array of time points at which the light curve
         shall be calculated.
        
     .. note:: time = 0 -> Planet is exactly in the line of sight (phase = 0).

     Returns
     -------
     Model : array
         The analytical light curve is stored in the property `lightcurve`.
    """

    # Translate the given parameters into an orbit and, finally,
    # into a projected, normalized distance (z-parameter)
    self._calcZList(time - self["T0"])
  
    # 'W' parameters corresponding to notation in Pal '08
    w  = 6.-2.*self["linLimb"]-self["quadLimb"]
    w0 = (6.-6.*self["linLimb"]-12.*self["quadLimb"])/w
    w1 = (6.*self["linLimb"]+12.*self["quadLimb"])/w
    w2 = 6.*self["quadLimb"]/w
    
    # Initialize flux decrease array
    df = numpy.zeros(len(time))
    # Get a list of 'cases' (according to Pal '08). Depends on radius ratio and 'z'-parameter along the orbit
    ca = self._selectCases()
    # Loop through z-list, and calculate the light curve at each point in z (->time)
    for i in self._intrans:
      # Calculate the coefficients to be substituted into the Pal '08 equation
      c = self._returnCoeff(ca[i].step, self._zlist[i])
      # Substitute the coefficients and get 'flux decrease'
      if ca[i].step!=12:
        # Calculate flux decrease only if there is an occultation
        if not self.useBoost:
          df[i] = w0*c[0] + w2*c[5] + w1*(c[1] + c[2]*mpmath.ellipk(c[6]**2) + c[3]*mpmath.ellipe(c[6]**2) + c[4]*mpmath.ellippi(c[7],c[6]**2) )
        else:
          df[i] = w0*c[0] + w2*c[5] + w1*(c[1] + c[2]*self.ell.ell1(c[6]) + c[3]*self.ell.ell2(c[6]) + c[4]*self.ell.ell3(c[7],c[6]) )
    self.lightcurve = (1.-df)*1./(1.+self["b"]) + self["b"]/(1.0+self["b"])
    return self.lightcurve





    
class PalLCKep(PalLC):
  """
    Analytical transit light-curves using the formulae provided by Pal 2008.
    
    More information on the Keplerian orbit can be found here: :ref:`keplerorbitpyasl`
    
    .. note :: The **evaluation of elliptical integrals** is essential
               in calculating the transit model. While both
               the *mpmath* module and the *Boost* libraries implement
               those integrals, it is the Boost library, which
               evaluates them by far more quickly. Yet, the support
               for Boost has to be added manually. 
    
    *Model parameters*:
      - `p` - Radius ratio between planet and star.
      - `a` - Semi-major axis of planetary orbit [stellar radii].
      - `w` - Longitude of periastron [deg].
      - `Omega` - Longitude of the ascending node [deg].
      - `e` - Orbital eccentricity.
      - `i` - Inclination of orbit in degrees (90 deg is *edge on* view).
      - `linLib` - Linear limb-darkening coefficient.
      - `quadLimb` - Quadratic limb-darkening coefficient.
      - `tau` - Time of periastron passage.
      - `per` - Period of planetary orbit.
      - `b` - Describes the flux ratio between a stellar companion and the main star (default is 0).

    This class inherits the functionality of funcFit's OneDFit object.
    You can, therefore, set and get the parameter using the brackets:
      e.g., pallc["p"] = 0.12345

    .. warning::
        Time units have to be consistent.
    
    Parameters
    ----------
    collisionCheck : boolean
        If set True, it will be checked whether the two bodies collide on the current orbit.
  """
  
  def __init__(self, collisionCheck=False):
    _ZList.__init__(self, "keplerian", collisionCheck)
    fuf.OneDFit.__init__(self, ["p", "a", "i", "linLimb", "quadLimb", "tau", "per", "b", \
                                "w", "Omega", "e"])
    self.freeze(["p", "a", "i", "linLimb", "quadLimb", "tau", "per", "b", \
                 "w", "Omega", "e"])
    self.setRootName("Pal08")
    self["per"] = 1.0
    self["w"] = -90
    
    self._zlist=None
    
    if (not mpmath_imported) and (not boostEll_imported):
      raise(PE.PyARequiredImport("Neither mpmath nor the elliptical integrals from the boost library could be found!", \
                                 where="PalLC::__init__", solution="Install mpmath or make Boost library available (more complicated - see documentation)"))
    self.useBoost = False
    if boostEll_imported:
      self.useBoost = True
      self.ell = ell.ell()

    self.collisionCheck = collisionCheck
    
    # Wrap get/setitem to inform about change in parameter name
    T0paE = PE.PyADeprecationError("The parameter 'T0pa' in PalLCKep had to be renamed 'tau'.", \
                                   solution="Use 'tau' instead of 'T0pa'.")
    
    def getitem(specifier, **kwargs):
      if specifier == "T0pa":
        PE.warn(T0paE)
      return PalLC.__getitem__(self, specifier, **kwargs)
    self.__getitem__ = getitem
    
    def setitem(specifier, value):
      if specifier == "T0pa":
        PE.warn(T0paE)
      PalLC.__setitem__(self, specifier, value)
    self.__setitem__ = setitem

  def evaluate(self, time):
    """ 
      Calculate a light curve according to the analytical models
      given by Pal 2008.
         
      Parameters
      ----------
      time : array
          An array of time points at which the light curve
          shall be calculated.
               
      Returns
      -------
      Model : array
          The analytical light curve is stored in the property `lightcurve`.
    """

    # Translate the given parameters into an orbit and, finally,
    # into a projected, normalized distance (z-parameter)
    self._calcZList(time)
  
    # 'W' parameters corresponding to notation in Pal '08
    w  = 6.-2.*self["linLimb"]-self["quadLimb"]
    w0 = (6.-6.*self["linLimb"]-12.*self["quadLimb"])/w
    w1 = (6.*self["linLimb"]+12.*self["quadLimb"])/w
    w2 = 6.*self["quadLimb"]/w
    
    # Initialize flux decrease array
    df = numpy.zeros(len(time))
    # Get a list of 'cases' (according to Pal '08). Depends on radius ratio and 'z'-parameter along the orbit
    ca = self._selectCases()
    # Loop through in-transit points in z-list,
    # and calculate the light curve at each point in z (->time)
    for i in self._intrans:
      # Calculate the coefficients to be substituted into the Pal '08 equation
      c = self._returnCoeff(ca[i].step, self._zlist[i])
      # Substitute the coefficients and get 'flux decrease'
      if ca[i].step!=12:
        # Calculate flux decrease only if there is an occultation
        if not self.useBoost:
          df[i] = w0*c[0] + w2*c[5] + w1*(c[1] + c[2]*mpmath.ellipk(c[6]**2) + c[3]*mpmath.ellipe(c[6]**2) + c[4]*mpmath.ellippi(c[7],c[6]**2) )
        else:
          df[i] = w0*c[0] + w2*c[5] + w1*(c[1] + c[2]*self.ell.ell1(c[6]) + c[3]*self.ell.ell2(c[6]) + c[4]*self.ell.ell3(c[7],c[6]) )
    self.lightcurve = (1.-df)*1./(1.+self["b"]) + self["b"]/(1.0+self["b"])
    return self.lightcurve


PalLC_Rebin = fuf.turnIntoRebin(PalLC)