import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE

def gyroAgeBarnes(p, bv):
  """
    Calculate gyrochronological age according to Barnes 2007.
    
    The gyrochronological age is calculated according to Eq. 3 in
    Barnes 2007 (ApJ 669, 1167). The derivation of the error follows Eq. 16. 
    
    Parameters
    ----------
    p : float
        Stellar rotation period [d].
    bv : float
        B-V color [mag]. Supported range is 0.4-1.6 mag.
    
    Returns
    -------
    Stellar age : float
        The gyrochronological age [Ga].
    Age error : float
        The error on the age [Ga].
  """
  # See Eq. 3
  n = 0.5189
  a = 0.7725
  b = 0.601
  x = bv-0.4
  
  if bv < 0.4:
    raise(PE.PyAValError("Relation not defined for B-V < 0.4. Value of " + str(bv) + " was given.", \
                         where="gyroAgeBarnes", \
                         solution="Use a value between 0.4 and 1.6 for B-V."))
  if bv > 1.6:
    PE.warn(PE.PyAValError("No calibration for B-V > 1.6. Value given is " + str(bv) + "." +\
                           "The result is an extrapolation.", \
                           where="gyroAgeBarnes"))
  
  # Logarithm of the age [Myr]
  lnt = 1./n * (np.log(p) - np.log(a) - b*np.log(x))
  # Age in Myr
  t = np.exp(lnt)
  # Relative error (Eq. 16)
  dtt = 0.02*np.sqrt(3. + 0.5*lnt**2 + 2.0*p**0.6 + (0.6/x)**2 + (2.4*np.log(x))**2)
  return t/1000., t*dtt/1000.


def chromoAgeRHK(log10RHKprime):
  """
    Calculate the chromospheric age according to Donahue 1998.
    
    Donahue 1998 (ASPC 154, 1235) give a relation between chromospheric
    activity as measured by the R'HK index and the age of late-type stars
    (Eq. 1).
    
    As the level of stellar activity undergoes continuous change, Donahue
    cautions that individual measurements of the activity level yield
    notoriously poor age estimates. As an example, the spread in
    chromospheric solar age resulting from the 11 yr activity cycle is given,
    which amounts to about 2.5 Ga. This gives an idea of the accuracy of
    the estimates.
    
    Parameters
    ----------
    log10RHKprime : float
        Chromospheric activity index log10(R'HK).
    
    Returns
    -------
    Age : float
        Stellar age [Ga].
  """ 
  RHK = 10.0**log10RHKprime
  R5 = 1e5*RHK
  
  logAge = 10.725 - 1.334*R5 + 0.4085*R5**2 - 0.0522*R5**3
  age = 10.0**logAge / 1e9
  return age