from __future__ import print_function, division
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyasl import AtomicNo
import os
from PyAstronomy.pyasl import _ic
import six

class FirstIonizationPot:
  """
    First ionization potentials of individual elements.
    
    The data have been obtained from NIST:
    http://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
  """

  def _decomposeFIP(self, s):
    """
      Decompose value string.
      
      Parameters
      ----------
      s : string
          The specification of the value. Values in
          round parentheses denote theoretical values,
          values in brackets are interpolated.
      
      Returns
      -------
      Value : float
          The first ionization energy [eV].
      Error: float
          The error of the first ionization energy [eV].
      isTheoretical : boolean
          True if value is theoretical.
      isInterpolated : boolean
          True if value has been interpolated.
    """
    isTheoretical = False
    isInterpolated = False
    
    s = s.strip()
    if s.startswith("("):
      isTheoretical = True
      s = s[1:-1]
    if s.startswith("["):
      isInterpolated = True
      s = s[1:-1]
    # Isolate error
    if s.find(")") != -1 :
      # There is an error given
      s = s.rstrip(")")
      r = s.split("(")
      error = float(r[1])
      # Determine number of significant digits
      ndig = len(r[0].split(".")[1])
      error *= 10.0**(-ndig)
    else:
      # No error was given
      error = None
      r = [s]
    val = float(r[0])
    return val, error, isTheoretical, isInterpolated

  def getFIP(self, atom):
    """
      Get the first ionization energy.
      
      Parameters
      ----------
      atom : string or integer
          Either a string specifying the elemental symbol (e.g., 'H') or
          an integer specifying the atomic number.
      
      Returns
      -------
      FIP : float
          First ionization potential [eV].
      FIP error : float
          Error of the FIP [eV]. May also be None
          if no error is available.
    """
    if isinstance(atom, six.string_types):
      an = self._an.getAtomicNo(atom)
    elif isinstance(atom, int):
      an = atom
    else:
      raise(PE.PyAValError("'atom' must be an integer specifying the atomic number are a string holding the elemental symbol.", \
                           where="FirstIonizationPot::getFIP"))
    return self._fip[an][0], self._fip[an][1]
    

  def __init__(self):
    self._fip = {}
    # Build FIP-data dictionary. Uses elemental number as key
    # and, for each, holds a tuple composed of value, error, a flag
    # describing whether the value is theoretical, and a flag specifying
    # whether the value has been interpolated (see NIST link for further)
    # details.
    import PyAstronomy as _PyA
    fn = os.path.join(os.path.dirname(_PyA.__file__), "pyasl", "resBased", "fip.dat")
    for l in open(fn):
      if l.startswith("#"):
        continue
      s = l.split("|")
      self._fip[int(s[0])] = self._decomposeFIP(s[8])
    # Convert between atomic number and elemental symbol
    self._an = AtomicNo()


def plotFIP():
  """
    Show a plot of first ionization energy vs. atomic number. 
  """
  if not _ic.check["matplotlib"]:
    raise(PE.PyARequiredImport("Could not import matplotlib."))
  
  fip = FirstIonizationPot()
  z = range(1,31,1)
  fips = []
  for an in z:
    fips.append(fip.getFIP(an)[0])
  
  import matplotlib.pylab as plt
  plt.xlabel("Atomic number")
  plt.ylabel("First ionization energy [eV]")
  plt.plot(z, fips, 'bp-')
  plt.show()
  
  
  