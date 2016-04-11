from __future__ import print_function, division
import unittest
from PyAstronomy import pyasl

class SanityOfPyaslExt1(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_pizzolatoExample(self):
    """
      Example of Pizzolato 2003 relations
    """
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    import numpy as np
    
    p = pyasl.Pizzolato2003()
    
    # Define array of rotation periods [days]
    prot = np.arange(0.2, 30, 0.1)
    
    lx = np.zeros(prot.size)
    lxlbol = np.zeros(prot.size)
    
    # B-V color of star
    bv = 0.7
    
    # Obtain ...
    for i in range(prot.size):
      # ... log10 of X-ray luminosity
      lx[i] = p.log10lxbv(bv, prot[i])[0]
      # ... and log10(Lx/Lbol)
      lxlbol[i] = p.log10lxlbolbv(bv, prot[i])[0]
    
    # Plot result
    plt.subplot(2,1,1)
    plt.plot(prot, lx, 'bp-')
    plt.subplot(2,1,2)
    plt.plot(prot, lxlbol, 'bp-')
#     plt.show()

  def sanity_pizzolato(self):
    """
      Pizzolato 2003 sanity
    """
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    import numpy as np
    
    p = pyasl.Pizzolato2003()
    
    x = p.log10lxbv(0.87, 1.0)
    self.assertAlmostEqual(x[0], 29.9, msg="Lx saturation level for b-v=0.87 does not match.", delta=1e-7)
    self.assertAlmostEqual(x[1], 0.3, msg="Lx error for saturation level for b-v=0.87 does not match.", delta=1e-7)
    
    x = p.log10lxlbolmass(1.15, 3.0)
    self.assertAlmostEqual(x[0], -4.1694, msg="Pizzolato relation (m=1.15, pr=3.0) failed (value="+str(x[0])+").", delta=1e-4)
    
  def sanity_expCorrRN_Example1(self):
    """
      Sanity of example 1 for expCorrRN
    """
    import matplotlib.pylab as plt
    from PyAstronomy import pyasl
    
    # Generate 200 exponentially correlated Gaussian
    # random numbers with a decay time of 5
    c1 = pyasl.expCorrRN(200, 5)
    
    # Generate 200 exponentially correlated Gaussian
    # random numbers with decay time 10, mean 4, and
    # standard deviation of 2.3.
    #
    # The results are: The correlated random numbers,
    # the uncorrelated numbers used as input, and the
    # correlated coefficient (exp(-1/tau)).
    c2, g, f = pyasl.expCorrRN(200, 10, mean=4.0, std=2.3, fullOut=True)
    
    plt.subplot(2,1,1)
    plt.plot(range(200), c1, 'bp-')
    plt.subplot(2,1,2)
    plt.plot(range(200), c2, 'bp-')
    plt.plot(range(200), g, 'g.')
#     plt.show()

  def sanity_expCorrRN_Example2(self):
    """
      Sanity of example 2 for expCorrRN
    """
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import pyasl
    
    # Generate n exponentially correlated Gaussian
    # random numbers with a decay time, tau
    n = 500
    tau = 5.
    c1 = pyasl.expCorrRN(n, tau)
    
    # Obtain autocorrelation function
    ac = np.correlate(c1, c1, mode="full")[n-1:]
    
    # Plot correlated random numbers and autocorrelation
    # function along with exponential model.
    x = np.arange(float(n))
    plt.subplot(2,1,1)
    plt.plot(x, c1, 'bp-')
    plt.subplot(2,1,2)
    plt.plot(x, ac, 'b.')
    plt.plot(x, np.exp(-x/tau)*ac.max(), 'r--')
  #   plt.show()

  def sanity_ramirez2005Example(self):
    """
      Check sanity of Ramirez/Melendez 2005 example
    """
    from PyAstronomy import pyasl
    
    # Create class instance
    r = pyasl.Ramirez2005()
    
    # Which color bands are available
    print("Available color bands: ", r.availableBands())
    
    # Convert B-V to effective temperature and back
    bv = 0.75
    feh = 0.0
    teff = r.colorToTeff("B-V", bv, feh)
    bv1 = r.teffToColor("B-V", teff, feh)
    # Watch out for differences between input bv and the output bv1
    print("B-V = ", bv, ", Teff = ", teff, ", bv1 = ", bv1, ", bv-bv1 = ", bv-bv1)

  def sanity_ramirez2005Example2(self):
    """
      Checking Ramirez, Ballesteros example
    """
    from PyAstronomy import pyasl
    
    b = pyasl.BallesterosBV_T()
    r = pyasl.Ramirez2005()
    
    # Convert B-V to effective temperature and back
    for bv in [0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45]:
      tr = r.colorToTeff("B-V", bv, 0.0)
      tb = b.bv2T(bv)
    
      print(("B-V [mag] = {3:4.2f} : Teff (R05) = {0:4.0f} K, " + \
              "Teff (B12) = {1:4.0f} K, dTeff = {2: 4.0f} K").format(tr, tb, tr - tb, bv))

  def sanity_ramirez2005(self):
    """
      Check sanity Ramirez/Melendez 2005 conversions
    """
    from PyAstronomy import pyasl
    import numpy as np
    from itertools import product
    
    feh = 0.0
    r = pyasl.Ramirez2005()
    # Compare with numbers given in Allen (astrophysical quantities, 4th edition)
    for teff, bv in zip([7000.,5950.,5310.,4410.],[0.35,0.58,0.74,1.15]):
      self.assertAlmostEqual(bv, r.teffToColor("B-V", teff, feh), delta=0.03, \
                             msg="B-V conversion (Teff " + str(teff) + ") out of sync with Allen")

    for teff, bv in zip([6500.,6000.,5000.],[0.286,0.36,0.535]):
      self.assertAlmostEqual(bv, r.teffToColor("b-y", teff, feh), delta=0.01, \
                             msg="B-V conversion (Teff " + str(teff) + ") out of sync with Allen")

    for teff, bv in zip([5000.,4500.,4000.],[0.433,0.51,0.735]):
      self.assertAlmostEqual(bv, r.teffToColor("R_C-I_C", teff, feh, 'g'), delta=0.04, \
                       msg="R_C-I_C conversion (Teff " + str(teff) + ") out of sync with Allen (giants)")
  
    for teff, bv in zip([5050.,4660.,4050.],[0.86,1.,1.5]):
      self.assertAlmostEqual(bv, r.teffToColor("B-V", teff, feh, 'g'), delta=0.1, \
                       msg="B-V conversion (Teff " + str(teff) + ") out of sync with Allen (giants)")
    
    # Check consistency
    bv = np.arange(0.4, 1.01, 0.05)
    feh = np.arange(-2., 0.3, 0.2)
    for x in product(bv, feh):
      teff = r.colorToTeff("B-V", x[0], x[1], 'ms')
      col = r.teffToColor("B-V", teff, x[1], 'ms')
      self.assertAlmostEqual(x[0], col, delta=1e-2, \
                             msg="Two-way color-temperature conversion is inconsistent. bv = %.4e, feh = %.4e" % x)
  
  def sanity_BallesterosExample(self):
    """
      Checking example for Ballesteros
    """
    from PyAstronomy import pyasl
    
    b = pyasl.BallesterosBV_T()
    
    bv = 0.65
    
    # Convert B-V into effective temperature
    teff = b.bv2T(0.65)
    print("B-V = {0:4.2f} mag -> Teff = {1:4.0f} K".format(bv, teff))
    
    # Convert effective temperature into B-V color
    teff = 4568.0
    bv = b.t2bv(teff)
    print("Teff = {0:4.0f} K -> B-V = {1:4.2f} mag".format(teff, bv))
  
  def sanity_gyroAgeBarnes(self):
    """
      Check Gyro age of Barnes 2007
    """
    # Eq 17, 18 (B-V, period, relative error)
    dat = [[0.5, 7, 20], [0.65, 12, 15], [1.0, 20, 13], [1.5, 30, 13]]
    for d in dat:
      a, e = pyasl.gyroAgeBarnes(d[1], d[0])
      self.assertAlmostEqual(e/a*100., d[2], delta=1.0, \
                             msg="Relative errors of Barnes 2007 cannot be reproduced.")
  
  def sanity_sanity_gyroAgeBarnesExample(self):
    """
      Example gyro age (Barnes 2007)
    """
    from PyAstronomy import pyasl
    
    # Parameters of the Sun (Barnes 2007, p 1174)
    bv = 0.642
    p = 26.09
    
    # Obtain solar age ...
    age = pyasl.gyroAgeBarnes(p, bv)
    # ... and print it
    print("Solar age: {0:4.2f} +/- {1:4.2f} Ga".format(*age))

  def sanity_chromoAgeRHK(self):
    """
      Checking sanity of chromospheric age (Donahue)
    """
    import numpy as np
    # R'HK vs log10(age) (estimated from Fig. 1 in Donahue 1998)
    dat = [[-4.25, 7.], [-4.6, 9.1], [-5.0, 9.7]]
    for d in dat:
      a = np.log10(pyasl.chromoAgeRHK(d[0])*1e9)
      print(a, d[1])
      self.assertAlmostEqual(a, d[1], delta=0.3,
                             msg="Cannot reproduce age estimate from Donahue fig 1.")
  
  def sanity_chromoAgeRHKExample(self):
    """
      Example chromospheric age (Donahue)
    """
    from PyAstronomy import pyasl
    
    # Approximate chromospheric age of the Sun
    print("Solar age: {0:4.2f} Ga".format(pyasl.chromoAgeRHK(-4.95)))
  
  def sanity_abundancepatternsExample(self):
    """
      Checking example of AbundancePatterns
    """
    from PyAstronomy import pyasl
    
    ap = pyasl.AbundancePatterns()
    
    print("Names of the available abundance patterns:")
    print(ap.availablePatterns())
    
    print()
    print("Get the Asplund et al. pattern (aspl) as a dictionary using")
    print("atomic number as a key:")
    print(ap.pattern("aspl", form="dict", key="number"))
    
    print()
    print("Get (relative) number abundance of oxygen using elemental symbol:")
    print(ap.abundance("O", pat="wilm"))
    print("or atomic number")
    print(ap.abundance(8, pat="wilm"))
  
  def sanity_abundancePatterns(self):
    """
      Checking abundance pattern
    """
    from PyAstronomy import pyasl

    ap = pyasl.AbundancePatterns()
    
    self.assertAlmostEqual(2.82e-07, ap.abundance("P", pat="feld"), delta=1e-12, \
                           msg="Abundance of P (feld) incorrect.")
    p = ap.pattern("wilm", form="dict", key="symbol")
    self.assertAlmostEqual(p["Ti"], 6.46e-08, delta=1e-12, \
                           msg="Ti abundance broken (wilm)")
    