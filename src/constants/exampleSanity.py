from __future__ import print_function, division
import unittest

class ConstantsExampleSanity(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_example1(self):
    from PyAstronomy import constants as c
    
    # Print a summary of available constants
    # on screen
    c.summary()
    
    # Which unit system is in use?
    print()
    print("Current unit system: ", c.getSystem())
    
    # Access a constant
    print()
    print("Gravitational constant: ", c.G)
    # The 'f_' prefix is used as a convention. These
    # attributes hold `Quantity` objects as defined in
    # the `quantities` package. These encapsulate value
    # and unit. The prefixless attribute holds only the
    # number.
    print("             with unit: ", c.f_G)
    print("       error with unit: ", c.f_G_err)
    
    #Change the unit system
    print()
    print("Change the unit system")
    c.setSystem('SI')
    
    # Which unit system is in use?
    print("Current unit system: ", c.getSystem())
    
    # Access a constant again...
    print()
    print("Gravitational constant: ", c.G)
    print("             with unit: ", c.f_G)
    print("       error with unit: ", c.f_G_err)
    
    # Separate value and unit
    print()
    print("Value: ", c.f_G.magnitude, ", units: ", c.f_G.units)
    
    # Look up details
    print()
    print("What exactly was G?")
    c.constantDetails("G")
    
    # Apply unit conversion
    print()
    print("Use some other units")
    G_InFeet = c.inUnitsOf("G", "ft**3/(kg * s**2)")
    print("G with feet [ft**3/(kg * s**2)]: ", G_InFeet)

  def sanity_example2(self):
    from PyAstronomy.constants import PyAConstants
    
    c = PyAConstants()
    
    # Which unit system is in use?
    print()
    print("Current unit system: ", c.getSystem())
    
    # Access a constant
    print()
    print("Gravitational constant: ", c.G)
    # The 'f_' prefix is used as a convention. These
    # attributes hold `Quantity` objects as defined in
    # the `quantities` package. These encapsulate value
    # and unit. The prefixless attribute holds only the
    # number.
    print("             with unit: ", c.f_G)
    print("       error with unit: ", c.f_G_err)
