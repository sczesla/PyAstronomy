PyA's "constants" package
==========================

.. currentmodule:: PyAstronomy.constants

There is probably no one who did not spend some time with
looking up constants for some calculation. PyA's `constants`
package has been introduced to help a little with this task.

This package uses the `quantities` package to manage units
and unit conversion. As most people are, however, interested
in having quick access to the numbers, the package focuses
on "easy access".

.. note:: This package requires the `quantities` package. 

Global vs. class scope
----------------------

The constants package provides the same functionality in a
*global* (module) scope and a class scope. The global scope has been
implemented to allow the easiest possible access. Using the global
scope within a function or class is, however, not advisable, because
it may cause or suffer from side effects. For this purpose, the
constants are also available from an object interface (the class scope),
which should be used instead. Both are demonstrated in the examples.

.. warning:: Using constants from global scope in a function/object may
             cause or suffer from side effects. Use class scope.


Examples
------------------

Below, we give an example of the usage of the `constants`
package. In this example, the module-scope constants are
used.

::
    
    from __future__ import print_function, division
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

 
::
      
    from __future__ import print_function, division
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


Custom constants
--------------------

In order to use a custom set of constants, you
need to set up a file holding them in appropriate format. This
package uses standard INI-style configuration files. An example
of a valid file can look like:

::

  [ArbitrarySectionName]
  descr   = Description of my custom constant
  symbol  = XXConst
  valueSI = 1.7656e58
  errSI   = 2e58
  unitSI  = W
  unitcgs = erg * s^-1
  source  = Nonsense for demonstration 
  
  [RFootball]
  descr   = Radius of a professional football
  symbol  = RFb
  valueSI = 0.11
  errSI   = 0.002
  unitSI  = m
  unitcgs = cm
  source  = Television

Then use :py:func:`load` with the name of your custom constants file. `load` will
loop through all sections defined in the file and read the definitions. The names
of the sections are arbitrary and will be ignored. Each section must contain the
following entries:
 
  - **descr**: A description of the constant.
  - **symbol**: The symbol (and attribute name) used to represent the constant. Must be unique.
  - **valueSI**: The value in SI units.
  - **errSI**: The error in SI units (if available, use 0 otherwise).
  - **unitSI**: The SI units (used by the `qunatities` package). E.g.: "m^2 * s^-1".
  - **unitcgs**: The unit in the cgs system.
  - **source**: From where did you get the numbers?

In case you do not want to predefined constants, you can use the :py:func:`cleanUp`
function defined in the package. This will delete all currently loaded constants and
remove the associated attributes from the package namespace.


Implementation details?
------------------------

Here we briefly describe how information is managed in this package.

On creation, the module loads a basic constants data-set from a default
file. The file is a simple Python configuration (INI format) file. The
format is described in the documentation of :py:func:`load`, which can
also be used to specify and load files holding information about
further constants.

All information about the constants is saved in the `inventory` attribute,
which is a dictionary mapping "constant symbol", i.e., the name of attribute
used to represent the constant, to a dictionary holding the details.

When a new file has been loaded, the constants now specified in the
`inventory` are mapped to module attributes to make them easily accessible.
The name of the attribute is the "symbol" defined for the constant. This
attribute will only contain the value of the constant. As a convention,
attributes with a prefix "f_" will store also the unit. The uncertainty
is mapped to an attribute of the for "f_SYMBOL_err".

The same rules apply to objects of type `PyAConstants`, which provide
the same functionality in the scope of a class.

The package API
------------------------

.. automodule:: PyAstronomy.constants
   :members:

The `PyAConstants` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: PyAstronomy.constants.PyAConstants
   :members: