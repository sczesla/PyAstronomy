Elemental abundances
================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

Abundance patterns
-----------------------

.. autoclass:: AbundancePatterns
   :members:


Example of usage
~~~~~~~~~~~~~~~~~

::
    
    from __future__ import print_function, division
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
