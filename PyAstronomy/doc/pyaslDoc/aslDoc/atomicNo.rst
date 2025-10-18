Atomic number and elemental symbols
========================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: AtomicNo
   :members:

Example
----------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    an = pyasl.AtomicNo()
    
    # Show table with all atomic numbers, elemental
    # symbols, and the names of the elements
    an.showAll()
    
    print(an.getElSymbol(26))
    
    print(an.getAtomicNo("He"))
    
    print(an.getElementName(25))
