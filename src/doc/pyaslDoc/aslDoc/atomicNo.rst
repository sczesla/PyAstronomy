Atomic number and elemental symbols
========================================

.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: AtomicNo
   :members:

Example
----------

::

  from PyAstronomy import pyasl
  
  an = pyasl.AtomicNo()
  
  # Show table with all atomic numbers, elemental
  # symbols, and the names of the elements
  an.showAll()
  
  print an.getElSymbol(26)
  
  print an.getAtomicNo("He")
  
  print an.getElementName(25)