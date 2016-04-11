First ionization potential
============================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: FirstIonizationPot
   :members:

.. autofunction:: plotFIP

Example
--------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    fip = pyasl.FirstIonizationPot()
    
    print("First ionization energy of Li = %4.2e +/- %4.2e eV" % fip.getFIP(3))
    print("First ionization energy of Protactinium = %4.2e +/- %4.2e eV" % fip.getFIP(91))
    
    # And the same using elemental symbols
    print()
    print("First ionization energy of Li = %4.2e +/- %4.2e eV" % fip.getFIP("Li"))
    print("First ionization energy of Protactinium = %4.2e +/- %4.2e eV" % fip.getFIP("Pa"))
    
    # Plot the first ionization energy as a function of atomic number
    pyasl.plotFIP()
