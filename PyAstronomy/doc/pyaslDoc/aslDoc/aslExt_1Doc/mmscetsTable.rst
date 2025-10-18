A Modern Mean Stellar Color and Effective Temperature Sequence for O9V-Y0V Dwarf Stars
===========================================================================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

A table entitled "A Modern Mean Stellar Color and Effective Temperature Sequence for O9V-Y0V Dwarf Stars" is
kindly made available by `E. Mamajek <http://www.pas.rochester.edu/~emamajek/>`_. It holds a wealth of
information on various stellar parameters as a function of spectral type. 

.. autoclass:: MMSCETSTable
   :members:


Example
----------

::
    
    from __future__ import print_function
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    
    # Instantiate class
    m = pyasl.MMSCETSTable()
    
    # Print the entire data file
    for l in m.getContent():
        print(l, end='')
    
    print()
    print("Available columns: ", ", ".join(m.availableColumns()))
    
    # Get the tabulated data as an ASCII table
    td = m.getTable()
    
    # Plot absolute visual brightness vs. effective temperature
    plt.plot(td["Teff"], td["Mv"], 'b.-')
    plt.xlabel("Teff [K]")
    plt.ylabel("Mv [mag]")
    plt.show()
