The stringlength method
===========================

The *stringlength* method is a technique to search for potentially nonsinusoidal, periodic
variation in a data set. The idea consists to fold the data with a number of trial periods.
Consecutive points in the phased data set are connected by hypothetical lines and the total
length (the string length) is calculated. The strong length statistic assumes a minimum, when
the folding produces a well ordered set. The stringlength method is described by in detail
by Dworetsky 1983, MNRAS 203, 917-924.


Example: A string length analysis
-------------------------------------

::
    
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import pyTiming as pyt
    
    period = 1.75
    
    x = np.linspace(0, 10, 100)
    y = 20 * np.sin(2*np.pi*x/period)
    y += np.random.normal(0, 1, len(x))
    
    # Trial periods to be tested (200 trial periods between 0.5 and 4.5;
    # same units as x-axis)
    tps = (0.5, 4.5, 200)
    
    # Calculate string length
    p, sl = pyt.stringlength_dat(x, y, tps)
    
    # Show the string length. An alias at the double period
    # (half frequency) is obvious.
    plt.plot(p, sl, 'b.-')
    plt.ylabel("String length")
    plt.xlabel("Trial period")
    plt.show()
    

Function documentation
---------------------------

.. currentmodule:: PyAstronomy.pyTiming
.. autofunction:: stringlength_dat
.. autofunction:: stringlength_pm
.. autofunction:: stringlength_norm
