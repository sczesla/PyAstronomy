Find zero crossings in discrete data set
==========================================

.. currentmodule:: PyAstronomy.pyaC
.. autofunction:: zerocross1d

Example
~~~~~~~~

::

    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import pyaC
    
    # Generate some 'data'
    x = np.arange(100.)**2
    y = np.sin(x)
    
    # Set the last data point to zero.
    # It will not be counted as a zero crossing!
    y[-1] = 0
    
    # Set point to zero. This will be counted as a
    # zero crossing
    y[10] = 0.0
    
    # Get coordinates and indices of zero crossings
    xc, xi = pyaC.zerocross1d(x, y, getIndices=True)
    
    # Plot the data
    plt.plot(x, y, 'b.-')
    # Add black points where the zero line is crossed
    plt.plot(xc, np.zeros(len(xc)), 'kp')
    # Add green points at data points preceding an actual
    # zero crossing.
    plt.plot(x[xi], y[xi], 'gp')
    plt.show()
