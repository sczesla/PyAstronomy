Smoothing data
================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: smooth


Example of usage
-----------------

::
  
    from PyAstronomy import pyasl
    import numpy as np
    import matplotlib.pylab as plt
    
    # Create same "data" using a sine
    # and random noise
    x = np.linspace(0,10,100)
    y = 0.5*np.sin(x/3.0*2.0*np.pi + 1.7)
    y += np.random.normal(0.0, 0.2, len(y))
    
    # Use two smoothing windows with the same
    # window size
    sm1 = pyasl.smooth(y, 11, 'flat')
    sm2 = pyasl.smooth(y, 11, 'hamming')
    
    # Plot the outcome
    plt.title("Data and smoothed curves: flat (blue) and hamming window (red)")
    plt.plot(x, y, 'bp')
    plt.plot(x, sm1, 'b--')
    plt.plot(x, sm2, 'r--')
    plt.show()
