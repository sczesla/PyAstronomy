Estimate the Signal to Noise Ratio (SNR)
===========================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

The :py:func:`estimateSNR` function provides an algorithm
to estimate the SNR in a data set.

.. autofunction:: estimateSNR

Example---Estimating SNR ratio
--------------------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import numpy as np
    
    # Number of data points
    N = 10000
    # Signal to noise ratio
    SNR = 50.0
    
    # Create some data with noise and a sinusoidal
    # variation.
    x = np.arange(N)
    y = np.random.normal(0.0, 1.0/SNR, N) + 1.0
    y += np.sin(x/500.0*2.*np.pi)*0.1
    
    # Estimate the signal to noise ratio. Check whether the
    # estimate fits the input...
    # Use a chunk length of 20 data points, a polynomial of degree
    # one, and produce a "control plot".
    snrEsti = pyasl.estimateSNR(x, y, 20, deg=1, controlPlot=True)
    print("Estimate of the SNR: ", snrEsti["SNR-Estimate"])
    
    # Use a chunks with a length of 27, a polynomial of degree
    # two, and produce a "control plot".
    snrEsti = pyasl.estimateSNR(x, y, 27, deg=2, controlPlot=False, xlenMode="excerpt")
    print("Estimate of the SNR: ", snrEsti["SNR-Estimate"])
