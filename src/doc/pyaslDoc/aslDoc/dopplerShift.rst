Doppler shifting a spectrum
============================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: dopplerShift

Note on linear interpolation and noise
------------------------------------------

Linear interpolation can (apparently) affect the noise in the data set.
Say we start with a number of iid (identically and independently
distributed) residuals, which obey a Gaussian distribution with some
fixed standard deviation characterizing the noise. This is what is
typically and often implicitly assumed in noise estimation.
Now we generate
a new data set by averaging adjacent points, i.e., we apply a special case
of linear interpolation with the new data points half way between the
initial ones.
If we now estimate the noise in the new data set under the
same iid assumption, we will find that it decreased by a factor
of sqrt(2). Now, it seems unlikely that we can increase the information content
in our data by interpolation alone. In fact, we must acknowledge that the iid
assumption does no longer hold, because the information previously contained in
a single data point is now, by construction, contained in two of the new, interpolated
data points, so that these points must no longer be considered independent.

While the effect is strongest for the case outlined above, it also holds for the
more general case, where interpolated data points are not half way between the
original ones. 

Example
--------

::

    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    import numpy as np
    
    # Create a "spectrum" with 0.01 A binning ...
    wvl = np.linspace(6000., 6100., 10000)
    # ... a gradient in the continuum ...
    flux = np.ones(len(wvl)) + (wvl/wvl.min())*0.05
    # ... and a Gaussian absoption line
    flux -= np.exp( -(wvl-6050.)**2/(2.*0.5**2) )*0.05
    
    # Shift that spectrum redward by 20 km/s using
    # "firstlast" as edge handling method.
    nflux1, wlprime1 = pyasl.dopplerShift(wvl, flux, 20., edgeHandling="firstlast")
    
    # Shift the red-shifted spectrum blueward by 20 km/s, i.e.,
    # back on the initial spectrum.
    nflux2, wlprime = pyasl.dopplerShift(wvl, nflux1, -20., \
                            edgeHandling="fillValue", fillValue=1.0)
    
    # Check the maximum difference in the central part
    indi = np.arange(len(flux)-200) + 100
    print("Maximal difference (without outer 100 bins): ", \
                    max(np.abs(flux[indi]-nflux2[indi])))
    
    # Plot the outcome
    plt.title("Initial (blue), shifted (red), and back-shifted (green) spectrum")
    plt.plot(wvl, flux, 'b.-')
    plt.plot(wvl, nflux1, 'r.-')
    plt.plot(wvl, nflux2, 'g.-')
    plt.show()
