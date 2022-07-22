Doppler shifting a spectrum
============================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: dopplerShift


Example: Shifting
----------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    import numpy as np
    
    # Create a "spectrum" with 0.01 A binning ...
    wvl = np.linspace(6000., 6100., 10000)
    # ... a gradient in the continuum ...
    flux = np.ones(len(wvl)) + (wvl/wvl.min())*0.05
    # ... and a Gaussian absorption line
    flux -= np.exp(-(wvl-6050.)**2/(2.*0.5**2))*0.05
    
    # Shift that spectrum redward by 20 km/s using
    # "firstlast" as edge handling method.
    nflux1, wlprime1 = pyasl.dopplerShift(wvl, flux, 20., edgeHandling="firstlast")
    
    # Shift the red-shifted spectrum blueward by 20 km/s, i.e.,
    # back on the initial spectrum.
    nflux2, wlprime = pyasl.dopplerShift(wvl, nflux1, -20.,
                                         edgeHandling="fillValue", fillValue=1.0)
    
    # Check the maximum difference in the central part
    indi = np.arange(len(flux)-200) + 100
    print("Maximal difference (without outer 100 bins): ",
          max(np.abs(flux[indi]-nflux2[indi])))
    
    # Plot the outcome
    plt.title("Initial (blue), shifted (red), and back-shifted (green) spectrum")
    plt.plot(wvl, flux, 'b.-')
    plt.plot(wvl, nflux1, 'r.-')
    plt.plot(wvl, nflux2, 'g.-')
    plt.show()
  
  
Example: Including uncertainties
------------------------------------
    
::
    
	from __future__ import print_function, division
	from PyAstronomy import pyasl
	import matplotlib.pylab as plt
	import numpy as np
	
	# Create a "spectrum" with 0.01 A binning ...
	wvl = np.linspace(6000., 6100., 1000)
	# ... and a Gaussian absorption line
	flux = 1 - 0.7 * np.exp(-(wvl-6050.)**2/(2.*1.5**2))
	# Add some noise
	err = np.ones_like(flux) * 0.05
	# Some points with unusually large error
	err[500] = 0.3
	err[[351, 497, 766, 787]] = 0.3
	flux += np.random.normal(0, err, len(flux))
	
	# Shift that spectrum to the blue by 17 km/s including errors
	nflux1, wlprime1, nerr1 = pyasl.dopplerShift(wvl, flux, -17., edgeHandling="firstlast", err=err)
	
	plt.errorbar(wvl, flux, yerr=err, fmt='b+', label="Original")
	plt.errorbar(wvl, nflux1+0.2, yerr=nerr1, fmt='r+', label="Shifted")
	plt.legend()
	plt.show()


  
Note on linear interpolation and noise
------------------------------------------

Say, we represent the spectrum by a series of independent random variables,
:math:`F_i`, normally distributed according to 

.. math::

	F_i \sim N(\mu_i,\sigma_i^2)
	
so that the :math:`\mu_i` are the mean values and the
:math:`\sigma_i` their standard deviations. A measurement of the spectrum with data points
:math:`f_1 \ldots f_n` is then a realization of these random variables. When we obtain a
shifted spectrum by linear interpolation between adjacent data points, :math:`f_i` and
:math:`f_{i+1}`, we obtain updated data points according to

.. math::

	g_i = a_i f_i + (1-a_i) f_{i+1} \sim N(a_i \mu_i + (1-a_i) \mu_{i+1}, a_i^2 \sigma_i^2 + (1-a_i)^2 \sigma_{i+1}^2 )
	
where the :math:`a_i` are the interpolation weightings so that :math:`0 \le a \le 1`. The sequence
:math:`g_i \ldots g_n` is a realization of respective random variables :math:`G_i`.
If the original spectrum is flat and zero (:math:`\mu_i=\mu=0`) and the noise is the same
throughout (:math:`\sigma_i=\sigma`), the usual unbiased estimator of the variance of the series
of data points

.. math::

	E[s^2(f_i \ldots f_n)] = E[F_i^2] = E\left[\frac{1}{N-1} \sum (f_i - \bar{f})^2\right] = \sigma^2
	
where :math:`\bar{f}` is the mean value has expectation :math:`\sigma^2`, which is here identical to the
variance of the individual variables.
The same
estimate for the shifted spectrum reads 

.. math::

	E[s^2(g_i \ldots g_n)] = E[G_i] = E\left[a F_i + (1-a) F_{i+1}\right] = (a^2 + (1-a)^2) \sigma^2 \le \sigma^2
	 
where we assume equal weights :math:`a_i=a`. The variance would thus be underestimated (using this estimator).
**Therefore, by shifting a flat spectrum by half a bin, its sample variance is cut in half.**
The explanation for this behavior is correlation between consecutive data points. The covariance
between adjacent points reads

.. math::

	COV(G_i, G_{i+1}) &= E[G_i G_{i+1}] = E\left[(a F_i + (1-a) F_{i+1}) (a F_{i+1} + (1-a)F_{i+2}) \right] = a(1-a) \sigma^2 \\
	COV(G_i, G_{i+j}) &= 0  \;\; \mbox{for} \;\; j > 1. 

Therefore, the autocorrelation function (ACF) of the resulting spectrum is non-zero only for lag one, where
one finds

.. math::

	\rho_1 = \frac{a(1-a)}{a^2 + (1-a)^2}





