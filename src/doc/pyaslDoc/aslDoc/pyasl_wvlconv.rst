Wavelength conversion (air and vacuum)
========================================

PyA offers the function to transform the wavelength axis:

 * :py:func:`airtovac <PyAstronomy.pyasl.airtovac>`
 * :py:func:`vactoair <PyAstronomy.pyasl.vactoair>`

Moreover, you can use:

 * :py:func:`specAirVacConvert <PyAstronomy.pyasl.specAirVacConvert>`
 
to convert from vacuum into air conditions or vice versa, but
retain the input wavelength axis.

`airtovac` and `vactoair`
---------------------------

Both functions have been ported from IDL's astrolib.

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: airtovac
.. autofunction:: vactoair


Example: From air to vacuum and back
-----------------------------------------

::

  from PyAstronomy import pyasl
  import numpy
  
  # Define wavelength array
  wvl = numpy.arange(10) + 5000.0
  print "Input wavelengths: ", wvl
  
  # Convert wavelength in air to wavelength
  # in vacuum
  wvlVac = pyasl.airtovac(wvl)
  print "Wavelength in vacuum: ", wvlVac
  
  # Convert wavelength from vacuum to air
  wvlAir = pyasl.vactoair(wvlVac)
  print "Wavelength in air: ", wvlAir

----------------------------------------------------

Transform spectrum, but retain wavelength axis
-----------------------------------------------------------

.. autofunction:: specAirVacConvert

Example: Convert spectrum from vacuum into air conditions
-----------------------------------------------------------

::

  from PyAstronomy import pyasl
  from PyAstronomy import funcFit as fuf
  import numpy as np
  import matplotlib.pylab as plt
  
  g = fuf.GaussFit1d()
  g["A"] = -0.1
  g["sig"] = 0.03
  g["mu"] = 5004.4752
  g["off"] = 1.0
  
  # Create a "spectrum" ...
  wvl = np.linspace(5000.,5010.,1000.)
  flux = g.evaluate(wvl)
  # ... and add some noise
  flux += np.random.normal(0., 0.02, len(wvl))
  
  # Assume that this spectrum is one observed
  # or modeled in vacuum conditions. We want
  # to convert it into a spectrum in air, but
  # retain the old wavelength axis.
  airflux, vind = pyasl.specAirVacConvert(wvl, flux, \
                  direction="vactoair")
  
  # Plot the result
  plt.plot(wvl, flux, "b.-")
  plt.plot(wvl, airflux, "r.-")
  plt.show()