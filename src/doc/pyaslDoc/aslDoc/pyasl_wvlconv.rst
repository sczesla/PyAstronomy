Wavelength conversion (air and vacuum)
========================================

.. p23ready

PyA provides tools to convert wavelengths from air into vacuum conditions and vice versa:

 * :py:func:`airtovac2 <PyAstronomy.pyasl.airtovac2>`
 * :py:func:`vactoair2 <PyAstronomy.pyasl.vactoair2>`

These functions allow to use the conversions specified by:
  - Edlen 1953, J. Opt. Soc. Am 43 no. 5
  - Peck and Reeder 1972, J. Opt. Soc. 62 no. 8
  - Ciddor 1996, Applied Optics 35 no. 9

.. warning::
   The function `airtovac` and `vactoair` are based on the formula by Edlen 1953.
   It seems that air wavelengths (or wavenumbers) are used, where vacuum wavenumbers
   should be used.

   * :py:func:`airtovac <PyAstronomy.pyasl.airtovac>`
   * :py:func:`vactoair <PyAstronomy.pyasl.vactoair>`

Moreover, you can use:

 * :py:func:`specAirVacConvert <PyAstronomy.pyasl.specAirVacConvert>`
 
to convert from vacuum into air conditions or vice versa, but
retain the input wavelength axis.

Usage examples for `airtovac2` and `vactoair2`
----------------------------------------------------

Example: From air to vacuum and back
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import numpy as np
    
    # Define wavelength array
    wvl = np.arange(10) + 5000.0
    print("Input wavelengths: ", wvl)
    
    # Convert wavelength in air to wavelength
    # in vacuum. By default, the conversion
    # specified by Ciddor 1996 are used.
    wvlVac = pyasl.airtovac2(wvl)
    print("Wavelength in vacuum: ", wvlVac)
    
    # Convert wavelength from vacuum to air
    wvlAir = pyasl.vactoair2(wvlVac)
    print("Wavelength in air: ", wvlAir)


Example: Compare Edlen and Ciddor conversions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::
  
    from PyAstronomy import pyasl
    import numpy
    import matplotlib.pylab as plt
    
    # Compare the Edlen and Ciddor conversions
    wvl = numpy.arange(3000, 10000, 1.0)
    wvlVacCiddor = pyasl.airtovac2(wvl, mode="ciddor")
    wvlVacEdlen = pyasl.airtovac2(wvl, mode="edlen53")
    
    plt.subplot(2,1,1)
    plt.title("Difference in air wavelength (Ciddor-Edlen)")
    plt.ylabel("dWvl [A]")
    plt.xlabel("Vacuum wvl [A]")
    plt.plot(wvl, wvlVacCiddor-wvlVacEdlen, 'b.-')
    plt.subplot(2,1,2)
    plt.title("Difference in air wavelength (Ciddor-Edlen, in RV)")
    plt.ylabel("dRV [m/s]")
    plt.xlabel("Vacuum wvl [A]")
    plt.plot(wvl, (wvlVacCiddor-wvlVacEdlen)/wvlVacCiddor*299792458., 'b.-')
    plt.show()


API documentation: `airtovac2` and `vactoair2`
----------------------------------------------------

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: airtovac2
.. autofunction:: vactoair2


Obtaining the refractive index
----------------------------------------------------

Although the conversion of wavelength between air and vacuum conditions
is the primary application targeted here, what is really needed to carry
out the conversion is the refractive index.

.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: RefractiveIndexAV
   :members:
   :private-members:
   
Example: Obtain refractive index of standard air
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import numpy as np
    
    # Define wavelength array
    wvl = np.arange(10) + 5000.0
    
    # Obtain refractive index according to Ciddor 1996
    ri = pyasl.RefractiveIndexAV(mode="ciddor")
    n = ri.refractiveIndex(wvl)
    
    print("Wavelength and 1.0 - Refractive index of 'standard air':")
    for w, nc in zip(wvl, n):
      print("{0:5.1f}  {1:10.8e}".format(w, nc-1.0))


Transform spectrum, but retain wavelength axis
-----------------------------------------------------------

.. autofunction:: specAirVacConvert

Example: Convert spectrum from vacuum into air conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

------------------------------------------------  

API documentation: `airtovac` and `vactoair`
----------------------------------------------------

.. note:: These function are now considered deprecated.

Both functions have been ported from IDL's astrolib.

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: airtovac
.. autofunction:: vactoair