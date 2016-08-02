Rotational broadening
======================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

The functions here implement rotational broadening as described
in D.F. Gray's book "The Observation and Analysis of Stellar
Photospheres". Stellar limb-darkening is accounted for using
the linear law.

To apply rotational broadening, either :py:func:`rotBroad` or
:py:func:`fastRotBroad` can be used. The "fast" algorithm uses
a single broadening kernel, which is appropriate as long as
the wavelength range remains small. Otherwise the slower but
accurate algorithm implemented in :py:func:`rotBroad` should be
used.

.. autofunction:: rotBroad

Example
---------

Apply rotational broadening to a Gaussian and check the
result equivalent widths.

::

    from __future__ import print_function, division
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import funcFit as fuf
    from PyAstronomy import pyasl
    import scipy.integrate as sci
    
    # Create a spectrum with a single Gaussian
    # line using funcFit's GaussFit1d object.
    # Note that this object is not used for
    # fitting here, but only a calculate a
    # Gaussian.
    g = fuf.GaussFit1d()
    g["mu"] = 5005.
    g["A"] = -0.1
    g["sig"] = 0.1
    g["off"] = 1.0
    
    # Evaluate the spectrum with 0.01 A bin size
    wvl = np.linspace(5003., 5007., 400)
    flux = g.evaluate(wvl)
    
    # Obtain the broadened spectrum using
    # vsini = 13.3 km/s and no limb-darkening
    rflux = pyasl.rotBroad(wvl, flux, 0.0, 13.3)
    
    # Obtain the broadened spectrum using
    # vsini = 13.3 km/s and strong limb-darkening
    lflux = pyasl.rotBroad(wvl, flux, 0.9, 13.3)
    
    # Check that the area of the line did not change
    # in response to the broadening
    print("Initial EW [A]: ", 4. - sci.trapz(flux, wvl))
    print("After broadening without LD: ", 4. - sci.trapz(rflux, wvl))
    print("After broadening with LD: ", 4. - sci.trapz(lflux, wvl))
    
    # Plot the results
    plt.title("Rotational broadening")
    plt.xlabel("Wavelength [A]")
    plt.ylabel("Normalized flux")
    plt.plot(wvl, flux, 'b-')
    plt.plot(wvl, rflux, 'r-')
    plt.plot(wvl, lflux, 'g-')
    plt.show()


`fastRotBroad`---a faster algorithm
----------------------------------------


.. autofunction:: fastRotBroad

Example: Using the fast algorithm
----------------------------------

::
  
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import pyasl
    
    # Create data with a Gaussian absoprtion line
    wvl = np.arange(4999., 5011., 0.04)
    flux = np.zeros(len(wvl))
    
    # The Gaussian
    A = -0.05
    s = 0.1
    mu = 5004.1635788
    flux += A/np.sqrt(2.*np.pi*s**2) * \
      np.exp(-(wvl-mu)**2/(2.*s**2))
    
    # Apply the fast algorithm and ...
    bfast = pyasl.fastRotBroad(wvl, flux, 0.81, 11.37)
    # ... the slower one
    bslow = pyasl.rotBroad(wvl, flux, 0.81, 11.37)
    
    plt.xlabel("Wvl [A]")
    plt.ylabel("Flux [au]")
    plt.title("Initial spectrum (black), fast (blue), slow (red, shifted)")
    plt.plot(wvl, flux, 'k.-')
    plt.plot(wvl, bfast, 'b.-')
    plt.plot(wvl, bslow+0.01, 'r.-')
    plt.show()