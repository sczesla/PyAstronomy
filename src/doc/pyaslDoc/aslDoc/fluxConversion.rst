Flux-photon conversion
=======================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: flux2photons
.. autofunction:: photons2flux

Example
--------

::

    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Wavelength in Angstrom
    wvl = 4000.
    # Flux in erg/s
    flux = 1.5e-14
    
    # Convert into photons
    photons = pyasl.flux2photons(wvl, flux)
    
    # How many photons is this?
    print("%g erg/s at %g A correspond to %g photons/s" \
            % (flux, wvl, photons))
    
    # Converting back
    flux2 = pyasl.photons2flux(wvl, photons)
    
    print("%g photons/s at %g A correspond to %g erg/s" \
            % (photons, wvl, flux2))