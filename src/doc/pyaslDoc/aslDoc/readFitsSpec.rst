Read/write 1d spectrum from/to fits file
===========================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

The routines :py:func:`read1dFitsSpec` and :py:func:`write1dFitsSpec`
provide simple interfaces for reading and writing one-dimensional fits
spectra with pyfits (astropy.io.fits).

Reading fits spectrum
-------------------------

.. autofunction:: read1dFitsSpec

Writing fits spectrum
------------------------

.. autofunction:: write1dFitsSpec

Example of usage
~~~~~~~~~~~~~~~~~~~

::

    import numpy as np
    from PyAstronomy import pyasl
    
    # Generate a "spectrum"
    wvl = np.arange(5000., 5010., 0.01)
    flux = np.random.normal(1.0, 0.01, wvl.size)
    
    # Write spectrum providing wavelength array
    pyasl.write1dFitsSpec("test1.fits", flux, wvl=wvl, clobber=True)
    
    # Write spectrum specifying wavelength-related header keywords
    # manually
    wp = {"CRVAL1":5000., "CDELT1":0.01, "CRPIX1":1}
    pyasl.write1dFitsSpec("test2.fits", flux, waveParams=wp, clobber=True)
