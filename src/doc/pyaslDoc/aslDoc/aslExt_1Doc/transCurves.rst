Photometric filter transmission curves
=========================================

.. currentmodule:: PyAstronomy.pyasl

Broad-band photometry is an important tool to study variability and
band transmission curves are, e.g., required to derive photometric fluxes
from spectra. The :py:func:`TransmissionCurves` class provides access to
photometric transmission curves.

Further information on photometric bands:

  ==================  ===========================================================================================
  Band                Destination
  ==================  ===========================================================================================
  Bessel              Bessel 1990, PASP 102, 1181B (`helpful link
                      <http://spiff.rit.edu/classes/phys440/lectures/filters/filters.html>`_)
                      
  Johnson             Johnson and Morgan 1951, ApJ 114, 522 (`helpful link
                      <http://obswww.unige.ch/gcpd/filters/fil01.html>`_)
                      
  Spitzer IRAC        `Transmission curves for Spitzer IRAC instrument
                      <http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/spectralresponse/>`_
                      
                      Call **addSpitzerIRACPassbands** to add to passband inventory.                     
  ==================  ===========================================================================================


Example of usage
-----------------

::

    from __future__ import print_function
    from PyAstronomy import pyasl
    import numpy as np
    import matplotlib.pylab as plt
    # Import six for Python 2/3 compatibility
    import six
    
    # Get transmission curve object
    tcs = pyasl.TransmissionCurves()
    
    # Add passbands from Spitzer IRAC
    tcs.addSpitzerIRACPassbands()
    
    print("Available bands: ", tcs.availableBands())
    
    # Wavelength axis
    wvl = np.linspace(3000, 10000, 10000)
    
    # Plot transmission curves for Bessel b, v, and r bands
    for (b, c) in six.iteritems({"b":"b", "v":"k", "r":"r"}):
      tc = tcs.getTransCurve("Bessel " + b)
      trans = tc(wvl)
      plt.plot(wvl, trans, c+'-', label="Bessel " + b)
    
    # Plot transmission curves for Johnson U, B, and V bands
    for (b, c) in six.iteritems({"U":"m", "B":"b", "V":"k"}):
      tc = tcs.getTransCurve("Johnson " + b)
      trans = tc(wvl)
      plt.plot(wvl, trans, c+'--', label="Johnson " + b)
    
    plt.legend()
    plt.xlabel("Wavelength [$\AA$]")
    plt.ylabel("Transmission")
    plt.show()
    
    
    # Create Planck spectrum ...
    wvl = np.arange(3000., 10000., 1.0)
    spec = pyasl.planck(T=5777., lam=wvl*1e-10)
    # ... and convolve with Johnson V band
    vbs = tcs.convolveWith(wvl, spec, "Johnson V")
    
    plt.plot(wvl, spec, 'b-', label='Input spectrum')
    plt.plot(wvl, vbs, 'r--', label='Convolution with Johnson V band')
    plt.legend()
    plt.show()


API
---------

.. autoclass:: TransmissionCurves
   :members:

