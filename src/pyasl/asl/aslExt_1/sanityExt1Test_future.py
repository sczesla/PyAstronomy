from __future__ import print_function, division
import unittest

class SanityOfPyaslExt1(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_ExampleTransmissionCurves(self):
    """
      Checking example for 'TransmissionCurves'
    """
    #from __future__ import print_function
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
    #plt.show()
    
    
    # Create Planck spectrum ...
    wvl = np.arange(3000., 10000., 1.0)
    spec = pyasl.planck(T=5777., lam=wvl*1e-10)
    # ... and convolve with Johnson V band
    vbs = tcs.convolveWith(wvl, spec, "Johnson V")
    
    plt.plot(wvl, spec, 'b-', label='Input spectrum')
    plt.plot(wvl, vbs, 'r--', label='Convolution with Johnson V band')
    plt.legend()
    #plt.show()