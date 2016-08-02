S-Index and RHK
=================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: SMW_RHK
   :members:

Convert Mount-Wilson S-index into RHK
----------------------------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    ss = pyasl.SMW_RHK()
    
    bv = 0.8
    teff = 5100.0
    s = 0.4
    
    print("Convert S-index to RHK assuming a giant")
    ss.SMWtoRHK(s, teff, bv, lc="g", verbose=True)
    
    print()
    print()
    print("Convert S-index to RHK assuming a main-sequence star")
    ss.SMWtoRHK(s, teff, bv, lc="ms", verbose=True)


Show the Ccf conversion factor
-------------------------------

::

    from PyAstronomy import pyasl
    import numpy as np
    import matplotlib.pylab as plt
    
    ss = pyasl.SMW_RHK()
    
    bv = np.arange(0.4,0.9,0.05)
    ccfn = bv * 0.0
    ccfr = bv * 0.0
    ccfrg = bv * 0.0
    
    for i in range(len(bv)):
      ccfn[i] = ss.log10ccfNoyes(bv[i])
      ccfr[i] = ss.log10ccfRutten(bv[i])
      ccfrg[i] = ss.log10ccfRutten(bv[i], lc="g")
    
    plt.plot(bv, ccfn, 'b.-', label="Noyes")
    plt.plot(bv, ccfr, 'r.-', label="Rutten (ms)")
    plt.plot(bv, ccfrg, 'g.-', label="Rutten (g)")
    plt.xlabel("B - V [mag]")
    plt.ylabel("Ccf")
    plt.legend()
    plt.show()
