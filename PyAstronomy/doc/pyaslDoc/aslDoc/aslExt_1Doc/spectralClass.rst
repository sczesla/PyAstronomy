Spectral type vs. Teff and luminosity
==================================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

The class :py:func:`SpecTypeDeJager` implements the relation given by
de Jager and Nieuwenhuijzen 1987, A&A 177, 217-227.

Example: Basic Usage
--------------------------

::

    from __future__ import print_function
    from PyAstronomy import pyasl
    
    # Instantiate class object
    sdj = pyasl.SpecTypeDeJager()
    
    llum, lteff = sdj.lumAndTeff("K0", "V")
    
    print("Luminosity = {0:4.2f} Lsun".format(10.0**llum))
    print("Effective temperature = {0:6.1f} K".format(10.0**lteff))


Example: Teff and luminosity as a function of spectral type
--------------------------------------------------------------

::
    
    from __future__ import print_function
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    
    # Instantiate class object
    sdj = pyasl.SpecTypeDeJager()
    
    # Set luminosity class
    lk = "V"
    
    # Save spectral types, log(teff), and log(luminosity)
    spts = []
    lteffs = []
    llums = []
    
    # Save information to annotate abscissa
    xt = []
    xtl = []
    
    for t in "OBAFGKM":
        for n in range(10):
            if (t == "O") and (n == 0):
                # Skip the invalid "O0" type
                continue
    
            # Save the spectral type
            spts.append(t + str(n))
    
            # Get log10 of luminosity and effective temperature
            ll, lt = sdj.lumAndTeff(spts[-1], lk)
            # and save to lists
            llums.append(ll)
            lteffs.append(lt)
    
            # Save location (i.e., number in the list) and
            # spectral for annotating the abscissa
            if (n == 0) or (n == 5):
                xt.append(len(spts)-1)
                xtl.append(spts[-1])
    
    
    ax1 = plt.subplot(2, 1, 1)
    # Plot log10(effective temperature)
    plt.plot(lteffs)
    plt.ylabel("$\log_{10}$(T$_{eff}$)")
    plt.setp(ax1, xticks=xt, xticklabels=xtl)
    ax2 = plt.subplot(2, 1, 2)
    # Plot log10(luminosity)
    plt.plot(llums)
    plt.ylabel("$\log_{10}$($L/L_{\odot}$)")
    plt.setp(ax2, xticks=xt, xticklabels=xtl)
    plt.xlabel("Spectral type")
    plt.show()

API
-------

.. autoclass:: SpecTypeDeJager
   :members:
   :private-members:
   
