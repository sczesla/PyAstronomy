Poisson source count-rate posterior with background
======================================================

.. currentmodule:: PyAstronomy.pyasl

.. autofunction:: ps_pdf_lams

Calculate factor Ci
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: PyAstronomy.pyasl.asl.aslExt_1.poissonPosterior._cirob


Example
-----------

::

    from __future__ import print_function, division
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import pyasl
    import scipy.integrate as scinteg
    
    # No. of counts in source and background region
    ns = 20
    nb = 100
    # Scaling factor between source and background region
    f = 0.1
    
    # Source rates at which to evaluate posterior
    lss = np.arange(0., 35., 0.1)
    # Posterior density
    post = pyasl.ps_pdf_lams(lss, ns, nb, f)
    
    # Check normalization of posterior
    norm = scinteg.trapz(post, lss)
    print("(Approximate) Normlization of posterior: ", norm)
    
    # Expectation value
    exv = scinteg.trapz(lss*post, lss)
    print("Posterior expectation value of lam_s: ", exv)
    
    plt.title("Posterior density")
    plt.plot(lss, post, 'b.-')
    plt.xlabel("$\\lambda_s$")
    plt.show()

    