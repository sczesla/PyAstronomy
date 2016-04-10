The Singular Value Decomposition (SVD)
========================================

.. p23ready

PyA's `SVD` class provides the means to carry out a
singular value decomposition as is, e.g., useful to
recover line profiles.

.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: SVD
   :members:

Example: Delta functions and rotational broadening
--------------------------------------------------

::
  
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import pyasl
    
    # Get some "data"
    wvl = np.arange(5000., 5010., 0.02)
    template = np.ones(len(wvl))
    
    # There are two sharp lines in the template
    template[100] = 0.5
    template[115] = 0.3
    # Apply rotational broadening to the delta spectrum
    nflux = pyasl.rotBroad(wvl, template, 0.5, 23.45, edgeHandling="firstlast")
    
    # Carry out decomposition
    svd = pyasl.SVD()
    # Use 51 bins for the width of the broadening function.
    # Needs to be an odd number.
    svd.decompose(template, 51)
    # Obtain the broadening function needed to
    # recover "observed" spectrum. Note that the
    # edges (51/2 bins) will actually be neglected.
    b = svd.getBroadeningFunction(nflux)
    # Get the model, which results from the broadening
    # function and the template; obtain the indices
    # where it applies, too.
    m, mind = svd.getModel(b, modelIndices=True)
    
    # Plot the outcome
    plt.plot(b, 'bp-')
    plt.plot(mind, m, 'r.')
    plt.plot(nflux, 'g--')
    plt.show()


Example: Adding noise and neglecting singular values
------------------------------------------------------

::
    
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import pyasl
    
    # Get some "data"
    wvl = np.arange(5000., 5010., 0.02)
    template = np.ones(len(wvl))
    
    # There are two sharp lines in the template
    template[100] = 0.5
    template[115] = 0.3
    # Apply rotational broadening to the delta spectrum
    nflux = pyasl.rotBroad(wvl, template, 0.5, 23.45, edgeHandling="firstlast")
    nflux += np.random.normal(0., 0.005, len(nflux))
    
    # Carry out decomposition
    svd = pyasl.SVD()
    svd.decompose(template, 51)
    
    # Access the singular values
    sv = svd.getSingularValues()
    
    # Calculate the reduced chi square as a function of the number
    # of singular values neglected in the calculation of the
    # model.
    chi = []
    for i in range(1, len(sv), 5):
      b = svd.getBroadeningFunction(nflux, wlimit=sorted(sv)[i])
      m, mind = svd.getModel(b, modelIndices=True, asarray=True)
      chi.append( ((nflux[mind] - m)**2/0.005**2).sum() / len(mind) )
    
    plt.title("Reduced $\chi^2$ vs. number of neglected singular values")
    plt.plot(range(1, len(sv), 5), chi, 'bp-')
    plt.show()
