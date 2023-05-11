Cross-correlation
===================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: crosscorrRV

.. seealso:: An algorithm for finding the extreme points by
             parabolic approximation (:py:func:`quadExtreme`).

Example: Cross-correlation with a Gaussian
-------------------------------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import numpy as np
    import matplotlib.pylab as plt
    
    # Create the template
    tw = np.linspace(5000, 5010, 1000)
    tf = np.exp(-(tw-5004.0)**2/(2.*0.1**2))
    
    # Create data, which are not that well sampled
    dw = np.linspace(5000, 5010, 200)
    df = np.exp(-(dw-5004.17)**2/(2.*0.1**2))
    
    # Plot template and data
    plt.title("Template (blue) and data (red)")
    plt.plot(tw, tf, 'b.-')
    plt.plot(dw, df, 'r.-')
    plt.show()
    
    # Carry out the cross-correlation.
    # The RV-range is -30 - +30 km/s in steps of 0.6 km/s.
    # The first and last 20 points of the data are skipped.
    rv, cc = pyasl.crosscorrRV(dw, df, tw, tf, -30., 30., 30./50., skipedge=20)
    
    # Find the index of maximum cross-correlation function
    maxind = np.argmax(cc)
    
    print("Cross-correlation function is maximized at dRV = ", rv[maxind], " km/s")
    if rv[maxind] > 0.0:
        print("  A red-shift with respect to the template")
    else:
        print("  A blue-shift with respect to the template")
    
    plt.plot(rv, cc, 'bp-')
    plt.plot(rv[maxind], cc[maxind], 'ro')
    plt.show()
    

Example: Cross-correlation including weights
------------------------------------------------
    
::
	
	from __future__ import print_function, division
	from PyAstronomy import pyasl
	import numpy as np
	import matplotlib.pylab as plt
	
	# Create the template
	tw = np.linspace(5000, 5010, 1000)
	tf = np.exp(-(tw-5004.0)**2/(2.*0.1**2))
	
	# Create data, which are not that well sampled
	dw = np.linspace(5000, 5010, 200)
	df = np.exp(-(dw-5004.17)**2/(2.*0.1**2))
	df += 0.2*np.exp(-(dw-5004.67)**2/(2.*0.1**2))
	# Define errors, use large errors for strong peak
	err = np.ones_like(dw) * 0.01
	err[(dw > 5004) & (dw < 5004.4)] *= 40
	
	# Plot template and data
	plt.subplot(2,1,1)
	plt.title("Template (blue) and data (red)")
	plt.plot(tw, tf, 'b.-', label="Template")
	plt.errorbar(dw, df, yerr=err, fmt='r+-', label="Date")
	plt.legend()
	
	# Carry out the cross-correlation.
	# The RV-range is -30 - +100 km/s in steps of 0.6 km/s.
	# The first and last 50 points of the data are skipped.
	rv, cc = pyasl.crosscorrRV(dw, df, tw, tf, -30., 100., 1., skipedge=50)
	
	# Include weights (1/err**2)
	rvw, ccw = pyasl.crosscorrRV(dw, df, tw, tf, -30., 100., 1., skipedge=50, \
	    weights=1/err**2)
	
	plt.subplot(2,1,2)
	plt.plot(rv, cc/max(cc), 'b-', label="normalized unweighted CCF")
	plt.plot(rvw, ccw/max(ccw), 'r.-', label="normalized weighted CCF")
	plt.legend()
	plt.show()

