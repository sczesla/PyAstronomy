SysRem
=================

.. p23ready
.. currentmodule:: PyAstronomy.pyasl

SysRem is an algorithm to remove systematic effects, defined as those which occur similarly (linearly) across a set of
observations such as sets of light curves or spectra. The algorithm
was described by
`Tamuz et al. 2005 (MNRAS 356, 1466) <https://ui.adsabs.harvard.edu/abs/2005MNRAS.356.1466T/abstract>`_
in the context of correcting systematic effects in samples of
light curves. Repeated application allows to remove more complex effects but the danger of removing actual
signal tends to increase with the number of iteration.

.. note::

    The data (or residuals) are stored such that X[::,0] gives the first observation, i.e.,
    the observation is stored in the first column. The data matrix is, therefore, stored here as
    the transpose of that adopted for the PCA etc..

.. note::
    
    Compared to the presentation by Tamuz et al., the roles of `a` and `c` are exchanged.


Example: Application to mock data
-------------------------------------

::

	from PyAstronomy import pyasl
	import numpy as np
	import matplotlib.pylab as plt
	
	# No. of data sets (e.g., light curves)
	nds = 50
	# No. of data points per data set
	nobs = 200
	
	obs, sigs = [], []
	
	# Generate mock data
	x = np.linspace(0,1,nobs)-0.5
	for i in range(nds):
	    # Add noise
	    y = np.random.normal(0,0.01+0.0001*i,nobs)
	    # Add 3rd degree polynomial
	    p = (i, 0.2*i, 0.3*i)
	    y += np.polyval(p,x)
	    # Add moving Gaussian signal
	    y -= 0.02 * np.exp(-(x-(i/nds-0.5))**2/(2*0.02**2))
	    # Add growing but stationary Gaussian signal
	    y -= (i**2*0.05) * np.exp(-x**2/(2*0.06**2))
	    obs.append(y)
	    sigs.append(np.ones_like(y)*0.01+0.0001*i)
	    
	    plt.plot(x, y, '.-')
	plt.title("Mock data")
	plt.show()
	
	sr = pyasl.SysRem(obs, sigs)
	# First iteration
	r, a, c = sr.iterate()
	plt.subplot(2,1,1)
	plt.title("Residuals after first (top) and second iteration (bottom)")
	plt.imshow(r, origin='lower', aspect="auto")
	# Second iteration
	r, a, c = sr.iterate()
	plt.subplot(2,1,2)
	plt.imshow(r, origin='lower', aspect="auto")
	plt.show()


Example: SYSREM, PCA, and eigenvectors
----------------------------------------

* :doc:`ex_sysrem_pca` :download:`(Download notebook) <ex_sysrem_pca.ipynb>`


API
----

Interface class
~~~~~~~~~~~~~~~~~

.. autoclass:: SysRem
   :members:

Helper functions
~~~~~~~~~~~~~~~~~~

.. autofunction:: sysrem_iter_c
.. autofunction:: sysrem_iter_a
.. autofunction:: sysrem_data_prepare
.. autofunction:: sysrem_update_rij
