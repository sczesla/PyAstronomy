Equidistant interpolation
============================

Interpolating tabulated data (x, y) onto an evenly sampled, equidistant axis
is a frequent problem. 

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: equidistantInterpolation


Example
----------

::

	import numpy as np
	import matplotlib.pylab as plt
	from PyAstronomy import pyasl
	
	x = np.array([0.1, 0.2, 0.5,0.87,1.5,2])
	y = np.array([1,2,3,1,4,1.1])
	
	w1, f1 = pyasl.equidistantInterpolation(x, y, "2x")
	w2, f2 = pyasl.equidistantInterpolation(x, y, "mean")
	w3, f3 = pyasl.equidistantInterpolation(x, y, 0.05)
	
	plt.plot(x,y,'bo', label="Data")
	plt.plot(w1, f1, 'r.-', label="2x")
	plt.plot(w2, f2, 'g.-', label="mean")
	plt.plot(w3, f3, 'm.-', label="0.02")
	plt.legend()
	plt.show()


Example (list of y axes and explicit x-axis)
-------------------------------------------------

::

	x = np.array([0.1, 0.2, 0.5,0.87,1.5,2])
	yy = [np.array([1,2,3,1,4,1.1]), np.array([10,20,30,1,10,-4.5])]
	
	# Apply interpolation to list of arrays for y
	w, ff = pyasl.equidistantInterpolation(x, yy, "2x")
	
	for f in ff:
	    plt.plot(w, f, '.-')
	
	# Specify new x-axis explicitly
	_, g = pyasl.equidistantInterpolation(x, yy[0], w)
	plt.plot(w, g, 'r--')
	plt.show()
