Parallactic angle
===================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: parallacticAngle

Example: Parallactic angle for objects observed by Keck
-----------------------------------------------------------

::

	import numpy as np
	import matplotlib.pylab as plt
	from PyAstronomy import pyasl
	
	# Geo-latitude of Keck [deg]
	gl = 19.83
	
	# Declinations to be considered
	decs = np.arange(-80, 81, 20)
	# Hour angles to be considered
	t = np.linspace(0.0, 8*15, 1000)
	
	for dec in decs:
	    q = pyasl.parallacticAngle(t, gl, dec)
	    plt.plot(t/15, q, '-', label="$\delta = $" + f"{dec} deg")
	
	plt.ylabel("Parallactic angle [deg]")
	plt.xlabel("Hour angle [h]")
	plt.legend()
	plt.show()

	