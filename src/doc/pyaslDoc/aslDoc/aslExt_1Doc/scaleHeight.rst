Atmospheric scale height
=============================

.. currentmodule:: PyAstronomy.pyasl

The atmospheric scale height, H, characterizes the extent of the atmosphere.
It is defined as

.. math::
	
	H = \frac{k_B T}{\mu m_u g}
	
where :math:`k_B` is the Boltzmann constant, T is the atmospheric temperature, :math:`\mu` is the mean
molecular weight, :math:`m_u` is the unified atomic mass unit and g is the gravitational
acceleration. PyAstronomy provides a function accepting SI units (:py:func:`atmosphericScaleHeight`) and an alternative accepting
Earth or Jovian units :py:func:`atmosphericScaleHeight_MR`). Required conversion constants are adopted from PyA's
constants package.


Example: 
~~~~~~~~~~~~~~~~~

::

	from PyAstronomy import pyasl
	
	T, mu, g = 290, 28.97, 9.8
	she = pyasl.atmosphericScaleHeight(T, mu, g) 
	
	print("Earth")
	print(f"T, mu, g = {T} K, {mu}, {g} m/s**2")
	print(f"Scale height = {she:4.1f} [km]")
	
	
	T, mu, mp, rp = 165, 2.2, 1, 1
	shj = pyasl.atmosphericScaleHeight_MR(T, mu, mp, rp, "J")
	
	print("Jupiter")
	print(f"T, mu, mp, rp = {T} K, {mu}, {mp} [MJ], {rp} [RJ]")
	print(f"Scale height = {shj:4.1f} [km]")

	
API documentation 
------------------------------------
  
.. autofunction:: atmosphericScaleHeight
.. autofunction:: atmosphericScaleHeight_MR
  