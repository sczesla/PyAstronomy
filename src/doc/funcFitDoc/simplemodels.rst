Common fitting models
=====================
.. currentmodule:: PyAstronomy.funcFit

Several often needed models are shipped with funcFit.

:Common fitting model available:
  - Gaussian profile (class :py:class:`GaussFit1d`)
  - Multicomponent Gaussian (class :py:class:`MultiGauss1d`)
  - Cauchy-Lorentz profile (class :py:class:`CauchyLorentz1d`)
  - Voigt profile (class :py:class:`Voigt1d`)
  - Multicomponent Voigt (class :py:class:`MultiVoigt1d`)
  - Sine wave (class :py:class:`SinusFit1d`)
  - Exponential decay (class :py:class:`ExpDecayFit1d`)
  - Polynomial of degree *n* (class :py:class:`PolyFit1d`)
  - Constant (class :py:class:`ConstantFit1d`)
  - 2d Gaussian (class :py:class:`GaussFit2d`)
  - Multicomponent 2d Gaussian (class :py:class:`MultiGauss2d`)

The use of these models is demonstrated in the tutorial.

.. _GaussFit1d:

Single Gaussian with linear continuum
-------------------------------------
.. autoclass:: GaussFit1d

Multicomponent Gaussian with linear continuum
---------------------------------------------
.. autoclass:: MultiGauss1d

.. _CauchyLorentz1d:

Cauchy-Lorentz with linear continuum
---------------------------------------------
.. autoclass:: CauchyLorentz1d

.. _Voigt1d:

Voigt with linear continuum
---------------------------------------------
.. autoclass:: Voigt1d

.. _MultiVoigt1d

Multicomponent Voigt profile with linear continuum
---------------------------------------------------

.. autoclass:: MultiVoigt1d
   :members:

Sine fit
---------------------------------------------
.. autoclass:: SinusFit1d

Exponential decay
---------------------------------------------
.. autoclass:: ExpDecayFit1d

Polynomial of degree n
---------------------------------------------
.. autoclass:: PolyFit1d

Constant
---------------------------------------------
.. autoclass:: ConstantFit1d

Two-dimensional Gaussian
---------------------------------------------
.. autoclass:: GaussFit2d

Multicomponent two-dimensional Gaussian
---------------------------------------------
.. autoclass:: MultiGauss2d

Demonstration models
=====================

Models mainly used for demonstration purposes. 

Simple circular orbit
-----------------------

.. autoclass:: Circle2d