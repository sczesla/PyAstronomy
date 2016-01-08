List of internal and external optimizers
==========================================

.. currentmodule:: PyAstronomy.funcFit

Internal optimizers
-----------------------

Base class for internal optimizers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: IFitterBase
   :members:

ScipyFMIN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: ScipyFMIN
   :members:


FuFNM
~~~~~~~~~~~

FuFNM is based on the Nelder-Mead-Simplex implemenetation of funcFit (see :py:func:`NelderMead`).

.. autoclass:: FuFNM
   :members:
   


External optimizers
-------------------------

.. toctree::
    :maxdepth: 2
    
    nelderMead.rst


