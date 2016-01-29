Numerical integration
==========================

.. currentmodule:: PyAstronomy.pyaC

Trapezoid rule with interpolated boundaries
------------------------------------------------

.. autofunction:: ibtrapz

Example
~~~~~~~~~

::

    from __future__ import print_function
    from PyAstronomy.pyaC import mtools
    import numpy as np
    
    x = np.arange(-2.,2.01,0.1)
    y = x**3 + 1.7
    
    x0 = -1.375
    x1 = +1.943
    
    # Analytical value of integral
    analyt = 0.25*(x1**4 - x0**4) + 1.7*(x1-x0)
    
    print("Analytical value: ", analyt)
    print("ibtrapz: ", mtools.ibtrapz(x, y, x0, x1))
