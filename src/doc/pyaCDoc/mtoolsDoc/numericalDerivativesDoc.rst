Compute numerical derivatives
===================================

.. currentmodule:: PyAstronomy.pyaC
.. autofunction:: diffCFD

Example
~~~~~~~~

::

  from PyAstronomy import pyaC
  import matplotlib.pylab as plt
  import numpy as np
  
  x = np.linspace(-10,10,1000)
  
  # Computer polynomial and its derivatives
  # (quasi analytically)
  y = [np.poly1d((0.03, -0.31, 0.4, 0.35, 1.4))]
  for i in range(4):
    y.append(y[-1].deriv())
  
  # Compute derivates numerically and compare to
  # analytic solution
  erro = 2
  for i in range(1,5):
    indi, der = pyaC.diffCFD(x, np.polyval(y[0], x), i, erro)
    plt.plot(x[indi], np.polyval(y[i], x[indi]), 'b.')
    plt.plot(x[indi], der, 'r--')
    plt.show()
