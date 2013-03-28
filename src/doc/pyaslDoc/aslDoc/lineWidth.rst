Convert damping constant into line width
==========================================

Convert damping constant into natural line-width.

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: convertDampingConstant

Example
---------

::

  from PyAstronomy import pyasl
  
  # Einstein coefficient relevant for hydrogen LyA
  gLya = 6.258085e8
  
  print "Width of H LyA line at 1215.67 A = %e cm" % \
        pyasl.convertDampingConstant(gLya, 1215.67)