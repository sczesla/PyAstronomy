Thermal broadening
====================

.. p23ready

.. currentModule:: PyAstronomy.pyasl
.. autofunction:: thermalBroadeningWidth
.. autofunction:: tempFromthermalBroadeningWidth

Example
--------

::

    from PyAstronomy import pyasl
    
    w0 = 6564.0
    T = 9567.0
    
    linefwhm = pyasl.thermalBroadeningWidth(w0, T)
    tbroad = pyasl.tempFromthermalBroadeningWidth(w0, linefwhm, awidth=0.0)
    
    print("Line width [FWHM]: %5.2f" % linefwhm)
    print("Thermal broadening temperature: %6.1f" % tbroad)
    
    
