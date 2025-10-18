Calculate the position angle
============================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: positionAngle

Example: Position angle of Alcor and Mizar
------------------------------------------------

::

    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Position of Mizar: 200.98141867 +54.92535197
    # Position of Alcor: 201.30640764 +54.98795966
    
    # Calculate position angle
    r = pyasl.positionAngle(200.98141867, +54.92535197, 201.30640764,+54.98795966)
    
    print("Position angle of Alcor (from Mizar): %4.2f deg" % r)
