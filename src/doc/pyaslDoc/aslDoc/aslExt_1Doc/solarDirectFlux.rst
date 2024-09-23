Solar direct flux
==================================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: solarDirectFluxMeinel


Example
~~~~~~~~

::

    from PyAstronomy import pyasl
    import numpy as np
    import matplotlib.pylab as plt
    
    z = np.linspace(0,89.9,200)
    
    for h in [0,2,4]:
    
        f = pyasl.solarDirectFluxMeinel(z, height=h)
    
        plt.plot(z, f, '-', label=f"Height = {h} km")
        
    plt.legend()
    plt.xlabel("Zenith angle [deg]")
    plt.ylabel("Direct solar flux [W/m**2]")
    plt.show()