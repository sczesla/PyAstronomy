Rotational broadening profile
================================

.. p23ready

.. currentmodule:: PyAstronomy.modelSuite
.. autoclass:: RotBroadProfile
   :members:
  
Example of usage
-----------------
   
::
     
    import numpy as np
    import matplotlib.pylab as plt
    from PyAstronomy import modelSuite as ms
    
    # Get an instance of the model ...
    x = ms.RotBroadProfile()
    # ... and define some starting value
    x["xmax"] = 60.0
    x["A"] = 1.0
    x["eps"] = 0.8
    x["off"] = 0.0
    
    # Define a radial velocity axis
    vv = np.linspace(-90.,90.,200)
    
    # Construct some "data" and ...
    data = x.evaluate(vv)
    # ... add noise
    data += np.random.normal(0.0, 1e-3, data.size)
    
    # Fit the model using A, xmax, and eps as free
    # parameters ...
    x.thaw(["A", "xmax", "eps"])
    x.fit(vv, data)
    # ... and show the resulting parameter values.
    x.parameterSummary()
    
    # Plot the data and the model
    plt.plot(vv, data, 'bp')
    plt.plot(vv, x.model, 'r--')
    plt.show()