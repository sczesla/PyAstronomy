
Thin shell transit-model
===============================

.. p23ready

Transit light-curves usually have one minimum and a "U" shape. This can be
different, e.g., if optically thin chromospheric emission-lines are considered.
In this case, most emission may come from the
stellar limb resulting in a "W"-shaped transit-profile.

The model has been presented by Schlawin et al. [#fSchlawin]_.

The model class
----------------

.. currentmodule:: PyAstronomy.modelSuite.XTran.limBrightTrans
.. autoclass:: LimBrightTrans
   :members:




Example code - Calculate model light curve
------------------------------------------

::

  
    # Import some unrelated modules
    import matplotlib.pylab as mpl
    from numpy import arange
    # ... and now the LimBrightTrans module
    from PyAstronomy.modelSuite import LimBrightTrans
  
  
    # Create LimBrightTransit object
    lbt = LimBrightTrans()
    # Set parameters
    lbt["p"] = 0.08
    lbt["a"] = 6.70
    lbt["i"] = 87.84
    lbt["T0"] = 4.0
    lbt["per"] = 10.
  
    # Choose some time axis and calculate model
    time=arange(3.,5.,0.001)
    y = lbt.evaluate(time)
  
    # Let's see what happened...
    mpl.ylabel("Relative Flux")
    mpl.xlabel("Time")
    mpl.plot(time, y, 'r-')
    mpl.show()


.. [#fSchlawin] Schlawin et al. 2010, "Exoplanetary Transits of Limb-brightened Lines:
                Tentative Si IV Absorption by HD 209458b", 2010ApJ...722L..75S