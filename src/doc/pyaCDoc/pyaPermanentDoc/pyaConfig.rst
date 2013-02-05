PyA's root configuration
=========================

.. currentmodule:: PyAstronomy.pyaC.pyaPermanent
.. autoclass:: PyAConfig
   :members:

Example session
----------------

When PyAConfig is called without a configured data path,
it will ask for a location where PyA is allowed to store
its data and configuration. This will look something like
this:

::

  >>> from PyAstronomy.pyaC import pyaPermanent
  >>> pc = pyaPermanent.PyAConfig()
  Please provide a directory where PyA can store data (may already exist):
  Press enter to accept default.
    Path (default = /pathToHome/PyAData): another/Path/I/Like 
  PyA data path configured successfully. Using path: 
    another/Path/I/Like 
  
.. note:: The data path is **not** a subdirectory of the path you specify, but
          exactly the path you specify.