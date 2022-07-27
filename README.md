# PyAstronomy

What is it?
-----------

  PyAstronomy is a collection of astronomy-related
  packages written in Python.

  Currently, the following subpackages are available:

    funcFit:    A convenient fitting package providing support
                for minimization and MCMC sampling.

    modelSuite: A Set of astrophysical models (e.g., transit
                light-curve modeling), which can be used
                stand-alone or with funcFit.

    AstroLib:   A set of useful routines including a number
                of ports from IDL's astrolib.

    Constants:  The package provides a number of often-needed
                constants.

    Timing:     Provides algorithms for timing analysis such as
                the Lomb-Scargle and the Generalized Lomb-Scargle
                periodogram

    pyaGUI:     A collection of GUI tools for interactive work.

Installation
------------

  To install the latest release via pip from PyPI use
  
    pip install PyAstronomy[occult]
    
  or
  
    pip install git+https://github.com/sczesla/PyAstronomy.git#egg=PyAstronomy[occult]
    
  to install the current state from github including non-Python dependencies. Remove [occult]
  to skip these dependencies.
  
  Alternatively, download the source and use
  
    python setup.py [--with-ext] install
    
  If you specify the `--with-ext` flag, setup will try to
  compile non-Python modules.

Documentation and further information
-------------------------------------

  View the latest documentation on [Read the
  Docs](https://pyastronomy.readthedocs.org/en/latest/)

  Visit the documentation of the latest release:
  
  https://pyastronomy.readthedocs.io/en/v_0-18-0/

Licensing
---------

  Where not stated otherwise, PyAstronomy is released under the
  MIT license (see also documentation).
