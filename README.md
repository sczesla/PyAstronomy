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

  Use `python setup.py [--with-ext] install`

  If you specify the `--with-ext` flag, setup will try to
  compile non-Python modules. This is not mandatory.

Documentation and further information
-------------------------------------

  View the latest documentation on [Read the
  Docs](https://pyastronomy.readthedocs.org/en/latest/)

  Visit the documentation of the latest release:
  
  http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/index.html

Licensing
---------

  Where not stated otherwise, PyAstronomy is released under the
  MIT license (see also documentation).
