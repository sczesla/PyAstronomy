Installation of PyAstronomy
==============================

PyAstronomy is most easily obtained from the `PYthon Package Index (PYPI) <https://pypi.org/project/PyAstronomy/>`_ using pip.
Alternatively, it can be obtained as a source package.

Installation using PYPI and pip
---------------------------------

The easiest way to install PyAstronomy is by using pip (numpy may be left out; if already
installed the command triggers, if anything, an update of numpy).

.. note:: Depending on the setup of your Python installation, you may need administrator (root)
          privileges to install a package. 

::

    pip install numpy PyAstronomy


It is advisable to install a few other packages to work more efficiently with PyAstronomy
(and python). All these packages can also be obtained from PYPI. They can be installed using
the command

::

    
    pip install scipy matplotlib quantities emcee


Note that any package may be left out. The overhead of having more package installed is, however, rather
small.

Occasionally, PyAstronomy (and other packages) are updated. To obtain the latest version use

::

    pip install --upgrade PyAstronomy

Installing from source code
-------------------------------

PyAstronomy (PyA) uses numpy's distutils to support installation.
Start by downloading the source distribution, which is most easily obtained from
`github <https://github.com/sczesla/PyAstronomy>`_.
Save it to whatever place
you prefer on your system, extract the files, and change into the thus created
directory; on linux use:

::
  
  tar xfv PyAstronomy.tar.gz
  cd PyAstronomy

.. note:: The package and directory name usually contain a version number.

In the directory created by unpacking the tar-ball, you find a script called *setup.py*.
This script will do the work of installing the package for you. Execute it by typing:

::
  
  python setup.py [--with-ext] install

.. note:: **--with-ext** is an optional flag. If specified, the installer will try to build
          non-Python extension. As of version 0.4.0, this only comprises the transit modeling
          routines from Mandel and Agol.  Building the extensions requires a fortran compiler. 

.. note:: Depending on the setup of your Python installation, you may need administrator (root)
          privileges to install a package. The default path for installing packages is the
          *site-packages* directory of your Python installation. You can modify this target location
          by using "python setup.py install --home=XYZ", where XYZ is your preferred installation
          path. Note that this path has to be added to your `PYTHONPATH` environment variable if
          XYZ is a nonstandard path.

Building the documentation
-----------------------------

PyAstronomy is distributed including documentation. Note that the most convenient way to access the
latest documentation is via `readthedocs <https://pyastronomy.readthedocs.io/en/latest/index.html>`_.
To build the documentation yourself, change
into the directory where you installed (or unpacked) PyA. Change into the subdirectory named *doc*
(may not be the first level). In this directory, you find a Makefile, which is responsible for
building the documentation.

.. _Sphinx: http://sphinx.pocoo.org/

.. note:: To build the documentation you need to have installed the Sphinx_ package. In addition,
          PyAstronomy must be installed. 

The HTML documentation is built by using:

::
  
  make html


  
