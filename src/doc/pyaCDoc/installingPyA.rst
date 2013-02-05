Installation of PyAstronomy
==============================

PyAstronomy is most easily obtained as a source package. The installation is supported by
the numpy's distutils.

Installing from source code
-------------------------------

PyAstronomy (PyA) uses numpy's distutils to support installation.
Start by downloading the source distribution and save it to whatever place
you prefer on your system. Extract the files and change into the thus created
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
          routines from Mandel and Agol. 

.. note:: Depending on the setup of your Python installation, you may need administrator (root)
          privileges to install a package. The default path for installing packages is the
          *site-packages* directory of your Python installation. You can modify this target location
          by using "python setup.py install --home=XYZ", where XYZ is your preferred installation
          path. Note that this path has to be added to your `PYTHONPATH` environment variable if
          XYZ is a nonstandard path.

Building the documentation
-----------------------------

PyAstronomy is distributed including documentation. To build the documentation yourself, change
into the directory where you installed (or unpacked) PyA. Change into the subdirectory named *doc*
(may not be the first level). In this directory, you find a Makefile, which is responsible for
building the documentation.

.. _Sphinx: http://sphinx.pocoo.org/

.. note:: To build the documentation you need to have installed the Sphinx_ package. In addition,
          PyAstronomy must be installed. 

The HTML documentation is built by using:

::
  
  make html


  