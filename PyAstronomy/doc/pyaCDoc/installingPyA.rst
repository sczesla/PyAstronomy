Installation and source code of PyAstronomy
================================================

PyAstronomy is available from the `Python Package Index (PyPI) <https://pypi.org/project/PyAstronomy/>`_ using pip.

Sources, bug tracking, and opportunities for contributions are available on
`github <https://github.com/sczesla/PyAstronomy>`_.

Installation with PyPI and pip
---------------------------------

PyAstronomy can be installed via pip

.. note:: Depending on the setup of your Python installation, you may need administrator (root)
          privileges to install a package. 

::

    pip install PyAstronomy

or (if non-Python dependencies are required)

::    

    pip install PyAstronomy[occult]



Installation from github via pip
-----------------------------------

The current development branch can be installed from `github <https://github.com/sczesla/PyAstronomy>`_ via

::

    pip install git+https://github.com/sczesla/PyAstronomy.git#egg=PyAstronomy[occult]


Installing from source code
-------------------------------

PyAstronomy can be installed from the source.
Source distributions can be obtained from
`github <https://github.com/sczesla/PyAstronomy/releases>`_.
Save it to whatever place you prefer on your system, extract the files, and change into the thus created
directory; on linux use, e.g.,:

::
  
  tar xfv PyAstronomy.tar.gz
  cd PyAstronomy

.. note:: The package and directory name usually contain a version number.

In the directory created by unpacking the tar-ball, you find a script called *setup.py*.
This script will do the work of installing the package for you. Execute it by typing:

::
  
  python setup.py [--with-ext] install

.. note:: **--with-ext** is an optional flag. If specified, the installer will try to build
          non-Python extension. Building the extensions requires a fortran compiler. 

.. note:: Depending on the setup of your Python installation, you may need administrator (root)
          privileges to install a package. The default path for installing packages is the
          *site-packages* directory of your Python installation. You can modify this target location
          by using "python setup.py install --home=XYZ", where XYZ is your preferred installation
          path. Note that this path has to be added to your `PYTHONPATH` environment variable if
          XYZ is a nonstandard path.

Building the documentation
-----------------------------

PyAstronomy is distributed including documentation. The
latest documentation is available via `readthedocs <https://pyastronomy.readthedocs.io/en/latest/index.html>`_.
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


  
