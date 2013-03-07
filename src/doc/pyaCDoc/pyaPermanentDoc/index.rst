PyA permanent --- Configuration and Data
=========================================

The *pyaPermanent* package is a core module of PyAstronomy. Its purpose is
to manage content, which shall be permanently available to PyA, but may
be created or added dynamically; examples comprise configuration and
newly downloaded data.

What gets stored where?
---------------------------

To store the data, PyA requires a directory (in the following the "data directory")
to which it has write access.
Such a directory can be assigned (and created) at the first call to a method,
which demands access to it.
By default, PyA will suggest to use a subdirectory of your home directory,
but this choice can be modified. If an appropriate directory
has been assigned to PyA, a
single file in your home directory will be created (named ".pyaConfigWhere"), which
contains a single line, which points PyA to its data directory; unfortunately, this
cannot be avoided, because PyA needs one point, which it can find without any additional
information.

In PyA's data directory, the file "pyaConfig.cfg" holds the "root configuration"
for PyA. Other modules may add own content, usually in subdirectories.

How to use it?
---------------

Of course, you can directly access the files in PyA's data path. However, PyA
provides some useful classes, which offer some additional convenience and maybe
"bug security". PyA packages will use these to access the data directory:

.. toctree::
   :maxdepth: 1
   
   pyaConfig.rst
   pyaFS.rst
   updateCycler.rst