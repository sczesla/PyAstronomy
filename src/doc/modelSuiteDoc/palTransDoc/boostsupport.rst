Adding Boost support
=====================

.. _Boost: http://www.boost.org/

The Boost_ libraries are a set of high-quality peer-reviewed c++ libraries written and
maintained by excellent programmers.

Boost offers both a convenient interface to create c++/c written Python modules and
an implementation of the elliptical integrals of the third kind. 

Boost support is mainly an issue of speed. A major part of the evaluation of the transit model
is evaluating elliptic integrals. This can conveniently be done, using the *mpmath* module, which
is written in pure Python. Yet, the evaluation is much faster, if the Boost library is applied.

How can Boost support be activated?
------------------------------------

To enable boost support, unfortunately, a little bit of manual work is needed.
This may require some knowledge about
compiling c++ code and the Boost libraries in particular. In the *palTrans* directory of the source
distribution of PyAstronomy, you will find the *ellint* subdirectory, which contains two files ell.cpp and
a makefile_template.

First you have to adapt the makefile_template to contain the library paths appropriate for your system;
then rename the makefile_template and call in "makefile".
Compiling the ell.cpp file using the makefile, produces a shared object library *ell.so*, if
successful. This library file has, finally, to be copied to the installation directory of
PyAstronomy.modelSuite.palTrans so that it can be imported by this module.

You can use the *whichEllInts()* method of the PalLC class to check whether Boost support is
present.