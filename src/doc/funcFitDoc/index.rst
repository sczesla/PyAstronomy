funcFit - A convenient fitting interface
=========================================

The funcFit package provides a convenient interface to the fitting algorithms
provided by the popular SciPy and pymc packages.
It implements a very flexible and simple parameter handling mechanism
making fitting in Python a much more enjoyable experience.

.. seealso:: PyAstronomy's *modelSuite*. The funcFit package itself comes with only a few basic fitting models such as
          a Gaussian. More complex models are provided in the frame of the *model suite*.

.. warning:: funcFit relies on **Python 2.7.x**. The 2.6.x series (and prior) has a bug affecting the copying of
             dynamically created class methods, which has not been (and will not be) corrected. This interferes with many of
             funcFit's algorithms.

The tutorial
-----------------

The funcFit tutorial gives you an introduction to the capabilities of the
package. The most important aspects are demonstrated by example. 

.. toctree::
    :maxdepth: 2

    tutorial.rst
    tutorial2.rst
    traceAnalysisTutorial.rst

External fitters
---------------------

.. toctree::
    :maxdepth: 2

    nelderMead.rst

Technical documentation
--------------------------

A more comprehensive (and technical) documentation of individual classes
and concepts can be found below:

.. toctree::
    :maxdepth: 2
    
    params.rst
    onedfit.rst
    modelRebin.rst
    simplemodels.rst
    syncFit.rst
    traceAnalysis.rst
    coordinateGrid.rst



