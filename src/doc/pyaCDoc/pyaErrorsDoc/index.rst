PyAstronomy exceptions and warnings
=========================================

Occasionally, it is necessary to raise an exception or, at least, show a
warning. In the ideal case,
this should inform the user exactly about what went wrong (if
anything) and provide her/him with a solution to the problem.

PyA contains an *exception template*, which is the base class of all more
specialized PyA exceptions. The template defines a set of information, which
can be given to support the user to solve the problem. An exception class
may also be used to show a warning by handing it to the *warn* function.

Contents:

.. toctree::
   :maxdepth: 1
   
   pyaExTemplate.rst
   showWarning.rst
   pyaExceptions.rst
   