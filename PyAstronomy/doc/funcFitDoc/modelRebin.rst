Using an overbinned model in calculations
=============================================

Calculating the model at more points than actually needed can be very useful, e.g., if
finite integration times are to be taken into account. The *turnIntoRebin* method
turns every "normal" model into one, which uses rebinning. See the tutorial for
an example of the application.

.. currentmodule:: PyAstronomy.funcFit
.. autofunction:: turnIntoRebin

The class object returned by the *turnIntoRebin* function is the received class object
with a modified *evaluate* method and some extra functionality; technically, this is achieved by
inheritance. The class object denoted by *_ModelRebinDocu* below inherits the functionality of the
incoming class object and is returned.

.. autoclass:: _ModelRebinDocu
   :members: