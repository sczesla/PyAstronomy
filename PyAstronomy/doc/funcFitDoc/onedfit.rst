Fitting object base class
==========================

The *OneDFit* class provides a convenient interface to fitting algorithms.

.. currentmodule:: PyAstronomy.funcFit
.. autoclass:: OneDFit
   :members:
   :inherited-members:


Variable naming
------------------------

.. autoclass:: ModelNameIdentBase
   :members:


Python decorators and *MiniFunc*
-------------------------------------------

In the following, we give a brief introduction to Python decorators and outline the use of the \
MiniFunc decorator.

:Decorators in Python:
  Loosely speaking, decorators can be used to add or modify the functionality of a callable object \
  such as a function. This can for example be useful if a whole set of function shares common functionality or \
  if a certain fraction of the required functionality is to be ``hidden'' from the user, as is the case here.
  
  Applying a decorator, which is in Python marked by the ``@'' sign, turns one callable object into another \
  one. In our case, it will be necessary to have a decorator with arguments. A minimalistic example \
  may look like this:

::

  class Deco:
    def __init__(self, argument):
      self.argument = argument
    def __call__(self, f):
      print "Deco.__call__ is speaking!"
      def newFunc(y):
        print "First argument was: ", self.argument
        print "But now I received: ", y
        print "Evaluating function at given point: "
        f(y)
      return newFunc
  
  
  def f1(x):
    print "x**2 = ", x**2
  
  @Deco(13)
  def f2(x):
    print "x**2 = ", x**2
  
  print "-------------------------------"
  
  f1(2)
  print "-------------------------------"
  f2(2)

Execution the above script results in the following output:

::

  Deco.__call__ is speaking!
  -------------------------------
  x**2 =  4
  -------------------------------
  First argument was:  13
  But now I received:  2
  Evaluating function at given point: 
  x**2 =  4

What happens?
At first sight, the functions *f1* and *f2* look exactly alike, yet the definition of *f2* is \
preceded by the *@Deco(13)* directive; obviously, a call to *f1* provides other output than a \
call of *f2*.

Above we defined a class called *Deco*. Instances of that class are callable objects, as they have
a *__call__* method. The *@Deco(13)* directive preceding *f2* triggers the following process:

1)  The constructor of Deco is called with a single argument (13),
2)  The __call__ method of Deco is called with a single argument,
    which is a function pointer pointing at *f2*.
3)  The __call__ method defines a new function, incorporating the functionality provided by *f2*,
    and returns it. This new function can combine the information given to the constructor and
    the arguments handed to the new function.

The whole process is, thus, replacing one callable object by another, possibly combining and extending the
functionality.

:MiniFunc:

The MiniFunc decorator can be used to apply user-defined minimization functions. In principle, such a \
function determines the quantity to be minimized (e.g., Chi square) for a given model and parameter vector. \
The implementation of the module requires to also carry out some other steps, which is done through the decorator, \
without making this visible for the user.

An example of a user defined objective function looks like:

::

  gf = GaussFit1d()
  
  @MiniObj(gf)
  def mini(odf, P):
    return sum(abs((odf.y - odf.model)/odf.yerr))

The first argument of the objective function is the fitting object, and the second is the parameter vector. \
Properties of the fitting objects such as `x`, `y`, or `model` may be accessed as usual. The return value is \
the quantity to be minimized.

.. note:: An example of the application of the decorator is given in the tutorial.

.. autoclass:: MiniFunc
   :members:


   