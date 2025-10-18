What makes a PyA exception? The PyA exception template
==============================================================

The PyA exception template is designed to provide information about
several issues, which, if the exception is not caught right away,
should help the user or developer to figure out the problem and,
hopefully, resolve it.

The questions to be addressed in the exception template are:

 - **What happened?**
 
   This point is the only mandatory field. The user should be provided with information \
   about what actually happened.
   
 - **Where did it happen?**
   
   The location where the exception occurred. For example, the module, class, and/or function.
 
 - **Why did it happen?**
 
   The reason or potential reasons for the exception to be raised. For example: The latitude may be \
   too large, because you are using degrees instead of radian. 
  
 - **How can it be solved?**
   
   If possible, provide a solution or a suggestion to solve the problem. For example: Set the \
   *deg* flag to True, if you use degrees.
 
 - **Additional information available?**
 
   For example: See web page abcde.com for further information on the algorithm.

Additionally, the template takes an argument specifying the *type* of the error to
give a rough classification of the exception. 

In many cases, not all of the fields can be addressed. The first point (what happened) \
is the only mandatory point; all others may be specified optionally.

The template implementation
----------------------------------------

.. currentModule:: PyAstronomy.pyaC.pyaErrors
.. autoclass:: PyaErrTemplate


Creating a custom exception class
-------------------------------------

All PyA exception classes should derive from the exception template class (*PyaErrTemplate*). \
The name of the exception raised is the first point where the user learns something \
about the nature of the error. Therefore, the classes should be named carefully.

A simple error class may look like this:

::
  
  class MyNewPyAError(PyaErrTemplate):
    
    def __init__(self, what, **keys):
      PyaErrTemplate.__init__(self, what, "My own error type", **keys)

The *keys* are then those given to the template, specifying it in more detail. Every exception should \
contain documentation, which must, at least, address the question: *When to be raised?*.