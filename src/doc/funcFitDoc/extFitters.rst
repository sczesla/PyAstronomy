Internal and external optimizers
====================================

A "fitter" or "optimizer" is the algorithm used to carry out the actual
optimization. In funcFit, there are "internal" and "external" optimizers.
The internal optimizers are those invoked by a call to the fit-method
of the OneDFit object. External optimizers take a OneDFit object as
an object and carry out the fitting without actually becoming a part
of the object.

While, quite generally, internal optimizers are very convenient, it may
be a little easier to exert control over the screws determining the
details of the algorithm, if it needs to be tweaked to solve a specific
problem.  

Internal optimizers
-----------------------

Internal optimizers are invoked by OndeDFit's fit method. Basically,
an internal optimizer is a callable object, which carries out the
optimization when called. A base class for internal optimizers
is provided by funcFit.

External optimizers
---------------------

An external fitter is a fitting algorithm not implemented as a part
of funcFit's OneDFit class. Further algorithms can thus be supported,
which may give, e.g., stricter control over individual steps or
internal settings.  
 
A convenient way to pass data to the external fitters is the
"funcFit data set" :doc:`fufDS`.


List of internal and external optimizers
-------------------------------------------

.. toctree::
    :maxdepth: 3
    
    fitters.rst
