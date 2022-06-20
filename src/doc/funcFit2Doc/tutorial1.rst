The funcFit2 introduction and tutorial
======================================

This tutorial is supposed to enable you to exploit the funcFit2 functionality
by presenting an anthology of examples.

Some knowledge of numpy, SciPy, and matplotlib is helpful to dive in.

.. contents:: Sections of this tutorial

.. _matplotlib: http://matplotlib.sourceforge.net/
.. _pymc: https://github.com/pymc-devs/pymc
.. _SciPy: www.scipy.org/
.. _numpy: numpy.scipy.org/
.. _XSPEC: http://heasarc.nasa.gov/xanadu/xspec/
.. _emcee: http://dan.iel.fm/emcee/current/

.. note:: Please mind the citation request for use of the scipy algorithms explained
          `here <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html>`_.
  

Prerequisites
-------------
To run the example in this tutorial you need to have installed the following packages:
 * numpy_
 * SciPy_
 * matplotlib_

After installing PyAstronomy (PyA), funcFit2
is available for being used. 
As a first step, let us import the
package and check the status of scipy and emcee in out installation.
To that end,

::

  from PyAstronomy import funcFit2 as fuf2
  
  fuf2.status()

Depending on your installation the output should look like:

::

    Status of funcFit2:
    --------------------------
    Is scipy.optimize available?  yes
      Version:  1.1.0
    Is emcee available?  yes
      Version:  2.2.1


Models, parameter, and basic fitting
------------------------------------

* :doc:`ex_parameter_access_1` :download:`(Download notebook) <ex_parameter_access_1.ipynb>`
* :doc:`ex_model_evaluation` :download:`(Download notebook) <ex_model_evaluation.ipynb>`
* :doc:`ex_gaussfit1` :download:`(Download notebook) <ex_gaussfit1.ipynb>`
* :doc:`ex_convenience_fit_1` :download:`(Download notebook) <ex_convenience_fit_1.ipynb>`

On the objective function
-------------------------

The objective function as the function, the value of which is to be minimized.
Generally, the objective function needs to be specified when a minimization is intended.
Models can but do not have to define objective functions. For example, the
GaussFit model defines a chi square objective functions, which is often useful, but
does not have to be used. 


Restrictions of parameter ranges
--------------------------------

A *restriction* in funcFit2 limits the valid range of values of a parameter. Restrictions are common as boundary
conditions in modeling and problems of optimization. 
For example, the width (standard deviation) of a Gaussian must be positive or certain spectral
lines must only occur in absorption or emission for physical reasons.

Restrictions can be handled in many ways. The possibilities to implement restrictions include:

- **Parameter transformations**
  The restriction can be absorbed in the definition of the model, e.g., by using the absolute value of the
  standard deviation in calculating a Gaussian curve. Also mapping of the real number onto a finite range
  such as the `sigmoid function <https://en.wikipedia.org/wiki/Sigmoid_function>`_ can be used.
- **Optimization algorithms** Some optimization algorithm allow to specify boundaries (or more general constraints)
  for the parameter values. One example of such an algorithm is scipy's
  `fmin_l_bfgs_b <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html>`_.
- **Penalty functions** Restrictions can be implemented by penalizing the objective function when the
  boundaries are violated.
  If combined with an optimization algorithm based on gradient descent, it is often helpful to implement
  "soft edges" for penalty, i.e., a strong but finite gradient in the objective, which allows the algorithm
  to "find its way back".


Use penalties to enforce restrictions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Penalties are funcFit2's default mechanism to account for restrictions. A restriction is specified
by a two-valued list with a lower and an upper limit for the parameter range. None can be used to
indicate that no boundary applies on one (or both) sides. Restrictions can be added
to parameters via the `setRestriction` method.

If a parameter value violates the specified boundary restrictions by some margin x, a value of
abs(x)*penaltyFactor is added to the value of the objective function. The default value of
penaltyFactor is 1e20.

* :doc:`ex_restrictviapenalties` :download:`(Download notebook) <ex_restrictviapenalties.ipynb>`


Algorithm allowing boundary conditions (fmin_l_bfgs_b)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here directly invoke the
`fmin_l_bfgs_b <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html>`_
as implemented in scipy to carry out an optimization with boundary conditions

* :doc:`ex_algorestrict` :download:`(Download notebook) <ex_algorestrict.ipynb>`

Use a convenience function to automatically channel the restrictions from the model to the
algorithm 

* :doc:`ex_algorestrict_con` :download:`(Download notebook) <ex_algorestrict_con.ipynb>`



Use a custom objective function
-------------------------------

Custom objective functions can be specified for any model. 

* :doc:`ex_custom_objfct1` :download:`(Download notebook) <ex_custom_objfct1.ipynb>`

Applying relations
------------------

Relations define functional dependences between different parameter values (e.g.,
it may be desirable sometimes to treat to parameters as equal).

* :doc:`ex_apply_relation` :download:`(Download notebook) <ex_apply_relation.ipynb>`


Custom models
-------------

Using custom models is easy.

* :doc:`ex_linmod1` :download:`(Download notebook) <ex_linmod1.ipynb>`
* :doc:`ex_linmod_jit` :download:`(Download notebook) <ex_linmod_jit.ipynb>`


Arithmetic combination of models
--------------------------------

Models can be combined by adding, subtracting, multiplying, or dividing them using the
conventional arithmetic operators. In funcFit2, the operation is actually applied to
the result of the 'evaluate' method of the models. This can be useful in many cases, but
it may fail if, e.g., the calling sequences of the evaluate methods differ or the model
does not have any such method.  

* :doc:`ex_adding_two_gaussians` :download:`(Download notebook) <ex_adding_two_gaussian.ipynb>`


MCMC sampling with emcee
------------------------

* :doc:`ex_emceesample_1` :download:`(Download notebook) <ex_emceesample_1.ipynb>`
