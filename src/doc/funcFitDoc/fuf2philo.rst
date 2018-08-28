The reasoning behind fuf2
=================================


Why have a successor to OneDFit?
------------------------------------

The model class of funcFit is OneDFit. The central element of OneDFit
is the `evaluate` method, which relates the model parameters to an actual model, :math:`f(x)`, such
as a Gaussian curve. This is very convenient if the goal is a minimization
of the sum of squares or chi square or, more generally, any function which
depends on :math:`f(x_i)` and :math:`y_i` alone.

A more inconvenient situation arises when the objective function also depends on the
the parameter values themselves. This can for example be the case, when some noise
model, such as additional jitter, is part of the optimization process. Such cases can
be managed by custom objective functions, but tend to require some additional effort.
Similar inconveniences arise when beyond point estimates credibility intervals or marginal
distributions under complex likelihood and prior functions are to be derived.

To handle such situations more efficiently and embed the model in a natural Bayesian
environment, a new model class is required.

