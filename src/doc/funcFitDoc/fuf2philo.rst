The reasoning behind fuf2
=================================

This gives a brief overview of why a new model class is desirable for funcFit.

Why have a successor to OneDFit?
------------------------------------

The model class of funcFit is OneDFit. The central element of OneDFit
is the `evaluate` method, which relates the model parameters to an actual model, :math:`y = f(x)`, such
as a Gaussian curve. This is very convenient if the goal is a minimization
of the sum of squares or chi square or, more generally, any function which
depends on :math:`f(x_i)` and :math:`y_i` alone. Consequently, OneDFit has been
efficiently and successfully used in a great number of application, which it handles
well.

One particular element of OneDFit, which stood the test of time is its concept
of parameter management. This involves addressing the parameters by 'name',
freeze and thaw them in the fitting (and sampling) process, restrict their ranges,
and so on. The fundamentals of this concept shall be preserved.

Unfortunately,
an inconvenient situation arises when the objective function also depends on the
parameter values themselves. This can for example be the case, when some noise
model, such as additional Gaussian jitter, is part of the optimization process. The inconvenience
ultimately arises as a result of the way the objective function is handled automatically
in OneDFit. In fact, such cases can
be managed by custom objective functions, but tend to require additional effort.
Similar inconveniences arise when beyond point estimates credibility intervals or marginal
distributions under complex likelihood and prior functions are to be derived. Yet other
difficulties emerge when data sets not corresponding to the default setup (x, y, error)
are considered.
Such cases were not the focus of funcFit, when the basic structure of OneDFit was developed. 

Probably many if not all cases mentioned can somehow be handled by OneDFit, the price being,
of course, ever more cumbersome and, thus, error prone coding.
To handle such situations more efficiently and embed the model in a natural Bayesian
environment, a new model class is therefore quite desirable.


How to remedy the shortcomings?
----------------------------------

Many of the above-mentioned inconveniences can be alleviated by a slight shift in paradigm,
in particular, by making the objective function a central element of the model class
rather than `evaluate`. This immediately resolves the complications brought about
by objective functions depending also on the parameter values. The ultimate reason is that the
step represented by the `evaluate` method in OneDFit is logically prior to the evaluation
of the objective function. Also, evaluation does not even always strictly require a complete
evaluation of the model.

This change naturally allows to more clearly separate optimization from the actual model class.
While the incorporation of the fitting (and sampling) into the model class was often convenient,
it also helped to create those structural problems, which ultimately led to the redesign of the
model class. Delegating optimization to convenience functions helps to maintain greater flexibility.





