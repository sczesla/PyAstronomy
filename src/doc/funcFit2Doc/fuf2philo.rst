Some philosophy behind funcFit2
===============================

FuncFit2 is thought to be a mostly independent successor of the original
funcFit implementation. A redesign of funcFit had become necessary because
of structural problems in funcFit, arising from early design decisions and
the scope of target problems. A number of related issued are discussed below.

Why have a successor to OneDFit from the original funcFit?
----------------------------------------------------------

The model class of funcFit is OneDFit. The central element of OneDFit
is the `evaluate` method, which relates the vector of model parameters values, :math:`P`,
to an actual model curve, so that :math:`y_i = f(x_i, P)`, where the :math:`y_i` denote the
values of the model curve at points :math:`x_i`. An example would be a Gaussian profile.
A common problem is to compare the model :math:`f(x_i, P)` with some measurements, :math:`m_i`,
with the goal of finding the set of parameter values, which fits best, i.e., the set, which
minimizes some objective function, :math:`g`.
The current formulation in OneDFit is very convenient if the goal is a minimization
of the sum of squares or chi square or, more generally, any function which
depends on :math:`f(x_i, P)` and :math:`m_i` alone, :math:`g(f(x_i, P), m_i)`.
In such cases, OneDFit can be used very efficiently.

Unfortunately,
an inconvenient situation arises when the objective function also depends directly on the
parameter values and not only through the model so that :math:`g'(f(x_i,P),P,m_i)`.
This can for example be the case, when some noise
model, such as additional Gaussian jitter, is part of the objective function. The inconvenience
ultimately arises as a result of the way the objective function is handled automatically
by OneDFit. Naturally, such cases can
be managed by custom objective functions, but this tends to require additional effort.
Similar inconveniences arise when beyond point estimates credibility intervals or marginal
distributions under complex likelihood and prior functions are to be derived. Yet other
difficulties emerge when data sets not corresponding to the default setup (x, y, error)
are considered.
Such cases were simply not the focus of funcFit, when the basic structure of OneDFit
emerged. 

Probably many if not all cases mentioned can somehow be handled by OneDFit, the price being,
of course, ever more cumbersome and, thus, error prone coding.
To handle such situations more efficiently and embed the model in a natural Bayesian
environment, a new model class is therefore quite desirable.

How to remedy the shortcomings?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many of the above-mentioned inconveniences can be alleviated by a slight shift in paradigm,
in particular, by making the objective function a central element of the model class
rather than the `evaluate` method. This immediately resolves the complications brought about
by objective functions depending also on the parameter values. The ultimate reason is that the
step represented by the `evaluate` method in OneDFit is logically prior to the evaluation
of the objective function. Actually evaluation does not even always strictly require a complete
evaluation of the model. The fundamental shift is therefore a shift of focus on to the objective function
rather than model evaluation.

This change naturally allows to more clearly separate optimization from the actual model class.
While the incorporation of the fitting (and sampling) into the model class was often convenient,
it also helped to create those structural problems, which ultimately led to the redesign of the
model class. Delegating optimization to convenience functions helps to maintain greater flexibility
and eases the incorporation with the infrastructure of libraries such as SciPy. 

One element of OneDFit, which stood the test of time, is its concept
of parameter management. This involves addressing the parameters by 'name',
freeze and thaw them in the fitting (and sampling) process, restrict their ranges,
and so on. One of the goals is, therefore, to preserve this design tenet.

