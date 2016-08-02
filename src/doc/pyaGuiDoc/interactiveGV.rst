Interactive Gauss/Voigt line fit
====================================

.. p23ready

The IAGVFit tool allows to interactively fit a series
of Gaussian or Voigt profiles to a given data set.

Once the data are on screen (see example below), a component
can be added by using the middle mouse button. In particular,
three points have to be specified from "left to right":
(1) the outer left points should approximately be placed where
the profile reaches half height/depth; (2) the middle points
should be placed at the "bottom" of the profile; and (3) the outer
right point should be placed at the right-side equivalent of (1).
Points (1) and (2) are used to estimate the width of the indicated
profile and point (3) is used to estimate the area.

The parameter values are shown in the "Parameters" panel. The "active
component" -- by default indicated by a black line in the figure --
is that whose parameter values are currently shown. The checkboxes preceding
"free" indicate whether the associated parameter is considered free in
the fit process.

Parameters can be fitted manually using the mouse wheel. Depending
on the choice in the "mouse wheel" panel, the parameter under
consideration is modified by multiplying with a given factor or
adding/subtracting the amount specified in the aforementioned panel.
The parameter affected by mouse-wheel action is determined using
the outer left checkbox in the "Parameters" panel.

Using the "Set fit range" button, the range of the data to be fitted
can be restricted. After clicking the button, to middle-mouse-button
clicks into the figure are required to set the "left" and "right"
edges of the git range. 

Example: Using the interactive fitter 
------------------------------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyaGui
    from PyAstronomy import funcFit as fuf
    import numpy as np
    
    # Data for the plot
    x = np.linspace(5000., 5010, 200)
    y = np.ones(len(x))
    yerr = np.ones(len(x)) * 0.01
    y += np.random.normal(0., 0.01, len(x))
    
    gf = fuf.GaussFit1d()
    gf["A"] = -0.3
    gf["mu"] = 5004.
    gf["sig"] = 0.2
    y += gf.evaluate(x)
    
    # Create interactive fitter
    igv = pyaGui.IAGVFit(x, y, yerr=yerr, mode="gauss")
    
    r = igv.interactiveFit()
    
    print("Parameters of the fit: ", r[0])
    print("Parameters of active component: ", r[1])
    print("No. of components: ", r[2])


Class API documentation
--------------------------

.. currentmodule:: PyAstronomy.pyaGui
.. autoclass:: IAGVFit
   :members: