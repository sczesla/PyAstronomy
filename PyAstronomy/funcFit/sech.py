# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
from .onedfit import OneDFit
from PyAstronomy.pyaC import pyaErrors as PE


class Sech1d(OneDFit):
    """
    Implements a one dimensional hyperbolic secant

    The functional form is:

    .. math:: \\frac{2 A}{e^{(x-mu)/w} + e^{-(x-mu)/w} } + x \\times lin + off

    Here, `lin` and `off` denote the linear and the offset term.

    .. note:: The area under the curve is given by :math:`\\pi \\times A`

    *Fit parameters*:
     - `A` - Amplitude (maximum/minimum of the curve, not area)
     - `mu` - Center of the hyperbolic secant
     - `w` - Width parameter
     - `off` - Offset
     - `lin` - Linear term
    """

    def __init__(self):
        OneDFit.__init__(self, ["A", "mu", "w", "off", "lin"])
        self.setRootName("Sech")

    def evaluate(self, x):
        """
        Evaluates the model for current parameter values.

        Parameters
        ----------
        x : array
            Specifies the points at which to evaluate the model.
        """
        if self["w"] == 0.0:
            raise(PE.PyAValError("Width of Sech must be larger than zero.",
                                 solution="Change width ('w')."))
        z = (x-self["mu"]) / self["w"]
        y = (2. * self["A"]) / (np.exp(z) + np.exp(-z))  + \
            + self["off"] + (self["lin"] * x)
        return y

