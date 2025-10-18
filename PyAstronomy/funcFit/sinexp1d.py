# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
from .onedfit import OneDFit


class SinusFit1d(OneDFit):
    """
    Implements a sinusoidal function of the form:
       A * sin[ 2*pi * (nu * t + phi) ] + off.

    *Fit parameters*:
     - `A` - Amplitude
     - `nu` - Frequency (= 1/Period)
     - `phi` - Phase
     - `off` - Offset
    """

    def __init__(self):
        OneDFit.__init__(self, ["A", "nu", "phi", "off"])
        self.setRootName("Sine")

    def evaluate(self, x):
        """
        Evaluates the model for current parameter values.

        Parameters
        ----------
        x : array
            Specifies the points at which to evaluate the model.
        """
        y = self["A"]*np.sin(2.0*np.pi*(self["nu"]*x+self["phi"])) \
            + self["off"]
        return y


class ExpDecayFit1d(OneDFit):
    """
    Implements an exponential decay function of the form
       A * exp[ - (t-t0)/tau] + off.

    *Fit parameters*:
     - `A` - Amplitude
     - `tau` - Mean lifetime (=1/decay rate)
     - `t0` - Onset time
     - `off` - Continuum offset
    """

    def __init__(self):
        OneDFit.__init__(self, ["A", "tau", "t0", "off"])
        self.setRootName("ExpDecay")

    def evaluate(self, x):
        """
        Evaluates the model for current parameter values.

        Parameters
        ----------
        x : array
            Specifies the points at which to evaluate the model.
        """
        y = np.zeros(len(x)) + self["off"]
        indi = x > self["t0"]
        y[indi] += self["A"]*np.exp(-(x[indi]-self["t0"])/self["tau"])
        return y
