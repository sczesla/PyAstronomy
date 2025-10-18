# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
from .onedfit import OneDFit
from PyAstronomy.pyaC import pyaErrors as PE


class GaussFit1d(OneDFit):
    """
    A one-dimensional Gaussian

    The functional form is:

    .. math:: \\frac{A}{\\sqrt{2\\pi\\sigma^2}}e^{-(x-\\mu)^2/(2\\sigma^2)} + x \\times lin + off

    Here, `lin` and `off` denote the linear and the offset term.

    *Fit parameters*:
     - `A` - Area of the Gaussian (formerly called Amplitude)
     - `mu` - Center of the Gaussian
     - `sig` - Standard deviation
     - `off` - Offset
     - `lin` - Linear term
    
    .. note::
        Other parameterizations using FWHM and 'height' of the curve can be used.
    
    Parameters
    ----------
    prm : optional tuple of two strings, {("A", "sig"), ("h", "sig"), ("A", "FWHM"), ("h", "FWHM")}
        Can be used to adapt the parameterization of the Gaussian. By default, the curve is
        parameterized by Area (A) and standard deviation (sig). Alternatively, also the height (h)
        of the curve (h = A/sqrt(2*pi*sig**2)) and the Full Width at Half Maximum (FWHM) can be
        used. The naming of the fit parameters is updated accordingly.
    """

    # Valid parameterizations
    _valid_prms = (("A", "sig"), ("h", "sig"), ("A", "FWHM"), ("h", "FWHM"))
    # FWHM -> STD
    _fwhm_to_sig = 2 * np.sqrt(2 * np.log(2))

    def __init__(self, prm=("A", "sig")):
        if not prm in GaussFit1d._valid_prms:
            raise (
                PE.PyAValError(
                    f"Invalid parameterization (f{str(prm)}).",
                    where="GaussFit1d",
                    solution="Choose any of: " + ", ".join([str(x) for x in GaussFit1d._valid_prms]),
                )
            )

        self._prm = prm
        OneDFit.__init__(self, ["mu", "off", "lin"] + list(prm))
        self.setRootName("Gaussian")

    def _check_zero_width(self):
        """
        Check if width (sig/FWHM) is zero. Raise exception if so.
        """
        if self[self._prm[1]] == 0.0:
            raise (
                PE.PyAValError(
                    "Width of Gaussian must be larger than zero.",
                    solution=f"Change width ('{self._prm[1]}').",
                )
            )

    def evaluate_h_sig(self, x, h, sig):
        """
        Evaluates the model for given height and standard deviation (sig)-

        Parameters
        ----------
        x : array
            Specifies the points at which to evaluate the model.
        
        Returns
        -------
        model : array
        """

        y = (
            h * np.exp(-((self["mu"] - x) ** 2) / (2.0 * sig**2))
            + self["off"]
            + (self["lin"] * x)
        )
        return y

    def evaluate(self, x):
        self._check_zero_width()
        if self._prm == ("h", "sig"):
                return self.evaluate_h_sig(x, self["h"], self["sig"])
        elif self._prm == ("A", "sig"):
                return self.evaluate_h_sig(
                    x, self["A"] / np.sqrt(2 * np.pi * self["sig"] ** 2), self["sig"]
                )
        elif self._prm == ("h", "FWHM"):
                return self.evaluate_h_sig(
                    x, self["h"], self["FWHM"] / GaussFit1d._fwhm_to_sig
                )
        elif self._prm == ("A", "FWHM"):
                return self.evaluate_h_sig(
                    x,
                    self["A"] / np.sqrt(2 * np.pi * (self["FWHM"]/self._fwhm_to_sig) ** 2),
                    self["FWHM"] / GaussFit1d._fwhm_to_sig,
                )


class MultiGauss1d(OneDFit):
    """
    A multicomponent Gaussian with a single linear continuum component.

    The parameters are the same as for the *GaussFit1d*,
    except that all receive a number specifying the Gaussian
    component to which they belong. Therefore, they are, e.g.,
    named `A1`, `mu2`, and so on, only `off` and `lin`
    remain unnumbered.

    Parameters
    ----------
    n : int
        The number if Gaussian components.
    """

    def __init__(self, n):
        self.n = n
        params = ["off", "lin"]
        for i in range(n):
            p = str(i + 1)
            params.extend(["A" + p, "mu" + p, "sig" + p])
        OneDFit.__init__(self, params)
        self.setRootName("MultiGauss")

    def evaluate(self, x):
        """
        Evaluates the model for current parameter values.

        Parameters
        ----------
        x : array
            Specifies the points at which to evaluate the model.
        """
        y = self["off"] + self["lin"] * x
        for i in range(self.n):
            p = str(i + 1)
            y += (
                self["A" + p]
                / np.sqrt(2.0 * np.pi * self["sig" + p] ** 2)
                * np.exp(-((self["mu" + p] - x) ** 2) / (2.0 * self["sig" + p] ** 2))
            )
        return y

    def evalComponent(self, x, p):
        """
        Evaluate the model considering only a single component.

        Parameters
        ----------
        x : array
            The abscissa.
        p : int
            Component number (starts with one).

        Returns
        -------
        single component model : array
            The model considering a single component.
        """
        if p > 0 and p <= self.n:
            p = str(p)
            y = self["off"] + self["lin"] * x
            y += (
                self["A" + p]
                / np.sqrt(2.0 * np.pi * self["sig" + p] ** 2)
                * np.exp(-((self["mu" + p] - x) ** 2) / (2.0 * self["sig" + p] ** 2))
            )
            return y
        else:
            raise (
                PE.PyAValError(
                    "No such component (no. " + str(p) + ")",
                    where="MultiGauss1d::evalComponent",
                    solution="Use value between 1 and " + str(self.n),
                )
            )
