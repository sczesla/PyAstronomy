from PyAstronomy.pyasl import KeplerEllipse
from PyAstronomy import funcFit as fuf
from numpy import pi
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE


class KeplerEllipseModel(fuf.OneDFit):
    """
    A model of a Keplerian orbit.

    This class uses the *KeplerEllipse* from the PyA's pyasl
    to calculate a Keplerian orbit. It may be used to fit
    complete 3d position or velocity information on the orbit;
    any individual axes may also be selected.

    The constructor allows to specify *relevant axes*, which
    are those axes considered in the calculation. The actual
    (technical) model is, however, only one-dimensional. The values
    returned by evaluate have the order
    a1, b1, c1, a2, b2, c3, ... . Where a, b, and c represent
    the first, second, and third axis and the number specifies
    the data point. Note that in this case, the resulting
    model has not the same number of points as the time
    axis.

    *Fit parameters*
      -  `a`    - The semi-major axis (same units as the data)
      - `per`   - The period (same time units as data)
      - `e`     - The eccentricity
      - `tau`   - Time of periapsis passage (same time units as data)
      - `Omega` - Longitude of the ascending node [deg]
      - `w`     - Argument of periapsis [deg]
      - `i`     - Inclination angle [deg]

    Parameters
    ----------
    relevantAxes : string
        A string containing any combination of x, y, and z.
        The string specifies the axes (and their order) to
        be considered in the calculations.
    mode : string, {"pos", "vel"}
        Determines whether the output is positions or
        velocities. In this case, the units are determined
        by the units of major axis and time (e.g., AU per day).
    """

    def __init__(self, relevantAxes="xyz", mode="pos"):
        self.ke = KeplerEllipse(1.0, 1.0)
        fuf.OneDFit.__init__(self, ["a", "per", "e", "tau", "Omega", "w", "i"])
        self["a"] = 1.0
        self["per"] = 1.0
        # Which axes to consider?
        # x=0, y=1, z=2
        self.axes = ()
        for i, axis in enumerate("xyz"):
            if axis in relevantAxes:
                self.axes += (i,)
        # Save the mode
        if not mode in ["pos", "vel"]:
            raise(PE.PyAValError("Unknown mode: " + str(mode), \
                                 where="KeplerEllipseModel", \
                                 solution="Choose either 'pos' or 'vel'."))
        self._mode = mode

    def evaluate(self, t):
        """
        Calculates and returns model according to the
        current parameter values.

        Although more than one axis may be relevant
        the output will be one dimensional. If, e.g.,
        the relevant axes are x and y, the order of
        the output will be x0, y0, x1, y1, ... .

        Parameters
        ----------
        t : array
            Times at which to evaluate the model.
        """
        self.ke.i = self["i"]
        self.ke.w = self["w"]
        self.ke.Omega = self["Omega"]
        self.ke.e = abs(self["e"])
        self.ke.a = self["a"]
        self.ke.per = self["per"]
        self.ke.tau = self["tau"]
        self.ke.n = 2.0 * pi / self["per"]
        # Get the data pertaining to the relevant axes
        if self._mode == "pos":
            result = self.ke.xyzPos(t)[::, self.axes]
        else:
            result = self.ke.xyzVel(t)[::, self.axes]
        # Reshape to 1d form
        # If the relevant axes are x and y, the order
        # of the output will be x0, y0, x1, y1, x2, y2, ...
        result = result.reshape(result.size)
        return result




class KeplerRVModel(fuf.OneDFit):
    """
    A model of a Keplerian orbit.

    This class uses the *KeplerEllipse* from the PyA's pyasl
    to calculate radial velocities for a Keplerian orbit. All calculations are
    based on the negligible companion mass hypothesis (i.e., this model has been
    implemented with planetary systems in mind).
    
    .. note:: Any planet with zero period will be ignored.

    .. note:: ? is a placeholder for an integer larger zero, indicating the number of
              the planet (the total number is controlled by the `mp` keyword).
              
    .. note:: RVs should be given in m/s and time stamps in days

    *Fit parameters*

      - `per?`   - The period in days
      - `e?`     - The eccentricity
      - `tau?`   - Time of periapsis passage 
      - `w?`     - Argument of periapsis [deg]
      - `K?`     - Semi-amplitude of radial velocity [m/s]
      - `mstar` - Stellar mass in solar masses. This parameter is usually not fitted.
                  It may be used to take into account the uncertainty on stellar mass
                  in Bayesian (MCMC) analysis.

    *Derived parameters (not to be fitted)*
      -  `a?`        - The semi-major axis in AU
      -  `MA?`       - Mean anomaly corresponding to time of first data point [deg]
      -  `msini?`    - Minimum mass (msini) in Jupiter masses 

    Parameters
    ----------
    mp : int, optional
        The number of planets considered in the model. Default is one. Note that
        signals are added, i.e., no interaction is taken into account.
    deg : int, optional
        Default is zero (i.e., a constant). The degree of a polynomial used to represent
        a systematic (non-periodic) evolution in the data.
    msun : float, optional
        Solar mass [kg]
    mJ : float, optional
        Jupiter mass [kg]
    au : float, optional
        Astronomical unit [m]
    """

    def __init__(self, mp=1, deg=0, msun=1.988547e30, mJ=1898.6e24, au=1.49597870700e11):
        self._msun = msun
        self._mJ = mJ
        self._au = au
        self._deg = deg
        self._mp = mp
        
        self.poly = fuf.PolyFit1d(deg)
        # Kepler models for individual planets
        self.kems = [KeplerEllipseModel(relevantAxes="z", mode="vel") for _ in range(1,mp+1)]
        
        pars = []
        for m in range(1, mp+1):
            # Independent parameters
            pars.extend([pn+str(m) for pn in ["K", "per", "e", "tau", "w"]])
            # Dependent parameters
            pars.extend([pn+str(m) for pn in ["msini", "a", "MA"]])
        # There is only one stellar mass
        pars.append("mstar")
        # Use parameters from polynomial
        pars.extend(self.poly.availableParameters())
        
        fuf.OneDFit.__init__(self, pars)
        self.setRootName("KeplerRVModel")
        # Dummy SMA and inclination
        for m in range(mp):
            self.kems[m]["a"] = 1.0
            self.kems[m]["i"] = 90.0
        # Default stellar mass = 0 solar mass (force thoughtful input)
        self["mstar"] = 0.0
        self.setRestriction({"mstar":[0,None]})
        
        for m in range(1, 1+mp):
            self.setRestriction({"K"+str(m):[0.,None], "per"+str(m):[0.,None], "e"+str(m):[0,1]})
    
    def _dtr(self, d):
        """ Convert degrees into rad """
        return d/180.*np.pi
    
    def evaluate(self, t):
        if self["mstar"] == 0.0:
            raise(PE.PyAValError("Stellar mass has to be specified before evaluation.", \
                                 where="KeplerRVModel"))
        
        # Take into account polynomial
        for i in range(self._deg+1):
            p = "c" + str(i)
            self.poly[p] = self[p]
        rvnew = self.poly.evaluate(t)
        
        for m in range(1, self._mp+1):
            # Name add
            nadd = str(m)
            if self["per"+nadd] == 0.0:
                continue
            for p in ["per", "e", "tau", "w"]:
                # Collect parameters from current model
                self.kems[m-1][p] = self[p+nadd]
        
                rvmodel = self.kems[m-1].evaluate(t)
            # Prevent dividing by zero
            imax = np.argmax(np.abs(rvmodel))
            cee = np.cos(self.kems[m-1].ke.trueAnomaly(t[imax]) + self._dtr(self.kems[m-1]["w"])) \
                  + self["e"+nadd]*np.cos(self._dtr(self["w"+nadd]))
            rvnew += rvmodel / rvmodel[imax] * cee * self["K"+nadd] 

            self["msini"+nadd] = self._getmsini(nadd)
            self["MA"+nadd] = self._MA(m)
            self["a"+nadd] = self._geta(nadd)
        
        return rvnew
                                                                                        
    def _getmsini(self, nadd):
        """
        Get msini
        """
        msini = self["K"+nadd] * ((self["per"+nadd]*86400.) * (self["mstar"]*self._msun)**2 / (2.*np.pi*6.67408e-11))**(1./3.) \
                * np.sqrt(1. - self["e"+nadd]**2)
        return msini/self._mJ
    
    def _MA(self, m):
        """
        Get mean anomaly corresponding to first data point [deg]
        """
        return self.kems[m-1].ke.meanAnomaly(0.0)/np.pi*180. % 360.0
    
    def _geta(self, nadd):
        """
        Get SMA in AU
        """
        return ((self["per"+nadd]*86400.)**2 * 6.67408e-11 * (self["mstar"] * self._msun) / (4.*np.pi**2) )**(1./3.) / self._au