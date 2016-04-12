from __future__ import print_function, division
from PyAstronomy import funcFit as fuf
from PyAstronomy import pyasl
from PyAstronomy.pyaC import pyaErrors as PE
import numpy as np
import scipy.interpolate as sci
from PyAstronomy.modelSuite import ic
import six.moves as smo

class LLGauss(fuf.OneDFit):
  """
    A spectral model based on Gaussian profiles and a line list.
    
    This class provides a simple spectral model based on a number
    of Gaussian lines, whose strength may be fitted individually.
    Note that the EW of the lines is given by:
    `A{n}`*`lineScale`, where `A{n}` is the area of the n-th
    Gaussian. The `scale` parameter does not influence the
    EW of the Gaussians.
    
    .. note:: The unit of the EWs given in the `lineList` needs
              to be the same as the wavelength units.
    
    *Fit parameters*:
      - `lineScale` - A common scaling of the area of *all* lines.
      - `scale` - A scaling of the *entire* spectrum.
      - `eps` - Linear limb-darkening coefficient.
      - `vrad` - Radial velocity [km/s].
      - `vsini` - Projected stellar rotational velocity [km/s].
      - `A{n}` - The amplitudes (area) parameters of the individual Gaussians.
      - `sig{n}` - The standard deviations of the individual Gaussian.
      - `mu{n}` - The position of the individual Gaussians.
    
    Parameters
    ----------
    onlyAbs : boolean, optional
        If True (default), restrictions will be applied, which
        prevent emission lines in the spectrum.
    lineList : array
        An array with either two or three columns. The first column
        given the position of the lines, the second gives the EW
        of the lines, and the third---if present---gives the depth
        of the lines. The depth is the maximal depression of the
        continuum, e.g., a value of 0.96 means that the center
        of the line of 4% below the continuum. If the depth is given,
        the width of the individual Gaussians is determined from it,
        unless the `uniSig` is specified.
    uniSig : float, optional
        Use "unified sigma", i.e., the same width for all lines.
        Note that this flag overrules the "depth" column in the
        `lineList`, if it has been specified.
    modelBinsize : float, optional
        Internally, the model should be calculated on a finer grid
        than the actual spectrum. This parameter specifies the used
        bin size, which is 0.005 by default.
    useFastRB : boolean, optional
        Use the "fast" rotational broadening algorithm. This algorithm
        uses a wavelength-independent broadening kernel, which is
        considerably faster than considering the wavelength dependence.
        Setting this flag to False is necessary if you use very
        long wavelength ranges; by default it is True.
    verbose : boolean, optional
        If True, the class print the current parameters during the
        evaluation.
  """

  def __init__(self, lineList, uniSig=None, modelBinsize=0.005, useFastRB=True, verbose=False, onlyAbs=True):
    
    if not ic.check["scipy"]:
      raise(PE.PyARequiredImport("Could not import scipy.", \
                                 solution="Install scipy."))
    # Check whether depths are given
    self._depthsGiven = (len(lineList[0,::]) == 3)
    if (not self._depthsGiven) and (uniSig is None):
      raise(PE.PyAValError("No width and no depth given.", \
                           solution="Specify line depth in `lineList` or use `uniSig`."))
    # Check whether unified sigma shall be used
    # Note: Line depth will be ignored then
    self._uniSig = uniSig 
    
    # Only absorption lines
    self._onlyAbs = onlyAbs
    # Verbosity
    self._verbose = verbose
    # Use fast rotational broadening?
    self._usefastRB = useFastRB
    # Save binsize for the model
    self._modelBinsize = modelBinsize
    # Copy line list
    self._lineList = lineList.copy()
    # Number of lines
    self._nl = len(lineList[::,0])
    # A MultiGaussFit object
    self._mg = fuf.MultiGauss1d(self._nl)
    # Store all parameter names from GaussFit
    pars = list(self._mg.parameters().keys())
    # Remove lin and off from multiGauss keys
    pars.remove("lin")
    pars.remove("off")
    
    pars.extend(["lineScale", "scale", "eps", "vrad", "vsini"])
    fuf.OneDFit.__init__(self, pars)
    # Name the model
    self.setRootName("LineListGauss")
    
    # Assign start values
    for i in range(self._nl):
      p = str(i+1)
      # Assign position
      self._mg["mu"+p] = lineList[i,0]
      self["mu"+p] = lineList[i,0]
      # Assign amplitude/EW
      # Note: Positive EWs correspond to absorption lines
      self._mg["A"+p] = -lineList[i,1]
      self["A"+p] = lineList[i,1]
      # Assign width (sigma)
      if self._depthsGiven and (self._uniSig is None):
        # Depth IS given and no unified sigma specified
        self._mg["sig"+p] = abs(lineList[i,1])/((1.-lineList[i,2]) * np.sqrt(2.*np.pi))
        self["sig"+p] = self._mg["sig"+p] 
      elif not (self._uniSig is None):
        # uniSig IS given
        self._mg["sig"+p] = self._uniSig
        self["sig"+p] = self._mg["sig"+p] 
    
    self["scale"] = 1.0
    self["lineScale"] = 1.0
    
    if self._onlyAbs:
      # Apply restrictions to prevent emission lines
      for i in range(self._nl):
        p = str(i+1)
        self.setRestriction({"A"+p:[0.0, None]})

  def thawLineStrengths(self, wlmin=None, wlmax=None):
    """
      Thaw line strengths.
      
      Thaws parameters of the from A{n}, where
      n is the number of the Gaussian component.
      By default all such parameters will be
      thawed. The selection may, however, be
      influenced by specifying `wlmin` and `wlmax`.
      
      Parameters
      ----------
      wlmin : float, optional
          If specified, only the strength of lines at
          wavelengths larger than this limits will be thawed.
      wlmax : float, optional
          If specified, only the strength of lines at
          wavelengths below this limit will be thawed.
      
      Returns
      -------
      Thawed parameters : list
          A list of thawed parameter names.
    """
    freeAs = []
    for i in smo.range(self._nl):
      p = str(i+1)
      mu = self["mu"+p]
      if wlmin is not None:
        if mu < wlmin: continue
      if wlmax is not None:
        if mu >= wlmax: continue
      freeAs.append("A"+p)
      
    self.thaw(freeAs)
    return freeAs

  def thawLineWidths(self, wlmin=None, wlmax=None):
    """
      Thaw line widths.
      
      Thaws parameters of the from sig{n}, where
      n is the number of the Gaussian component.
      By default all such parameters will be
      thawed. The selection may, however, be
      influenced by specifying `wlmin` and `wlmax`.
      
      Parameters
      ----------
      wlmin : float, optional
          If specified, only the strength of lines at
          wavelengths larger than this limits will be thawed.
      wlmax : float, optional
          If specified, only the strength of lines at
          wavelengths below this limit will be thawed.
      
      Returns
      -------
      Thawed parameters : list
          A list of thawed parameter names.
    """
    freeSigs = []
    for i in smo.range(self._nl):
      p = str(i+1)
      mu = self["mu"+p]
      if wlmin is not None:
        if mu < wlmin: continue
      if wlmax is not None:
        if mu >= wlmax: continue
      freeSigs.append("sig"+p)
      
    self.thaw(freeSigs)
    return freeSigs
    
  def numberOfLines(self):
    """
      Get the number of lines in the model.
      
      Returns
      -------
      Number of lines : int
          Number of Gaussian in the model.
    """
    return self._nl

  def evaluate(self, x):
    """
      Calculates the model for current parameters.

      The "model" is calculated on a wavelength axis with
      binning specified by the `modelBinsize` parameter.
      
      The line positions are Doppler shifted and the resulting
      model is rotationally broadened. Finally, the entire
      model is multiplied by the `scale` parameter to
      account for a global normalization. 

      Parameters
      ----------
      x : array
        The wavelengths at which to calculate the model.
      
      Returns
      -------
      model : array
          The model evaluated at the specified positions.
    """
    if self._verbose:
      print("Evaluating with parameters: ", self.parameters())
    # Calculate how much longer the model wavelength axis should
    # be to ensure that Doppler shift and rotational broadening
    # to not "violate" its edges. 
    maxV = abs(self["vsini"]) + abs(self["vrad"])
    deltaWvl = max(x) * (1.0 + maxV/299792.458) - max(x)
    # The additional summand `_modelBinsize` is needed to ensure
    # that the wavelength axis is long enough even if the velocities
    # are very small. 
    mwvl = np.arange(min(x)-deltaWvl, max(x)+deltaWvl+self._modelBinsize, self._modelBinsize)
    
    mflux = np.ones(len(mwvl))
    # Set parameters of multiGauss
    for i in smo.range(self._nl):
      p = str(i+1)
      # Apply global scaling of amplitudes
      self._mg["A"+p] = -self["A"+p] * self["lineScale"]
      self._mg["sig"+p] = self["sig"+p]
      # Apply doppler shift of lines
      self._mg["mu"+p] = (1.0 + self["vrad"]/299792.458) * self["mu"+p]
    mflux += self._mg.evaluate(mwvl)
    
    if abs(self["vsini"]) > 0.0:
      if self._usefastRB:
        mflux = pyasl.fastRotBroad(mwvl, mflux, self["eps"], np.abs(self["vsini"]))
      else:
        mflux = pyasl.rotBroad(mwvl, mflux, self["eps"], np.abs(self["vsini"]))
    
    mflux *= self["scale"]
    
    # Return model interpolated at input wavelength
    f = sci.interp1d(mwvl, mflux)
    return f(x)
