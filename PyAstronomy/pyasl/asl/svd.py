from __future__ import print_function, division
import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo

class SVD:
  """
    Use Singular Value Decomposition (SVD) to obtain broadening function.
    
    The technique is, e.g., described in Rucinski 1999, "Determination of
    Broadening Functions Using the Singular-Value Decomposition (SVD)
    Technique" and references therein.
  """
  
  def __init__(self):
    self.w = None
    self.des = None
    self.bn = None
  
  def _createDesignMatrix(self, template, bn):
    """
      Creates the "design matrix".
      
      The design matrix is written to the attribute `des`.
      
      Parameters
      ----------
      template : array
          The template spectrum.
      bn : int
          The length (number of bins) used for the
          broadening function.
    """
    # Number of bins used for the design matrix (edges
    # of the template are cut)
    n = len(template) - bn + 1
    self.des = np.matrix(np.zeros( (bn, n) ))
    for i in smo.range(bn):
      self.des[i,::] = template[i:i+n]
    self.des = self.des.T
  
  def getSingularValues(self):
    """
      Access the singular values.
      
      Returns
      -------
      Singular values : matrix
          The singular values.
    """
    if self.w is None:
      raise(PE.PyAOrderError("No singular values available yet.", \
                             "Use the `decompose` method first."))
    return self.w[:]
  
  def getModel(self, broad, asarray=False, modelIndices=False):
    """
      Calculates a model resulting from template and broadening function.
      
      Parameters
      ----------
      broad : array or matrix
          The broadening function
      asarray : bool, optional
          If True, the broadening function will be returned
          as an array instead of a matrix. The default is
          False.
      modelIndices : bool, optional
          If True, the method returns also an array of indices,
          which refer to the template, for which the model is
          calculated. Note that the model cannot be calculated
          at the edges. The default is False.
    """
    if self.des is None:
      raise(PE.PyAOrderError("You need to create a design matrix first.", \
                             solution="Use the `decompose` method."))
    if not modelIndices:
      if not asarray:
        return self.des * np.matrix(broad)
      else:
        return (self.des * np.matrix(broad)).getA()[::,0]
    else:
      modelInd = np.arange(self.bn//2, len(self.template)-self.bn//2)
      if not asarray:
        return (self.des * np.matrix(broad), modelInd)
      else:
        return ((self.des * np.matrix(broad)).getA()[::,0], modelInd)
  
  def getBroadeningFunction(self, flux, wlimit=0.0, asarray=False):
    """
      Calculate the broadening function.
      
      .. note:: The `decompose` method has to be called
                first. On this call, the template is
                specified.
      
      Parameters
      ----------
      flux : array or matrix
          The observed function (flux).
      wlimit : float, optional
          A limit to be applied to the singular
          values. Values smaller than the specified
          limit will be discarded in calculating
          the broadening function.
      asarray : bool, optional
          If True, the broadening function will be returned
          as an array instead of a matrix. The default is
          False.
      
      Returns
      -------
      Broadening function : matrix, array
          The broadening function. The return type can be set
          by the `asarray` flag.
    """
    if len(flux) == (len(self.template) - self.bn + 1):
      validIndi = np.arange(len(self.template) - self.bn + 1)
    elif len(flux) == len(self.template):
      validIndi = np.arange(self.bn//2, len(flux)-self.bn//2)
    else:
      raise(PE.PyAValError("Inappropriate length of flux array.\n" + \
                           "  Either len(template) or len(template)-bn+1 is accepted.", \
                           solution="Adjust the length of the input spectrum."))
    if self.w is None:
      raise(PE.PyAValError("There is no vector of singular values yet.", \
                           solution="Call `decompose` first."))
    w1 = 1.0/self.w
    indi = np.where(self.w < wlimit)[0]
    w1[indi] = 0.0
    b = np.dot(self.v, np.dot(np.diag(w1), np.dot(self.u.T, np.matrix(flux[validIndi]).T)))
    if not asarray:
      return b
    else:
      return b.getA()[::,0]
  
  def getRVAxis(self, binsize, refWvl):
    """
      Calculate a radial velocity axis for the broadening function.
      
      .. note:: Here, it is assumed that the broadening function
                refers to a spectrum in wavelength space.
      
      Parameters
      ----------
      binsize : float
          Size of spectral bins.
      refWvl : float
          The reference wavelength.
    
      Returns
      -------
      Radial velocity axis : array
          An array containing the radial velocity shift
          pertaining to the individual bins in the
          broadening function. The unit is km/s.
    """
    if self.bn is None:
      raise(PE.PyAOrderError("The length of the broadening function has not yet been specified.", \
                             solution="Call `decompose` first, where the length (bn) is specified."))
    result = -np.arange(self.bn) + float(self.bn//2)
    result *= binsize
    result /= refWvl
    result *= 299792.458
    return result
    
    
  
  def decompose(self, template, bn):
    """
      Carry out the SVD of the "design matrix".
      
      This method creates the "design matrix" by applying
      a bin-wise shift to the template and uses numpy's
      `svd` algorithm to carry out the decomposition.
      
      The design matrix, `des` is written in the form:
      "`des` = u * w * transpose(v)". The matrices des,
      w, and u are stored in homonymous attributes.
      
      Parameters
      ----------
      template : array or matrix
          The template with respect to which the broadening
          function shall be calculated.
      bn : int
          Width (number of elements) of the broadening function.
          Needs to be odd.
    """
    if bn % 2 != 1:
      raise(PE.PyAValError("`bn` needs to be an odd number."))
    self.bn = bn
    self.template = template.copy()
    self._createDesignMatrix(template, bn)
    # Get SVD deconvolution of the design matrix
    # Note that the `svd` method of numpy returns
    # v.H instead of v.
    self.u, self.w, v = np.linalg.svd(self.des, full_matrices=False)
    self.v = v.T
    