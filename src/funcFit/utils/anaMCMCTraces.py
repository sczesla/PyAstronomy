# -*- coding: utf-8 -*-
from __future__ import print_function, division
from numpy import mean, median, std
import numpy as np
import re
import os
from scipy.stats import spearmanr, pearsonr
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.funcFit.utils import ic
import itertools
from PyAstronomy import pyaC as PC
import six
import six.moves as smo

try:
  import pymc
except ImportError:
  pass

try:
  import matplotlib.pylab as plt
except ImportError:
  pass



def hpd(trace, cred):
  """
    Estimate the highest probability density interval.
    
    This function determines the shortest, continous interval
    containing the specified fraction (cred) of steps of
    the Markov chain. Note that multi-modal distribution
    may require further scrutiny.
    
    Parameters
    ----------
    trace : array
        The steps of the Markov chain.
    cred : float
        The probability mass to be included in the
        interval (between 0 and 1).
    
    Returns
    -------
    start, end : float
        The start and end points of the interval.
  """
  # Sort the trace steps in ascending order
  st = np.sort(trace)
  
  # Number of steps in the chain
  n = len(st)
  # Number of steps to be included in the interval
  nin = int(n * cred)
  
  # All potential intervals must be 1) continous and 2) cover
  # the given number of trace steps. Potential start and end
  # points of the HPD are given by
  starts = st[0:-nin]
  ends = st[nin:]
  # All possible widths are
  widths = ends - starts
  # The density is highest in the shortest one
  imin = np.argmin(widths)
  return starts[imin], ends[imin]
  

def quantiles(trace, qs):
  """
    Get quantiles for trace.
    
    Parameters
    ----------
    trace : array
        The steps of the Markov chain.
    qs : list or array
        The quantiles in *percent*.
    
    Returns
    -------
    Quantiles : dictionary
        For each quantile, the corresponding value.
  """
  # Sort the trace steps in ascending order
  st = np.sort(trace)
  n = len(st)
  # Convert percent into fractions 
  qfrac = np.array(qs) / 100.0
  # Evaluate quantiles
  result = {}
  for i in smo.range(len(qfrac)):
    result[qs[i]] = st[int(n*qfrac[i])]
  return result


class TraceAnalysis:
  """
    Support to analyze MCMC chains.

    This class provides a number of plotting methods. Note that
    **you still need to call *show()* from pylab** to see the
    result.

    Parameters
    ----------
    resource : string or pymc database object
        If string, it assumed to be the filename of the
        Markov Chain file. Otherwise, it is supposed
        to be a pymc database object. If the filename is
        of the form "*.emcee", it is assumed to be
        a trace produced by emcee. 
    
    Attributes
    ----------
    burn : int
        Applies a burn-in. All iterations
        earlier than `burn` will be neglected.
    
    thin : int
        Applies thinning to each chain.
        Retains every `k` th sample,
        where `k` is an integer value.  
  """

  def _parmCheck(self, parm):
    """
      Checks whether a trace is available for the given parameter.
      If not, it throws an exception.
      
      Parameters
      ----------
      parm : string
          Variable name.
    """
    if not parm in self.tracesDic:
      raise(PE.PyAValError("No trace available for parameter "+parm+".\n  Available parameters: "+', '.join(list(self.tracesDic.keys())), \
                           where="TraceAnalysis"))

  def __plotsizeHelper(self,size):
    """
      Helps to define the optimum plot size for large big-picture plots.
    """
    c=1; r=1
    while c*r<size:
      c+=1
      if c*r>=size: break
      else: r+=1
    return c, r

  def _loadEMCEEChain(self, fn=None, burn=0):
    """
      Load chains from file or adjust parameters.
      
      Based on the loaded data base, burn-in and walker selection
      is applied.
      
      Parameters
      ----------
      fn : string, optional
          If given, the filename of the emcee chain file.
          If not given, the previously saved data base will
          be used.
      burn : int, optional
          Number of burn-in steps.
    """
    if not fn is None: 
      self._emceedat = np.load(fn)
    self.emceepnames = list(self._emceedat["pnames"]) + ["deviance"]
    # Dummy tracesDic
    self.tracesDic = dict(zip( self.emceepnames, [None]*len(self.emceepnames)))
    # Build stateDic
    self.stateDic = {"stochastics":dict(zip( list(self._emceedat["pnames"]), [None]*len(list(self._emceedat["pnames"])))), "sampler":{}}
    
    # Flatten the chain
    s = self._emceedat["chain"].shape
    # s[0] = No. of chains, s[1] = Length of individual chain, s[2] = No. of parameters
    
    if not hasattr(self, "_selectedWalker"):
      # Select all walker ...
      selectedWalker = np.arange(s[0])
      # ... and save selection
      self._selectedWalker = selectedWalker
    # Use saved selection
    selectedWalker = self._selectedWalker
    
    # No. of walker considered in the analysis
    nchains = len(selectedWalker)
    
    # Check selected chains
    if np.any(selectedWalker > s[0]) or np.any(selectedWalker < 0):
      raise(PE.PyAValError("You selected at least one walker beyond the valid range (0 - " + str(s[0]-1) + ").",
                           solution="Adjust walker selection."))
    
    self.emceechain = np.zeros((nchains * (s[1]-burn), s[2]+1))
    self.emceechain[::,0:-1] = self._emceedat["chain"][selectedWalker,burn:,::].reshape(nchains * (s[1]-burn), s[2])
    self.emceelnp = self._emceedat["lnp"][selectedWalker,burn:]
    self.emceelnp = self.emceelnp.reshape(nchains * (s[1]-burn))
    
    # Incorporate "deviance" into the usual chain
    self.emceechain[::,-1] = -2.0 * self.emceelnp
    
    self.stateDic["sampler"]["_iter"] = self.emceechain.shape[0]
    self.stateDic["sampler"]["_burn"] = None
    self.stateDic["sampler"]["_thin"] = None
    
    
  def _checkPackage(self, p):
    """
      Check whether package is availabel and raise exception otherwise.
      
      Parameters
      ----------
      p : string
          Name of package
    """
    if not ic.check[p]:
      raise(PE.PyARequiredImport("The package '" + str(p) + "' is not currently installed.", \
                                 solution="Please install " + str(p)))
  
  def selectWalkers(self, ws):
    """
      Select walkers for emcee chains.
      
      Parameters
      ----------
      ws : list or array of integers
          The walker to be considered in the analysis. Counting
          starts at zero.
    """
    self._selectedWalker = np.array(ws, dtype=np.int)
    if self.dbtype == "emcee":
      # Apply burn-in to individual walkers
      self._loadEMCEEChain(burn=self.burn)
    else:
      raise(PE.PyAAlgorithmFailure("Walkers can only be selected for emcee data bases. Current data base type: " + str(self.dbtype)))
    
  
  def __init__(self, resource, db="pickle"):
    if isinstance(resource, six.string_types):
      # Resource is a filename
      if not os.path.isfile(resource):
        raise(PE.PyAFileError(resource, "ne"))
      self.file = resource
      
      if not (re.match(".*\.hdf5", resource) is None):
        db = "hdf5"
      if not (re.match(".*\.zlib", resource) is None):
        db = "hdf5"
      if not (re.match(".*\.emcee", resource) is None):
        db = "emcee"
        
      if db == "pickle":
        self.dbtype = "pymc"
        self._checkPackage("pymc")
        self.db = pymc.database.pickle.load(resource)
        self.stateDic = self.db.getstate()
        self.tracesDic = self.db._traces
      elif db == "hdf5":
        self.dbtype = "pymc"
        self._checkPackage("pymc")
        self.db = pymc.database.hdf5.load(resource)
        self.stateDic = self.db.getstate()
        self.tracesDic = self.db._traces
      elif db == "emcee":
        self.dbtype = "emcee"
        self._loadEMCEEChain(resource)
      else:
        raise(PE.PyAValError("Database type '"+db+"' is currently not supported."))
      
    elif isinstance(resource, pymc.database.base.Database):
      self.db = resource
    else:
      raise(PE.PyAValError("'resource' must be a filename or a pymc database object.", where="TraceAnalysis::__init__"))
    
    # Set default burn-in and thinning 
    self.burn = 0                             
    self.thin = 1                         

  def __getitem__(self, parm):
    """
      Returns the trace for parameter `parm`.
      
      Parameters
      ----------
      parm : string
          Name of the parameter.
      
      Returns
      -------
      trace : array
          The trace for the parameter.
    """
    self._parmCheck(parm)
    if self.dbtype == "pymc":
      return self.tracesDic[parm].gettrace()[self.burn::self.thin]
    elif self.dbtype == "emcee":
      index = self.emceepnames.index(parm)
      return self.emceechain[self.burn::self.thin,index]

  def __str__(self):
    """
      Prints basic information on the current MCMC sample file.
    """
    info =  "MCMC database - Basic information:\n"
    info += "----------------------------------\n\n"
    if hasattr(self, "file"):
      info += "MCMC sample file: "+self.file+"\n\n"
    info+="  Stochastics:  "+str(list(self.stateDic["stochastics"].keys()))+"\n\n"
    info+="  Sampler:      "+"Iterations: "+str(self.stateDic["sampler"]["_iter"])+"\n"
    info+="                "+"Burn:       "+str(self.stateDic["sampler"]["_burn"])+"\n"
    info+="                "+"Thin:       "+str(self.stateDic["sampler"]["_thin"])+"\n\n"
    info+="  Note: There is more information available using the state() method.\n"
    return info

  def availableParameters(self):
    """
      Returns list of available parameter names.
    """
    return list(self.stateDic["stochastics"].keys())

  def availableTraces(self):
    """
      Returns a list of available PyMC *Trace* objects
    """
    return list(self.tracesDic.values())

  def state(self):
    """
      Returns dictionary containing basic information
      on the sampling process.
    """
    return self.stateDic

  def plotTrace(self, parm, fmt='b-'):
    """
      Plots the trace.
      
      Parameters
      ----------
      parm : string
          The variable name.
      fmt : string, optional
          The matplotlib format string used to plot the trace.
          Default is 'b-'.
    """
    if not ic.check["matplotlib"]:
      PE.warn(PE.PyARequiredImport("To use 'plotTrace', matplotlib has to be installed.", \
                                   solution="Install matplotlib."))
      return
    self._parmCheck(parm)
    plt.plot(self[parm], fmt, label=parm + " trace")
    plt.legend()

  def plotTraceHist(self, parm):
    """
      Plots trace and histogram (distribution).

      Parameters
      ----------
      parm : string
          The variable name.
    """
    if not ic.check["matplotlib"]:
      PE.warn(PE.PyARequiredImport("To use 'plotTraceHist', matplotlib has to be installed.", \
                                   solution="Install matplotlib."))
      return
    self._parmCheck(parm)
    
    plt.subplot(1,2,1)
    self.plotTrace(parm)
    plt.subplot(1,2,2)
    self.plotHist(parm)

  def plotHist(self, parsList=None):
    """
      Plots distributions for a number of traces.

      Parameters
      ----------
      parsList : string or list of strings, optional,
          Refers to a parameter name or a list of parameter names.
          If None, all available parameters are plotted.
    """
    if not ic.check["matplotlib"]:
      PE.warn(PE.PyARequiredImport("To use 'plotHists', matplotlib has to be installed.", \
                                   solution="Install matplotlib."))
      return
    if isinstance(parsList, six.string_types):
      parsList = [parsList]
    tracesDic = {}
    if parsList is not None:
      for parm in parsList:
        self._parmCheck(parm)
        tracesDic[parm] = self[parm]
    else:
      # Use all available traces
      for parm in self.availableParameters():
        tracesDic[parm] = self[parm]

    cols, rows = self.__plotsizeHelper(len(tracesDic))

    for i,[pars,trace] in enumerate(tracesDic.items()):
      if len(tracesDic) > 1:
        plt.subplot(rows, cols, i+1)
      plt.hist(trace, label=pars + " hist")
      plt.legend()

  def plotDeviance(self, parsList=None):
    """
      Plots value of deviance over parameter values encountered during sampling.

      Parameters
      ----------
      parsList : string or list of strings, optional,
          Refers to a parameter name or a list of parameter names.
          If None, all available parameters are plotted.
    """
    if not ic.check["matplotlib"]:
      PE.warn(PE.PyARequiredImport("To use 'plotHists', matplotlib has to be installed.", \
                                   solution="Install matplotlib."))
      return
    if isinstance(parsList, six.string_types):
      parsList = [parsList]
    tracesDic = {}
    if parsList is not None:
      for parm in parsList:
        self._parmCheck(parm)
        tracesDic[parm] = self[parm]
    else:
      # Use all available traces
      for parm in self.availableParameters():
        tracesDic[parm] = self[parm]

    ps = self.__plotsizeHelper(len(tracesDic))

    for i,[pars,trace] in enumerate(six.iteritems(tracesDic), 1):
      plt.subplot(ps[0],ps[1],i)
      plt.xlabel(pars)
      plt.ylabel("Deviance")
      plt.plot(self[pars],self["deviance"],'.')


  def plotHists(self, parsList=None, **histArgs):
    print("TraceAnalysis::plotHists() - Warning! This function is deprecated. Use plotHist() instead.")
    self.plotHist(parsList, **histArgs)


  def hpd(self, parm, trace=None, cred=0.95):
    """
      Calculates highest probability density interval (HPD, minimum width BCI).
      
      Parameters
      ----------
      parm : string
          Name of parameter
      cred : float, optional
          Credibility level. Defaults to 0.95, i.e., the 95\% HPD will
          be calculated.
      trace : array, optional
          If a trace is given, it will be used in the calculation instead of the
          trace for `parm` stored in the class. Note that the parm will be
          ignored in this case!
      
      Returns
      -------
      HPD : tuple
          The lower and upper bound of the credibility interval.
    """
    if trace is None:
      self._parmCheck(parm)
      return hpd(self[parm], cred)
    else:
      return hpd(trace, cred)

  def quantiles(self, parm, qlist=None):
    """
      Quantiles for given trace.
      
      Parameters
      ----------
      parm : string
          Name of parameter
      qlist : list of floats (0-100), optional
          Specifies which quantiles shall be calculated. The default is
          2.5, 25, 50, 75, and 97.5 percent.
      
      Returns
      -------
      Quantiles : dictionary
          For each quantile (in percent) the corresponding value.
    """
    if qlist is None:
      qlist = [2.5, 25, 50, 75, 97.5]
    return quantiles(self[parm],qlist)

  def plotCorr(self, parsList=None, **plotArgs):
    """
      Produces correlation plots.
      
      Parameters
      ----------
      parsList : list of string,  optional
          If not given, all available traces are used.
          Otherwise a list of at least two parameters
          has to be specified.
      plotArgs : dict, optional
          Keyword arguments handed to plot procedure of
          pylab.
    """
    if not ic.check["matplotlib"]:
      PE.warn(PE.PyARequiredImport("To use 'plotCorr', matplotlib has to be installed.", \
                                   solution="Install matplotlib."))
      return
    tracesDic = {}
    if parsList is not None:
      for parm in parsList:
        self._parmCheck(parm)
        tracesDic[parm] = self[parm]
      if len(tracesDic) < 2:
        raise(PE.PyAValError("For plotting correlations, at least two valid parameters are needed.", \
                             where="TraceAnalysis::plotCorr"))
    else:
      # Use all available traces
      for parm in self.availableParameters():
        tracesDic[parm] = self[parm]

    pars = list(tracesDic.keys())
    traces = list(tracesDic.values())

    fontmap = {1:10, 2:8, 3:6, 4:5, 5:4}
    if not len(tracesDic)-1 in fontmap:
      fontmap[len(tracesDic)-1] = 3

    k = 1
    for j in range(len(tracesDic)):
      for i in range(len(tracesDic)):
        if i>j:
          plt.subplot(len(tracesDic)-1,len(tracesDic)-1,k)
          plt.title("Pearson's R: %1.5f" % self.pearsonr(pars[j],pars[i])[0], fontsize='x-small')
          plt.xlabel(pars[j], fontsize='x-small')
          plt.ylabel(pars[i], fontsize='x-small')
          tlabels = plt.gca().get_xticklabels()
          plt.setp(tlabels, 'fontsize', fontmap[len(tracesDic)-1])
          tlabels = plt.gca().get_yticklabels()
          plt.setp(tlabels, 'fontsize', fontmap[len(tracesDic)-1])
          plt.plot(traces[j],traces[i],'.',**plotArgs)
        if i!=j:
          k+=1

  def __hist2d(self,x,y,contour=False,bins=(200,200),cmap="Purples",interpolation='nearest',origin="lower", colors="k"):
    """
      Parameters
      ----------
      bins : tuple of two ints
    """
    H, xedges, yedges = np.histogram2d(x, y, bins)
    extent = [xedges.min(), xedges.max(), yedges.min(), yedges.max()]
    if not contour:
      plt.imshow(H, extent=extent, interpolation=interpolation,origin=origin, cmap=cmap, aspect="auto")
#      plt.colorbar(fontsize=4)
    else:
      CS = plt.contour(H, extent=extent,origin=origin,colors=colors)
      plt.clabel(CS, inline=1, fontsize=10)

  def correlationTable(self, parsList=None, coeff="pearson", noPrint=False):
    """
      Calculate and show parameter correlations
      
      Parameters
      ----------
      parsList : list of strings, optional
          A list of parameters for which to calculate the
          correlation. If None, all available parameters
          will be used.
      coeff : string, {"pearson", "spearman"}, optional
          The coefficient to be used. By default, Pearson's
          correlation coefficient will be used.
      noPrint : boolean, optional
          If True, the table output will be suppressed.
      
      Returns
      -------
      Coefficients : dictionary
          Maps each tuple of two parameter names to the associated
          correlation coefficient.
    """
    if parsList is None:
      pars = self.availableParameters()
    else:
      pars = parsList[:]
    if len(pars) < 1:
      raise(PE.PyAValError("You need to specify at least one parameter!"))
    
    if coeff == "pearson":
      coeFunc = self.pearsonr
    elif coeff == "spearman":
      coeFunc = self.spearmanr
    else:
      raise(PE.PyAValError("Unknown coefficient '" + str(coeff) + "'.", \
            solution="Use either of 'pearson' or 'spearman'."))
    
    corrs = {}
    cc = itertools.combinations_with_replacement(pars, 2)
    for c in cc:
      corrs[c] = coeFunc(c[0], c[1])[0]
    
    if not noPrint:
      # Find width of columns
      colWidth = {}
      for p in pars:
        colWidth[p] = max(len(p), 9)
      head = (" " * (max(colWidth.values()) + 1)) + "|"
      for p in pars:
        head += (" %" + str(colWidth[p]) + "s |") % p
      print(head)
      print("-" * len(head))
      for p2 in pars:
        line = ("%" + str(max(colWidth.values())) + "s |") % (p2)
        for p1 in pars:
          if (p2, p1) in corrs:
            coe = corrs[(p2, p1)]
          else:
            coe = corrs[(p1, p2)]
          line +=  " % 8.6f |" % coe
        print(line)
      print("-" * len(head))
    
    return corrs

  def plotCorrEnh(self, parsList=None, **plotArgs):
    """
      Produces enhanced correlation plots.
      
      Parameters
      ----------
      parsList : list of string,  optional
          If not given, all available traces are used.
          Otherwise a list of at least two parameters
          has to be specified.
      plotArgs : dict, optional
          Keyword arguments handed to plot procedures of
          pylab. The following keywords are available:
          contour,bins,cmap,origin,interpolation,colors
    """
    if not ic.check["matplotlib"]:
      PE.warn(PE.PyARequiredImport("To use 'plotCorr', matplotlib has to be installed.", \
                                   solution="Install matplotlib."))
      return
    tracesDic = {}
    if parsList is not None:
      for parm in parsList:
        self._parmCheck(parm)
        tracesDic[parm] = self[parm]
      if len(tracesDic) < 2:
        raise(PE.PyAValError("For plotting correlations, at least two valid parameters are needed.", \
                             where="TraceAnalysis::plotCorr"))
    else:
      # Use all available traces
      for parm in self.availableParameters():
        tracesDic[parm] = self[parm]

    pars = list(tracesDic.keys())
    traces = list(tracesDic.values())

    fontmap = {1:10, 2:9, 3:8, 4:8, 5:8}
    if not len(tracesDic)-1 in fontmap:
      fontmap[len(tracesDic)-1] = 8

    k = 1
    for j in range(len(tracesDic)):
      for i in range(len(tracesDic)):
        if i>j:
          plt.subplot(len(tracesDic)-1,len(tracesDic)-1,k)
#          plt.title("Pearson's R: %1.5f" % self.pearsonr(pars[j],pars[i])[0], fontsize='x-small')
          plt.xlabel(pars[j], fontsize=fontmap[len(tracesDic)-1])
          plt.ylabel(pars[i], fontsize=fontmap[len(tracesDic)-1])
          tlabels = plt.gca().get_xticklabels()
          plt.setp(tlabels, 'fontsize', fontmap[len(tracesDic)-1])
          tlabels = plt.gca().get_yticklabels()
          plt.setp(tlabels, 'fontsize', fontmap[len(tracesDic)-1])
#          plt.plot(traces[j],traces[i],'.',**plotArgs)
          self.__hist2d(traces[j],traces[i],**plotArgs)
        if i!=j:
          k+=1


  def correlationMatrix(self, toScreen=True, method="pearson", parList=None, covariance=False):
    """
      Calculates the correlation or covariance matrix.
      
      Parameters
      ----------
      parList : list of strings, optional
          The list of parameters used in the calculation.
          If not given, all available parameters will be
          used.
      toScreen : boolean, optional
          If True, the result will be printed to stdout
      method : string, {"pearson", "spearman"}
          The correlation coefficient to be used.
      covariance : boolean, optional
          If True, the covariance will be returned instead
          of the correlation. The default is False.
      
      Returns
      -------
      Parlist : list
          Parameter names in the order used in the calculation.
      Correlation matrix : 2d array
          The correlation matrix
      lines : list of strings
          Formatted version of the correlation matrix in the form
          of a list of strings.
    """
    if parList is None:
      parList = self.availableParameters()
    corFunc = None
    if method == "pearson": corFunc = self.pearsonr
    if method == "spearman": corFunc = self.spearmanr
    if corFunc is None:
      raise(PE.PyAValError("The method "+str(method)+" is currently not supported.", \
            solution="Change method argument e.g. to 'pearson'."))
    for p in parList:
      self._parmCheck(p)
    if covariance:
      # The covariance is requested. In this case, the
      # correlation coefficient has to be multiplied by the
      # standard deviation(s).
      stds = {}
      for p in parList:
        stds[p] = self.std(p)
    # Calculate the matrix
    n = len(parList)
    matrix = np.zeros( (n, n) )
    for i in smo.range(n):
      for j in smo.range(n):
        matrix[i, j] = corFunc(parList[i], parList[j])[0]
        if covariance:
          matrix[i, j] *= (stds[parList[i]] * stds[parList[j]])

    # Format the output
    lines = PC.matrix2doutput(matrix, colNames=parList, rowNames=parList, toScreen=toScreen)
 
    return parList, matrix, lines
    

  def pearsonr(self, parm1, parm2):
    """
      Calculates a Pearson correlation coefficient and the
      p-value for testing non-correlation.
      
      Parameters
      ----------
      parm1, parm2 : string
          The names of the two parameters used in the evaluation. 
      
      Returns
      -------
      Pearson correlation coefficient : float
      p-value : float  
      
      Notes
      -----
      Uses SciPy's *scipy.stats.pearsonr* to evaluate.
      
      The SciPy documentation of scipy.stats.pearsonr:
      
        The Pearson correlation coefficient measures the linear
        relationship between two data sets. Strictly speaking, Pearson's
        correlation requires that each data set be normally distributed. 
        Like other correlation coefficients, this one varies between
        -1 and +1 with 0 implying no correlation. Correlations of
        -1 or +1 imply an exact linear relationship. Positive
        correlations imply that as x increases, so does y. Negative
        correlations imply that as x increases, y decreases.
        The p-value roughly indicates the probability of an uncorrelated
        system producing data sets that have a Pearson correlation at
        least as extreme as the one computed from these data sets.
        The p-values are not entirely reliable but are probably reasonable
        for data sets larger than 500 or so.
    """
    self._parmCheck(parm1)
    self._parmCheck(parm2)
    return pearsonr(self[parm1], self[parm2])

  def spearmanr(self, parm1, parm2):
    """
      Calculates a Spearman rank-order correlation coefficient
      and the p-value to test for non-correlation.
      
      Parameters
      ----------
      parm1, parm2 : string
          The names of the two parameters used in the evaluation.
      
      Returns
      -------
      Spearman rank-order correlation coefficient : float
      p-value : float
      
      Notes
      -----
      Uses SciPy's *scipy.stats.spearmanr* to evaluate.
      
      The SciPy documentation of scipy.stats.spearmanr:
      
        The Spearman correlation is a nonparametric measure of
        the monotonicity of the relationship
        between two data sets. Unlike the Pearson correlation,
        the Spearman correlation does not assume that both data
        sets are normally distributed. Like other correlation coefficients,
        this one varies between -1 and +1 with 0 implying no correlation.
        Correlations of -1 or +1 imply an exact monotonic relationship.
        Positive correlations imply that as x increases, so
        does y. Negative correlations imply that as x increases,
        y decreases. The p-value roughly indicates the probability of
        an uncorrelated system producing data sets that have a Spearman
        correlation at least as extreme as the one computed from these
        data sets. The p-values are not entirely reliable but are
        probably reasonable for data sets larger than 500 or so.
    """
    self._parmCheck(parm1)
    self._parmCheck(parm2)
    return spearmanr(self[parm1], self[parm2])

  def mean(self, parm):
    """
      Calculate mean.
      
      Parameters
      ----------
      parm : string
          Name of parameter.
      
      Returns
      -------
        The mean : float
    """
    self._parmCheck(parm)
    return mean(self[parm])

  def median(self, parm):
    """
      Calculate median.
      
      Parameters
      ----------
      parm : string
          Name of parameter.
      
      Returns
      -------
        The median : float
    """
    self._parmCheck(parm)
    return median(self[parm])

  def std(self, parm):
    """
      Calculate standard deviation.
      
      Parameters
      ----------
      parm : string
          Name of parameter.
      
      Returns
      -------
        The standard deviation : float
    """
    self._parmCheck(parm)
    return std(self[parm])

  def show(self):
    """
      Call *show()* from matplotlib to bring graphs to screen.
    """
    try:
      plt.show()
    except Exception as e:
      PE.warn(PE.PyAUnclassifiedError("Plot could not be shown. The following exception occurred:\n" \
                                      + str(e)))

  def setBurn(self, burn):
    """
      Change value of "post burn-in".
      
      In the case of an emcee trace, the "post burn-in" is
      applied to the trace of all walkers.
      
      Parameters
      ----------
      burn : int
          The number of samples to be neglected.
    """
    self.burn = burn
    if self.dbtype == "emcee":
      # Apply burn-in to individual walkers
      self._loadEMCEEChain(burn=self.burn)
  
  def setThin(self, thin):
    """
      Change value of "post thinning".
      
      Parameters
      ----------
      thin : int
           Applies thinning to each chain.
           Retains every `k` th sample,
           where `k` is an integer value.
      
      Notes
      -----
      Use the "post thinning" to thin out your chains.
    """
    self.thin = thin

  def parameterSet(self, prescription="lowestDev"):
    """
      Find parameter values for a particular prescription.
      
      Parameters
      ----------
      prescription : string, {"lowestDev", "mean", "median"}
          Which parameter set to find. If 'lowestDev' is used,
          the method will return the parameter set pertaining to
          the lowest deviance solution. If 'mean' or 'median' are
          specified, the mean or median parameter values
          determined from the Markov Chains are returned.
      
      Returns
      -------
      Parameter set : dictionary
          A dictionary mapping all parameter names to the value
          derived using the specified prescription.
      Lowest deviance index : int
          The index of the lowest deviance solution. Only
          returned if the prescription is 'lowestDev'
    """
    if prescription == "lowestDev":
      result = {}
      indi = np.argmin(self["deviance"])
      for par in self.availableParameters():
        result[par] = self[par][indi]
      return result, indi
    elif prescription == "mean":
      result = {}
      for par in self.availableParameters():
        result[par] = self.mean(par)
      return result
    elif prescription == "median":
      result = {}
      for par in self.availableParameters():
        result[par] = self.median(par)
      return result
    else:
      raise(PE.PyAValError("Unknown state '" + str(prescription) + "'.", \
            solution="Use either of 'lowestDev', 'mean', or 'median'.", \
            where="parameterSet"))

  def setToState(self, model, state="best", verbose=True):
    """
      Set the parameter values to a certain state.
      
      Parameters
      ----------
      model - fitting object 
          The fitting model object whose parameters will be updated.
      state : {"best", "mean"}, optional
          "best" : Set parameters to the "best fit" state as measured
                   by deviance. This is the default.
          "mean" : Set parameters to mean value of trace.
      verbose : bool, optional
          If False, no output about what is done will be generated
          (default is True).
    """
    if verbose:
      print("Setting model to state: ", state)
    if state == "best":
      # Setting to best state as measured by deviance
      indi = np.argmin(self["deviance"])
      if verbose:
        print("Lowest deviance of ", self["deviance"][indi], " at index ", indi)
      for par in self.availableParameters():
        if not par in six.iterkeys(model.parameters()):
          continue
        model[par] = self[par][indi]
        if verbose:
          print("Setting parameter: ", par, " to value: ", model[par])
    if state == "mean":
      # Setting to "mean" state
      for par in self.availableParameters():
        if not par in six.iterkeys(model.parameters()):
          continue
        model[par] = np.mean(self[par])
        if verbose:
          print("Setting parameter: ", par, " to mean value: ", model[par])
        