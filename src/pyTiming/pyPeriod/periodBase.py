# -*- coding: utf-8 -*-
from __future__ import print_function, division
_mplImport = True
try:
  import matplotlib.pylab as mpl
except ImportError:
  _mplImport = False
import numpy
from PyAstronomy.pyaC import pyaErrors as PE

class PeriodBase:
  """
    Base class for all periodograms.
  
    This class provides the framework
    for all periodograms within *pyPeriod* package.

    *PeriodBase* has a *plot* method, which can be used to provide
    a quick-look of the result. The significance of a
    feature with a power *Pn* can be assessed using the *prob* and
    *FAP* methods.

    Attributes
    ----------
    power : array
        The periodogram power.
    freq : array
        The frequencies at which the power are evaluated.
  """

  def __calcPeriodogram(self):
    pass

  def prob(self,Pn):
    raise(PE.PyANotImplemented("`prob` is not yet implemented."))

  def plot(self, *args,**kwargs):
    """
      Creates a matplotlib figure and axes class instance to visualize the result.
      
      Parameters:
        - `FAPlevels` - optional, List of false-alarm probability (FAP) levels
        - `*args` - optional, Arguments passed to plot method of axes class.
        - `**kwargs` - optional, Keyword arguments passed plot method to axes class.
      
      This method provides a quick and simple way to visualize the results of the a
      periodogram calculation.
      
      Returns:
        The created *Figure* and *Axes* class instances.
    """
    if not _mplImport:
      raise(PE.PyARequiredImport("To use this function matplotlib must be installed.", \
                                 where="PeriodBase::plot", solution=["Install matplotlib.", "Avoid call and use custom plotting interface."]))
    fig = mpl.figure()
    ax = fig.add_subplot(1,1,1)

    if "FAPlevels" in kwargs:
      fapLvs = numpy.array(kwargs["FAPlevels"])
      del kwargs["FAPlevels"]
      powLvs = self.powerLevel(fapLvs)
      for i,powLv in enumerate(powLvs):
        ax.plot([self.freq.min(),self.freq.max()],[powLv, powLv],'k--')
        ax.annotate("FAP "+("%2.3g" % fapLvs[i]), xy=(self.freq.max(),powLv), xytext=None, xycoords='data',
           textcoords='data', arrowprops=None, horizontalalignment="right")
      
    if self.label is None:
      ax.set_title("Periodogram")
      ax.set_xlabel("Frequency")
      ax.set_ylabel("Power")
    else:
      ax.set_title(self.label["title"])
      ax.set_xlabel(self.label["xlabel"])
      ax.set_ylabel(self.label["ylabel"])
    ax.plot(self.freq, self.power, *args, **kwargs)

    return fig, ax

  def FAP(self, Pn):
    """
      Obtain the false-alarm probability (FAP).

      The FAP denotes the probability that at least one out of M
      independent power values in a prescribed search band of a
      power spectrum computed from a white-noise time series is
      as large as or larger than the threshold, `Pn`.
      It is assessed through
      
      .. math:: FAP(Pn) = 1 - (1-Prob(P>Pn))^M \\; ,
      
      where "Prob(P>Pn)" depends on the type of periodogram
      and normalization and is
      calculated by using the *prob* method;
      *M* is the number of independent power
      values and is computed internally.

      Parameters
      ----------
      Pn : float
          Power threshold.
      
      Returns
      -------
      FAP : float
          False alarm probability.
    """
    prob = self.M * self.prob(Pn)
    if prob > 0.01:  return 1.-(1.-self.prob(Pn))**self.M
    return prob
  
  def probInv(self, Prob):
    raise(PE.PyANotImplemented("`probInv` is not yet implemented."))
  
  def powerLevel(self, FAPlevel):
    """
      Power threshold for FAP level.
    
      Parameters
      ----------
      FAPlevel : float or array
            "False Alarm Probability" threshold
    
      Returns
      -------
      Threshold : float or array
          The power threshold pertaining to a specified
          false-alarm probability (FAP). Powers exceeding this
          threshold have FAPs smaller than FAPlevel.
    """
    Prob = 1.-(1.-FAPlevel)**(1./self.M)
    return self.probInv(Prob)
    
  
  def stats(self, Pn):
    """
      Obtain basic statistics for power threshold.
      
      Parameters
      ----------
      Pn : float
          Power threshold.
    
      Returns
      -------
      Statistics : dictionary
          A dictionary containing {'Pn': *Pn*, 'FAP': *FAP(Pn)* ,
          'Prob': *Prob(Pn)*} for the specified power threshold, *Pn*.
    """
    return {'Pn': Pn, 'Prob': self.prob(Pn), 'FAP': self.FAP(Pn)}




class TimeSeries:
  """
    A container class that holds the observed light curve.
    
    Parameters
    ----------
    time : array
        The time array.
    flux : array
        The observed flux/data.
    error : array, optional
        The error of the data values.
  """

  def __init__(self, time, flux, error=None):
    self.time = numpy.array(time)
    self.flux = numpy.array(flux)
    self.error = error

    if len(self.time) != len(self.flux):
      raise(PE.PyAValError("Incompatible dimensions of input data arrays (time and flux).", \
                           where="TimeSeries::__init__", solution="Check the input arrays!"))

  def returnNyquist(self):
    """
      Calculate the average Nyquist frequency.
      
      Returns
      -------
      Nyquist frequency : float
          Half the sampling frequency of the time series.
    """
    return 0.5*1./numpy.mean(self.time[1::]-self.time[0:-1])