from PyAstronomy.funcFit import OneDFit
import numpy

class GeomPlanetBrightPhase(OneDFit):
  
  def __init__(self):
    """
      @TODO - docu
    """
    OneDFit.__init__(self, ["i", "per", "T0", "brat"])
  
  def evaluate(self, time):
    phase = (time - self["T0"])/self["per"]
    phase -= numpy.floor(phase)
    diskBright = 0.5 * ( 1.0 + numpy.cos( (180.0 - phase*360.0)/180.0 * numpy.pi ) * \
                 numpy.sin(self["i"]/180.0*numpy.pi) ) * self["brat"]
    return diskBright 
    