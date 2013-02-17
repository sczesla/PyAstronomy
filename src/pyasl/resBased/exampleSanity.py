import unittest

class SanityOfBaraffe98Tracks(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_example(self):
    from PyAstronomy.pyasl import resBased as rb
    import matplotlib.pylab as plt
    import numpy as np
    
    bt = rb.Baraffe98Tracks()
    
    print "Unique metallicity values: ", bt.getUniqueValues("Met")
    print "Unique Y values: ", bt.getUniqueValues("Y")
    print "Unique Lmix values: ", bt.getUniqueValues("Lmix")
    print "Unique mass values: ", bt.getUniqueValues("Mass")
    
    # Get model data and plot log10(age) versus effective temperature
    m = bt.getModelData((0.0, 0.275, 1.0, 0.3))
    plt.plot(np.log10(m.Age*1e9), m.Teff, 'b.-')
    
    # Find all models with metallicity 0.0, Y 0.275, Lmix 0.0,
    # and any mass.
    models = bt.findModels(Met=0.0, Mass=None, Y=0.275, Lmix=1.0)
    # Print out the model parameters.
    print
    print "Number of models found: ", len(models)
    for i, model in enumerate(models):
      print "Model no. %3d : Met = %3.1f, Y = %5.3f, Lmix = %3.1f, Mass = %4.2f" \
            % ((i+1,) + model)
    
    # Finally, show the plot
#    plt.show()

class SanityOfNasaExoplanetArchive(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_example(self):
    from PyAstronomy import pyasl
    import matplotlib.pylab as plt
    
    nexa = pyasl.NasaExoplanetArchive()
    
    # See what information is available
    cols = nexa.availableColumns()
    print
    
    # Get all information for planet 'wasp-12 b'
    # By default, the search is case-insensitive
    print "Entry of Wasp-12 b"
    print nexa.selectByPlanetName("Wasp-12 b")
    
    print
    # Get all data and plot ra vs. dec
    dat = nexa.getAllData()
    plt.plot(dat.ra, dat.dec, 'b.')
#    plt.show()