import unittest

class SanityOfPyasl(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def tearDown(self):
    pass
  
  def sanity_invertIndexSelectionExample(self):
    """
      Example of invertIndexSelection.
    """
    import numpy as np
    from PyAstronomy import pyaC as pc
    
    # Create "data" values and select some
    x = np.exp(np.arange(20.0)/20.0)
    indi = np.where(np.logical_and(x > 1.4, x < 1.7))
    
    print "Selected indices and values:"
    print "  indices: ", indi
    print "  values : ", x[indi]
    
    indiInv = pc.invertIndexSelection(x, indi)
    
    print
    print "Inverted selection:"
    print "  indices: ", indiInv
    print "  values : ", x[indiInv]
    
    # Check that the same result is obtained by simply
    # passing the length of the array `x`
    
    indiInvLen = pc.invertIndexSelection(len(x), indi)
    
    print
    print "Are indiInv and indiInvLen are the same? "
    print "  ", np.all(indiInvLen == indiInv)
  
  def sanity_invertIndexSelection(self):
    """
      Sanity of invertIndexSelection.
    """
    import numpy as np
    from PyAstronomy import pyaC as pc
    x = np.arange(20.0)
    indi = np.where(np.logical_or(x < 5.0, x > 17.0))[0]
    i1 = pc.invertIndexSelection(x, indi)
    self.assertEqual(len(x), len(indi)+len(i1), "indi and inverse do not make up the original array.")
    for i in indi:
      self.assertFalse(i in i1, "found an index of indi in inverse.")
    i2 = pc.invertIndexSelection(len(x), indi)
    self.assertTrue(np.all(i1 == i2), "i1 and i2 are not the same.")
    i3 = pc.invertIndexSelection(x, i1)
    self.assertTrue(np.all(i3 == indi), "re-inverse of i1 does not match indi.")