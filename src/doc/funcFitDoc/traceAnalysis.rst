.. _traceAnalysisClass:

Analysis of PyMC Markov-Chains
================================

The class *TraceAnalysis* provides some convenient methods
to carry out an analysis of the Markov-Chains resulting from
an MCMC "fit".

It is mainly a wrapper around functionality already present
more or less deep within PyMC. 

.. currentmodule:: PyAstronomy.funcFit
.. autoclass:: TraceAnalysis
   :members:

.. _traceAnalysisExample:

Example: Usage of the TraceAnalysis class
---------------------------------------------

This example demonstrates the usage of the class.

.. note:: It relies on the Markov-Chain produced by the tutorial example: :ref:`tutMCMCSampler`.
   
::

  from PyAstronomy import funcFit as fuf
  import matplotlib.pylab as mpl
  
  fn = "mcmcExample.tmp"
  
  ta = fuf.TraceAnalysis(fn)
  
  # Print out some information.
  print ta
  
  print "Mean, median, and standard deviation for 'A':"
  print ', '.join([str(s) for s in (ta.mean("A"), ta.median("A"), ta.std("A"))])
  
  print "Check whether amplitude and center are correlated: "
  print "Pearson rank correlation (val, p-value) : %12g, %12g" % ta.pearsonr("A", "mu")
  print "Spearman rank correlation (val, p-value): %12g, %12g" % ta.spearmanr("A", "mu")
  print
  
  # Get the trace of amplitude (as numpy array)
  print "Trace of A:"
  traceA = ta["A"]
  print traceA
  print
  
  print "Credibility intervals for amplitude:"
  for level in [0.68, 0.80, 0.95]:
    hpd = ta.hpd("A", level)
    print "HPD = %12.5f - %12.5f for level %4.2f" % (hpd[0], hpd[1], level)
  
  # Finally, plot trace and distribution of A.
  ta.plotTraceHist("A")
  mpl.show()