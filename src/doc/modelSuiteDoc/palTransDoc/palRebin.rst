Transit model with finite integration time
============================================

If the instrumental integration times are significant, the observed transit
light-curve can show an apparent distortion, which is the consequence of
integrating over a fraction of the transit (e.g., Kipping [#fKip]_). To take this effect into account,
the *PalFitRebin* class can be used, which calculates the model on a finer
grid and than averages to obtain the actual light curve:

.. currentmodule:: PyAstronomy.modelSuite.XTran.palTrans
.. autoclass:: PalLC_Rebin
   :members:
   
  
.. [#fKip] Kipping 2010, "Binning is sinning: morphological light-curve distortions
           due to finite integration time", 
           2010MNRAS.408.1758K