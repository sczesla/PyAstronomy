Signal to Noise Ratio (SNR) estimation (Beta Sigma Procedure) 
==============================================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

Estimating the amplitude of noise can often only be done after the
data have been taken, because not all noise sources may be known with sufficient
accuracy (if at all) prior to the measurement. Similarly, the details of
the signal to be measured are usually not known as well; otherwise a
measurement would be dispensable.

In the analysis of data, a reasonable estimate of the amplitude of noise is
usually required to carry out, e.g., a test of goodness-of-fit. Thus arises
the problem of estimating the amplitude (e.g., the variance or standard deviation)
of the noise contribution lacking knowledge about the actual signal.

A widely-used technique to estimate the amplitude of noise is fitting, e.g., a polynomial to
either all are a subsample of the measurements. The polynomial (or some other function) supposedly
represents the (unknown) signal and the resulting residuals can be studied to estimate the
amplitude of noise; an implementation of such a technique can be found here: :ref:`estimateSNRFCT`.

Noise estimation via the :math:`\beta\sigma` procedure
----------------------------------------------------------

The :math:`\beta\sigma` procedure is based on the analysis of the distribution of
numerical derivatives and is described in TBD. The idea is the following: If the signal can be
approximated by an N-th degree polynomial, its (N+1)-th derivative vanishes, leaving only noise
terms to contribute to the value of the numerical derivative. The same idea underlies the
`DER_SNR algorithm <http://www.stecf.org/software/ASTROsoft/DER_SNR/>`_; both are, in fact, identical
for specific settings of the :math:`\beta\sigma` procedure.

API documentation
~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2
   
   betaSigDoc/betaSigAPI.rst



     



