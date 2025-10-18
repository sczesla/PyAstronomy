Signal to Noise Ratio (SNR) estimation (Beta Sigma Procedure) 
==============================================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

Estimating the amplitude of noise can often only be done after the
data have been taken, because not all noise sources may be known with sufficient
accuracy (if at all) prior to the measurement. Similarly, the details of
the signal to be measured are usually not known as well, which is the point of
taking the measurements in the first place.

In the analysis of data, a reasonable estimate of the amplitude of noise is
often required to carry out, e.g., a test of goodness-of-fit. Thus arises
the problem of estimating the amplitude (e.g., the variance or standard deviation)
of the noise contribution lacking knowledge about both the noise *and* the actual signal.

A widely-used technique to estimate the amplitude of noise is fitting, e.g., a polynomial to
either all are a subsample of the measurements. The polynomial (or some other function) supposedly
represents the (unknown) signal and the resulting residuals can be studied to estimate the
amplitude of noise; an implementation of such a technique can be found here: :ref:`estimateSNRFCT`.
The :math:`\beta\sigma` procedure provides an approach to noise estimation without the explicit need
to evaluate an approximative function.


Noise estimation via the :math:`\beta\sigma` procedure
----------------------------------------------------------

The :math:`\beta\sigma` procedure is based on the analysis of the distribution of
numerical derivatives and is described in the **corresponding paper** entitled `"A posteriori
noise estimation in variable data sets" <http://cdsads.u-strasbg.fr/abs/2018A%26A...609A..39C>`_.
The idea is the following: If the signal can appropriately be
approximated by an N-th degree polynomial, its (N+1)-th derivative vanishes, leaving only noise
terms to contribute to the value of the numerical derivative. The same idea underlies the
`DER_SNR algorithm <http://www.stecf.org/software/ASTROsoft/DER_SNR/>`_; both are, in fact, identical
for specific settings of the :math:`\beta\sigma` procedure.

The basic assumptions made in the current implementation of the :math:`\beta\sigma` procedure
are (i) Gaussian noise with standard
deviation :math:`\sigma_0` in all analyzed bins and (ii) independent noise contributions in
individual measurements. The assumption of Gaussian noise can be considered a rather weak one. In practice,
a Gaussian distribution often provides a reasonable approximation of the noise properties, although
the noise itself may almost never be expected to be exactly Gaussian. If noise contributions, e.g.,
in adjacent measurements, are not independent, subsamples containing more distant measurements may be
selected (jump parameter) to estimate the noise properties.

Quick start
~~~~~~~~~~~~~~~

Estimating the noise contribution in a data set x, y:

::
    
    from __future__ import print_function
    from PyAstronomy import pyasl
    import numpy as np
    
    # Set up artificial data set
    x = np.linspace(0., 1., 200)
    y = 10.2 + x + 4 * x**2
    # Add some independent, Gaussian noise
    gstd = 0.1
    y += np.random.normal(0., gstd, len(x))
    
    # Estimate noise using robust estimate
    beq = pyasl.BSEqSamp()
    # Define order of approximation (use larger values such as 2,3, or 4 for
    # faster varying or less well sampled data sets; also 0 is a valid order)
    N = 1
    # Define 'jump parameter' (use larger values such as 2,3, or 4 if correlation
    # between adjacent data point is suspected)
    j = 1
    # Estimate noise assuming equidistant sampling (often a good approximation even
    # if data are not strictly equidistant) and robust estimation (often advantageous
    # in working with real data)
    nstd, nstdstd = beq.betaSigma(y, N, j, returnMAD=True)
    
    print("Input noise std: ", gstd)
    print("Estimated noise std = %5.3f +/- %5.3f" % (nstd, nstdstd))
    
.. note:: As the noise properties are usually not known a priori as in this example, it is
          advisable to compare two estimates, e.g., by increasing N by one to corroborate the
          validity of the result.


More detailed examples
-----------------------------

.. toctree::
   :maxdepth: 2
   
   betaSigDoc/betaSig.rst


API documentation
~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2
   
   betaSigDoc/betaSigAPI.rst



     



