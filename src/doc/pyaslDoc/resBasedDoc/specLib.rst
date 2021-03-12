Library of stellar spectra
==========================

.. currentmodule:: PyAstronomy.pyasl.resBased.spectralLib

Stellar spectra based on `Kurucz atmospheric models <http://kurucz.harvard.edu/grids.html>`_ and spectral
synthesis with the `stellar spectral synthesis program 'spectrum' <https://www.appstate.edu/~grayro/spectrum/spectrum.html>`_
by R.O. Gray. Please acknowledge both when using the spectra, but note that neither is responsible
for this library.  

The files can also be accessed manually `here <ftp://ftp.hs.uni-hamburg.de/pub/outgoing/czesla/spectrallib/>`_.

Available library:

* 'A' : solar metallicity, 3800-11000 A, Teff = 3500-20000 K, logg = 0-5, sampling 0.01 A


Demonstration of usage
----------------------

:doc:`ex_speclib` :download:`(Download notebook) <ex_speclib.ipynb>`

API
---

.. autoclass:: SpectralLib
   :members:
   
