Converting between effective temperature and stellar color
===============================================================

.. p23ready

Conversion between effective temperature and stellar color according
to :ref:`Ramirez and Melendez <Ramirez2005>` (several bands,
metallicity dependent) and :ref:`Ballesteros 2012 <Ballesteros2012>`
(black-body approximation).

.. _Ramirez2005:

Conversion according to Ramirez and Melendez 2005
-----------------------------------------------------

In their 2005 publication, Ramírez and Meléndez (ApJ 626, 465-485) present
metallicity-dependent relations between stellar effective temperature
and color. Based on these relations, the class Ramirez2005 allows
to convert between effective temperature and color. All 17 color
indices given by the authors can be used.

Example: 
~~~~~~~~~~~~~~~~~

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    # Create class instance
    r = pyasl.Ramirez2005()
    
    # Which color bands are available
    print("Available color bands: ", r.availableBands())
    
    # Convert B-V to effective temperature and back
    bv = 0.75
    feh = 0.0
    teff = r.colorToTeff("B-V", bv, feh)
    bv1 = r.teffToColor("B-V", teff, feh)
    # Watch out for differences between input bv and the output bv1
    print("B-V = ", bv, ", Teff = ", teff, ", bv1 = ", bv1, ", bv-bv1 = ", bv-bv1)



.. _Ballesteros2012:

Conversion according to Ballesteros 2012
---------------------------------------------

Ballesteros 2012 (EPL 97, 34008) present a conversion between
effective temperature and B-V color index based on a black body
spectrum and the filter functions.

Comparison to Ramirez and Mendelez 2005
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below, a comparison between the effective temperatures derived using
the Ballesteros 2012 and Ramirez and Mendelez 2005 procedures is given.
Solar metallicity main-sequence stars were assumed in the conversion.
Clearly, the difference reaches about 200 K for hot stars in the 7000 K
range and becomes smaller for cooler stars.

::
  
  from __future__ import print_function, division
  from PyAstronomy import pyasl
  
  b = pyasl.BallesterosBV_T()
  r = pyasl.Ramirez2005()
  
  # Convert B-V to effective temperature and back
  for bv in [0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45]:
    tr = r.colorToTeff("B-V", bv, 0.0)
    tb = b.bv2T(bv)
  
    print(("B-V [mag] = {3:4.2f} : Teff (R05) = {0:4.0f} K, " + \
            "Teff (B12) = {1:4.0f} K, dTeff = {2: 4.0f} K").format(tr, tb, tr - tb, bv))


  Output:
  -------
  
  B-V [mag] = 0.35 : Teff (R05) = 6952 K, Teff (B12) = 7158 K, dTeff = -206 K
  B-V [mag] = 0.45 : Teff (R05) = 6453 K, Teff (B12) = 6625 K, dTeff = -171 K
  B-V [mag] = 0.55 : Teff (R05) = 6033 K, Teff (B12) = 6170 K, dTeff = -138 K
  B-V [mag] = 0.65 : Teff (R05) = 5672 K, Teff (B12) = 5778 K, dTeff = -106 K
  B-V [mag] = 0.75 : Teff (R05) = 5358 K, Teff (B12) = 5436 K, dTeff =  -78 K
  B-V [mag] = 0.85 : Teff (R05) = 5082 K, Teff (B12) = 5134 K, dTeff =  -53 K
  B-V [mag] = 0.95 : Teff (R05) = 4835 K, Teff (B12) = 4866 K, dTeff =  -31 K
  B-V [mag] = 1.05 : Teff (R05) = 4612 K, Teff (B12) = 4626 K, dTeff =  -13 K
  B-V [mag] = 1.15 : Teff (R05) = 4410 K, Teff (B12) = 4409 K, dTeff =    1 K
  B-V [mag] = 1.25 : Teff (R05) = 4225 K, Teff (B12) = 4213 K, dTeff =   13 K
  B-V [mag] = 1.35 : Teff (R05) = 4055 K, Teff (B12) = 4034 K, dTeff =   21 K
  B-V [mag] = 1.45 : Teff (R05) = 3897 K, Teff (B12) = 3870 K, dTeff =   27 K



Example: 
~~~~~~~~~~~~~~~~~

::

    from __future__ import print_function, division
    from PyAstronomy import pyasl
    
    b = pyasl.BallesterosBV_T()
    
    bv = 0.65
    
    # Convert B-V into effective temperature
    teff = b.bv2T(0.65)
    print("B-V = {0:4.2f} mag -> Teff = {1:4.0f} K".format(bv, teff))
    
    # Convert effective temperature into B-V color
    teff = 4568.0
    bv = b.t2bv(teff)
    print("Teff = {0:4.0f} K -> B-V = {1:4.2f} mag".format(teff, bv))

  
API documentation (Ramirez2005)
------------------------------------
  
.. currentmodule:: PyAstronomy.pyasl
.. autoclass:: Ramirez2005
   :members:
   :private-members:
   
API documentation (BallesterosBV_T)
--------------------------------------

.. autoclass:: BallesterosBV_T
   :members:
   :private-members: