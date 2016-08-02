Observed and apparent altitude
================================

.. p23ready

.. currentmodule:: PyAstronomy.pyasl
.. autofunction:: co_refract_forward
.. autofunction:: co_refract

Example: `co_refract_forward` and `co_refract`
-------------------------------------------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import datetime
    import numpy as np
    
    # Assume, a star is observed at an altitude of 50 degrees
    alt = 50.
    # Now one wants to know the real altitude of the star, i.e.,
    # the altitude corrected for atmospheric refraction.
    print()
    print("Get apparent (real) altitude of a star with observed altitude of " + \
          str(alt) + " degrees")
    print("  ->  Apparent altitude = ", alt - pyasl.co_refract_forward(alt))
    
    print()
    print("You are not observing from sea level, but from an altitude of 5000 meter.")
    print(("Apparent altitude = %9.5f, estimated pressure [mbar] = %9.5f, " + \
          "estimated temperature [K] = %9.5f") % \
          pyasl.co_refract(alt, observer_alt=5000, convert_to_observed=False))
    
    print()
    print("Convert apparent (real) altitude into observed altitude.")
    print("Apparent altitude = " + str(alt) + " degrees", end=' ')
    print(" -> Observed altitude = " + str(pyasl.co_refract(alt, full_output=False,\
                                            convert_to_observed=True)[0]))
    
    print()
    print("The same object observed from different observer altitudes")
    apparentAltitudes = np.repeat(30.0, 10)
    obsalts = np.linspace(0.,5000.,len(apparentAltitudes))
    r = pyasl.co_refract(apparentAltitudes, observer_alt=obsalts, convert_to_observed=True)
    for i in range(len(r[0])):
      print("Observed altitude [deg] = %g, pressure [mbar] = %g, temperature [K] = %g" \
            % (r[0][i], r[1][i], r[2][i]))
