Examples
===============

.. p23ready

The following examples demonstrate the functionality of the PyPDM module.

Using the Scanner
--------------------

This example demonstrates the use of the *Scanner* class.

.. IsPyAExample
::

    from __future__ import print_function, division
    # Import PDM module
    from PyAstronomy.pyTiming import pyPDM
    
    # Get Scanner instance
    scanner = pyPDM.Scanner(minVal=0.5, maxVal=1.0, dVal=0.05, mode="period")
    # Print the periods covered by the scanner
    print("Periods: ", end=' ')
    for period in scanner:
      print(period, end=' ')


Carry out a PDM analysis
---------------------------

Here we demonstrate how to use the pyPDM class to carry out a
PDM analysis.

::

    import numpy
    import matplotlib.pylab as plt
    from PyAstronomy.pyTiming import pyPDM
    
    # Create artificial data with frequency = 3,
    # period = 1/3
    x = numpy.arange(100) / 100.0
    y = numpy.sin(x*2.0*numpy.pi*3.0 + 1.7)
    
    # Get a ``scanner'', which defines the frequency interval to be checked.
    # Alternatively, also periods could be used instead of frequency.
    S = pyPDM.Scanner(minVal=0.5, maxVal=5.0, dVal=0.01, mode="frequency")
    
    # Carry out PDM analysis. Get frequency array
    # (f, note that it is frequency, because the scanner's
    # mode is ``frequency'') and associated Theta statistic (t).
    # Use 10 phase bins and 3 covers (= phase-shifted set of bins).
    P = pyPDM.PyPDM(x, y)
    f1, t1 = P.pdmEquiBinCover(10, 3, S)
    # For comparison, carry out PDM analysis using 10 bins equidistant
    # bins (no covers).
    f2, t2 = P.pdmEquiBin(10, S)
    
    
    # Show the result
    plt.figure(facecolor='white')
    plt.title("Result of PDM analysis")
    plt.xlabel("Frequency")
    plt.ylabel("Theta")
    plt.plot(f1, t1, 'bp-')
    plt.plot(f2, t2, 'gp-')
    plt.legend(["pdmEquiBinCover", "pdmEquiBin"])
    plt.show()
