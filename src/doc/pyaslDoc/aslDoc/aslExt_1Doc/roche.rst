Roche potential
===============

.. p23ready

.. currentmodule:: PyAstronomy.pyasl

Let there be two masses :math:`m_1` and :math:`m_2` with :math:`m_1 \ge m_2` and
:math:`m_1` centered at :math:`(x,y,z) = (0,0,0)` and :math:`(x,y,z) = (a,0,0)`, where
:math:`a` is the separation between the masses.
Both objects are located in the x,y plane with x-axis counting distance along the
connecting line and z measuring distance perpendicular to the orbital plane.

For synchronous rotation and a circular orbit, the Roche potential
reads (e.g., Hilditch 2001 "An introduction to close binary stars")

.. math::

    \Phi &= -\frac{G m_1}{r_1} - \frac{G m_2}{r_2} - \frac{\omega^2}{2}\left[\left(x - a\frac{m_2}{m_1+m_2}\right)^2 + y^2 \right] \\

with

.. math::    

    r_1 &= \sqrt{x^2 + y^2 + z^2} \\
    r_2 &= \sqrt{(x-a)^2 + y^2 + z^2} \\
    q &= \frac{m_2}{m_1} \\
    M &= m_1 + m_2 \\
    m_1 &= \frac{M}{1+q} \\
    m_2 &= \frac{Mq}{1+q} \\
    \omega^2 &= \frac{GM}{a^3} \;\;\; \mbox{(Kepler's third law)} \\
    x', y', z', r_1', r_2' &= \frac{x}{a}, \frac{y}{a}, \frac{z}{a}, \frac{r_1}{a}, \frac{r_2}{a}
    
we write

.. math::

    \Phi &= -\frac{G M}{(1+q) r_1} - \frac{G M q}{(1+q) r_2} - \frac{GM}{2 a^3}\left[\left(x - \frac{a\,q}{1+q}\right)^2 + y^2 \right] \\
    \Phi &= -\frac{G M}{(1+q) r_1'\,a} - \frac{G M q}{(1+q) r_2'\,a} - \frac{GM}{2 a^3}\left[\left(x'\,a - \frac{a\,q}{1+q}\right)^2 + y'^2\,a^2 \right] \\
    \Phi &= -\frac{GM}{2 a} \left( \frac{2}{(1+q) r_1'} + \frac{2 q}{(1+q) r_2'} + \left(x' - \frac{q}{1+q}\right)^2 + y'^2   \right) \\
    \Phi &= -\frac{GM}{2 a} \; \Phi_n(x', y', z', q)
    
where :math:`\Phi_n` is the dimensionless Roche potential.

Note that moving material in the rotating frame is also subject to Coriolis forces, which are non-conservative
and, therefore, cannot be incorporated into the scalar Roche potential. 

Example of usage
----------------

::
    
    from __future__ import print_function, division
    from PyAstronomy import pyasl
    import numpy as np
    import matplotlib.pylab as plt
    
    x, y = np.linspace(-1.5, 2, 300), np.linspace(-1.6, 1.6, 300)
    xx, yy = np.meshgrid(x, y)
    # Coordinates in orbital plain
    z = 0
    
    # Mass ratio
    q = 0.2
    
    # Get dimensional values of Roche potential
    p = pyasl.rochepot_dl(xx, yy, z, q)
    
    # Positions (and potentials) of Lagrange points
    l1, l1pot = pyasl.get_lagrange_1(q)
    l2, l2pot = pyasl.get_lagrange_2(q)
    l3, l3pot = pyasl.get_lagrange_3(q)
    l4, l5 = pyasl.get_lagrange_4(), pyasl.get_lagrange_5()
    l4pot = pyasl.rochepot_dl(l4[0], l4[1], l4[2], q)
    l5pot = pyasl.rochepot_dl(l5[0], l5[1], l5[2], q)
    
    
    print("Effective (dimensionless) radii of first and second mass")
    print("According to the approximation of Eggleton 1983:")
    r1eff = pyasl.roche_lobe_radius_eggleton(q, 1)
    r2eff = pyasl.roche_lobe_radius_eggleton(q, 2)
    print("    Reff1: %5.3f" % r1eff)
    print("    Reff2: %5.3f" % r2eff)
    print()
    print("Roche volume and effective radius from Monte Carlo integration:")
    mcvol1 = pyasl.roche_vol_MC(q, 1)
    mcvol2 = pyasl.roche_vol_MC(q, 2)
    print("    MC Roche lobe volume 1: %6.4f +/- %6.4f" % (mcvol1[0:2]))
    print("    MC Roche lobe volume 2: %6.4f +/- %6.4f" % (mcvol2[0:2]))
    print("    MC effective radius 1: %6.4f +/- %6.4f" % (mcvol1[2:]))
    print("    MC effective radius 2: %6.4f +/- %6.4f" % (mcvol2[2:]))
    
    plt.contour(p, [l5pot*1.02, l3pot, l2pot, l1pot], colors=['g', 'c', 'b', 'r'], extent=[-1.5, 2, -1.6, 1.6])
    plt.text(l1, 0, 'L1', horizontalalignment='center')
    plt.text(l2, 0, 'L2', horizontalalignment='center')
    plt.text(l3, 0, 'L3', horizontalalignment='center')
    plt.text(l4[0], l4[1], 'L4', horizontalalignment='center')
    plt.text(l5[0], l5[1], 'L5', horizontalalignment='center')
    plt.show()

Functionality and API
---------------------

    - :ref:`rlpot`
    - :ref:`rlpade`
    - :ref:`rlvol`
    - :ref:`rllps`
    

.. _rlpot:

Roche lobe potential
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: rochepot_dl
.. autofunction:: rochepot

.. _rlpade:

Partial derivatives
~~~~~~~~~~~~~~~~~~~

.. autofunction:: ddx_rochepot_dl
.. autofunction:: ddy_rochepot_dl
.. autofunction:: ddz_rochepot_dl

.. _rlvol:

Roche lobe volume and radius
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: roche_lobe_radius_eggleton
.. autofunction:: roche_vol_MC
.. autofunction:: roche_yz_extent

.. _rllps:

Lagrange points
~~~~~~~~~~~~~~~

.. autofunction:: get_lagrange_1
.. autofunction:: get_lagrange_2
.. autofunction:: get_lagrange_3
.. autofunction:: get_lagrange_4
.. autofunction:: get_lagrange_5

Coriolis force
~~~~~~~~~~~~~~

.. autofunction:: coriolisAcceleration

