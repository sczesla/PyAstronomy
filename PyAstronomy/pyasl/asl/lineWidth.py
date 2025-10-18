import numpy as np

def convertDampingConstant(gamma, wavelength):
    """
      Convert damping constant into line-width in wavelength units.
      
      The inverse lifetime, :math:`1/\\tau`, is related to the damping constant,
      :math:`\\Gamma`, and the Einstein coefficients, A, according to  
      
      .. math::
       
          \\frac{1}{\\tau_{lu}} = \\Gamma_{lu} = \\sum_{j<u} A_{lj} + \sum_{j<l} A_{uj}\\,,
      
      where u and l denote the upper and lower state of the transition (see, e.g.,
      "Der neue Kosmos", Unsoeld und Baschek, 7th edition).
      The damping coefficient, `gamma`, determines the FWHM in units of
      circular frequency.
      The natural line width (FWHM of the Lorentzian profile) in wavelength
      units is then given by
      
      .. math::
       
          \\frac{\\lambda^2}{2\\pi c}\\Gamma\\,.

      Parameters
      ----------
      gamma : float 
          The damping constant [1/s]. The damping constant may be determined
          by summing all relevant Einstein coefficients as shown above.
      wavelength : float
          Wavelength of the transition in A.
      
      Returns
      -------
      Natural line-width : float
          The Full Width at Half Maximum of the Lorentzian [cm]. 
    """
    # The factor of 1e-8 transforms Angstrom into cm
    return gamma * (wavelength*1e-8)**2 / (2.0*np.pi*29979245800.)