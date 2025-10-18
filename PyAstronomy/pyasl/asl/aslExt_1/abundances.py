import numpy as np
import os
from PyAstronomy.pyaC import pyaErrors as PE
import PyAstronomy.pyasl


class AbundancePatterns:

    def _dataFN(self):
        """
        Returns the name of the data file.
        """
        return os.path.join(os.path.dirname(__file__), 'abundances.dat')

    def __init__(self):
        """
        This class provides widely adopted patterns of (number) abundances.

        Data are from https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/XSabund.html

        The following abundance patterns are available:

          - angr, from Anders E. & Grevesse N. (1989, Geochimica et Cosmochimica Acta 53, 197)
          - aspl, from Asplund M., Grevesse N., Sauval A.J. & Scott P. (2009, ARAA, 47, 481)
          - feld, from Feldman U.(1992, Physica Scripta 46, 202 except for elements not listed which are given grsa abundances),
          - aneb, from Anders E. & Ebihara (1982, Geochimica et Cosmochimica Acta 46, 2363),
          - grsa from Grevesse, N. & Sauval, A.J. (1998, Space Science Reviews 85, 161)
          - wilm from Wilms, Allen & McCray (2000, ApJ 542, 914 except for elements not listed which are given zero abundance)
          - lodd from Lodders, K (2003, ApJ 591, 1220)
          - lgpp from Lodders K., Palme H., Gail H.P. (2009, Landolt-Börnstein, New Series, vol VI/4B) (Photospheric, using Table 4)
          - lpgs from same as lgpp (Proto-solar, using Table 10)

        """
        # elements in the xspec vphabs models
        self.xspecVphabsAbundances = ["He", "C", "N", "O", "Ne", "Na", "Mg", "Al", "Si",
                                      "S", "Cl", "Ar", "Ca", "Cr", "Fe", "Co", "Ni"]
        self._an = PyAstronomy.pyasl.AtomicNo()
        
        # Atomic weight patterns {Z:atomic weight}
        self._awps = {}
        
        try:
            from mendeleev import element
        except ImportError as e:
            pass
        else:
            for z in range(1, 31):
                self._awps[z] = element(z).atomic_weight
    
    def getAtomicWeightPattern(self, awp):
        """
        Return atomic weight pattern
        
        Returns
        -------
        Pattern : dictionary
            Dictionary mapping atomic number to standard atomic weight.
        """
        return self._awps
            
    def availablePatterns(self):
        """
        Returns the names (abbreviations) of the available abundance patterns.

        Returns
        -------
        names : array of strings
            The names of the available abundance patterns.
        """
        names = np.genfromtxt(self._dataFN(), skip_footer=30, dtype=str)[1::]
        return names

    def pattern(self, name, form='array', key="symbol"):
        """
        Get the elemental abundance pattern

        Parameters
        ----------
        name : string
            The abbreviation of the abundance pattern.
        form : string, {array, dict}, optional
            Return the abundances as a plain array or as a dictionary.
            The default is 'array'.
        key : string, {symbol, number}, optional
            If return type is a dictionary, this parameter determined
            whether the elemental symbol or the atomic number is used
            as the key. The default is "symbol".

        Returns
        -------
        abundances : numpy array or dict (depending on ``form'')
            Number abundances (not mass) relative to H.
        """
        names = self.availablePatterns()
        if not name in names:
            raise(PE.PyAValError("No such pattern: '" + str(name) + "'",
                                 solution="Choose existing pattern: " +
                                 ', '.join(names),
                                 where="pyasl.AbundancePatterns"))
        if not form in ["dict", "array"]:
            raise(PE.PyAValError("No such format: '" + str(form) + "'",
                                 solution="Choose either 'dict' or 'array'",
                                 where="pyasl.AbundancePatterns"))
        if not key in ["symbol", "number"]:
            raise(PE.PyAValError("No such key: '" + str(key) + "'",
                                 solution="Choose either 'symbol' or 'number'",
                                 where="pyasl.AbundancePatterns"))

        gi = np.where(names == name)[0]
        # Read abundances
        dd = np.genfromtxt(self._dataFN(), skip_header=1, usecols=gi+1)
        # Read element symbols
        el = np.genfromtxt(self._dataFN(), skip_header=1, usecols=0, dtype=str)
        if form == 'array':
            return dd
        if form == 'dict':
            if key == "symbol":
                r = {k: v for (k, v) in zip(el, dd)}
            elif key == "number":
                r = {self._an.getAtomicNo(k): v for (k, v) in zip(el, dd)}
        return r

    def abundance(self, element, pat="angr"):
        """
        Get the abundance of a specific element.

        Parameters
        ----------
        element : int or string
            If an integer is given, it is interpreted as
            the atomic number. If a string is given, it
            is the elemental symbol (e.g., H, He).
        pat : string, optional
            The name of the abundance table.

        Returns
        -------
        abundance : float
            Number abundance of specified element relative to hydrogen.
        """
        if isinstance(element, int):
            p = self.pattern(pat, form='dict', key="number")
        else:
            p = self.pattern(pat, form='dict', key="symbol")
        return p[element]
        
    def patternByMass(self, name, awp=None, key="symbol"):
        """
        Get abundance of element in units of mass.
        
        Relative number abundances are converted into mass fractions
        based on standard atomic weights.
        
        Parameters
        ----------
        name : string
            The abbreviation of the abundance pattern.
        awp : dictionary, optional
            Atomic weight pattern. Dictionary mapping atomic number to
            standard atomic weight. If None, data for from the mendeleev
            package will be used.
        key : string, {symbol, number}, optional
            This parameter determines
            whether the elemental symbol or the atomic number is used
            as the key in the return dictionary. The default is "symbol".

        Returns
        -------
        Pattern : dictionary
            Dictionary mapping atomic number to the respective fraction of mass
        """
        if (len(self._awps) == 0) and (awp is None):
            raise(PE.PyAValError("No atomic weight pattern available.", \
                                 where="abundanceByMass", \
                                 solution=["Provide pattern via awp", "Install mendeleev package"]))
        
        p = self.pattern(name, form="dict", key="number")
        pms = {}
        for k, v in p.items():
            pms[k] = p[k]*self._awps[k]
        tm = sum(pms.values())
        for k in pms.keys():
            pms[k] = pms[k]/tm
        
        if key == "symbol":
            r = {self._an.getElSymbol(k):v for k, v in pms.items()}
            return r
            
        return pms 
        

        
