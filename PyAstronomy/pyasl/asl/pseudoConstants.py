import configparser
import os
import PyAstronomy as _PyA


class PseudoConstants:
    """
    A class that parses a constants configuration file and provides
    access to SI values via their symbols as dynamic properties.
    """
    
    def _checkSystem(self, s):
        if not s in ["SI", "cgs"]:
            raise Exception("PseudoConstants can only supports SI und cgs units!")
    
    def __init__(self, unitSystem="cgs"):
        self._checkSystem(unitSystem)
        self._unitSystem = unitSystem
        
        # Internal mapping of symbol -> valueSI
        self._load_constants()

    def setSystem(self, s):
        self._checkSystem(s)
        if self._unitSystem == s:
            return
        self._unitSystem = s
        self._load_constants()
            
    def _load_constants(self):
        """Iterates through sections and maps symbols to float values."""
        filepath = os.path.join(os.path.dirname(_PyA.__file__), "constants/cdat.dat")
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Could not find the file: {filepath}")

        self._config = configparser.ConfigParser()
        self._config.read(filepath)
        
        vkey = {"SI": "valueSI", "cgs":"valueCGS"}[self._unitSystem]
        
        self._constants = {}
        
        for section in self._config.sections():
            symbol = self._config.get(section, 'symbol', fallback=None)
            value = self._config.get(section, vkey, fallback=None)
            
            if symbol and value:
                try:
                    # Clean the value (some might have spaces or scientific notation)
                    self._constants[symbol] = float(value.split()[0])
                except ValueError:
                    # Skip values that cannot be converted to float
                    continue

    def __getattr__(self, name):
        """Allows access to constants via pc.symbol."""
        if name in self._constants:
            return self._constants[name]
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")

    def __repr__(self):
        return f"<PseudoConstants: {len(self._constants)} constants loaded>"
    
    
def getPyAConstantsObject(unitSystem="cgs"):
    """ Get the true PyAConstants object or the pseudo alternative.
        unitSystem defaults to cgs for consistency with true PyAConstants
    """
    try:
        from PyAstronomy import constants as PC
        return PC.PyAConstants(unitSystem=unitSystem)
    except:
        return PseudoConstants(unitSystem=unitSystem)

    