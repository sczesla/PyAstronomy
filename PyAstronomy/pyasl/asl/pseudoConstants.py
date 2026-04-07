import configparser
import os
import PyAstronomy as _PyA


class PseudoConstants:
    """
    A class that parses a constants configuration file and provides
    access to SI values via their symbols as dynamic properties.
    """
    def __init__(self):
        filepath = os.path.join(os.path.dirname(_PyA.__file__), "constants/cdat.dat")
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Could not find the file: {filepath}")

        self._config = configparser.ConfigParser()
        self._config.read(filepath)
        
        # Internal mapping of symbol -> valueSI
        self._constants = {}
        self._load_constants()

    def setSystem(self, s):
        if s != "SI":
            raise Exception("PseudoConstants can only do SI units!")

    def _load_constants(self):
        """Iterates through sections and maps symbols to float values."""
        for section in self._config.sections():
            symbol = self._config.get(section, 'symbol', fallback=None)
            value_si = self._config.get(section, 'valueSI', fallback=None)
            
            if symbol and value_si:
                try:
                    # Clean the value (some might have spaces or scientific notation)
                    self._constants[symbol] = float(value_si.split()[0])
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
    
    
def getPyAConstantsObject():
    """ Get the true PyAConstants object or the pseudo alternative """
    try:
        from PyAstronomy import constants as PC
        return PC.PyAConstants()
    except:
        return PseudoConstants()

    