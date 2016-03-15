from PyAstronomy import pyaC

_modules = ["baraffe98tracks", "nasaExoplanetArchive", "exoplanetEU", \
           "exoplanetsOrg", "kuruczModels", "fip", "sweet_cat"]

for m in _modules:
  pyaC.pyaimport(m, "PyAstronomy.pyasl.resBased", globals())

# from baraffe98tracks import *
# from nasaExoplanetArchive import *
# from exoplanetEU import *
# from exoplanetsOrg import *
# from kuruczModels import *
# from fip import *
# from sweet_cat import *