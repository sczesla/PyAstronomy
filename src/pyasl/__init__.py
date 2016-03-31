# -*- coding: utf-8 -*-
from __future__ import print_function, division
from PyAstronomy.pyaC import ImportCheck
from PyAstronomy.pyaC import pyaErrors as PE
import six

_ic = ImportCheck(["numpy", "quantities", "scipy", "matplotlib", "pyfits", "astropy.io.fits", "pandas", "six", "astropy"], required=["numpy", "scipy"])
_moduleImportStatus = {}

from PyAstronomy.pyasl.asl import *
from PyAstronomy.pyasl.resBased import *


def moduleImportStatus(mode="all"):
  """
    Provide information on model import status.
    
    Parameters
    ----------
    mode : string, {"all", "fail", "success"}, optional
        If 'all' is given, all modules will be shown. If
        'fail' is given, only failures will be shown and
        the opposite in case of 'success'.
  """
  mss = ("all", "fail", "success")
  if not mode in mss:
    raise(PE.PyAValError("Unknown mode: " + str(mode), \
                         solution="Use any of " + ', '.join(mss)))
  
  paks = set()
  mods = []
  for k in six.iterkeys(_moduleImportStatus):
    paks.add(k[0])
    mods.append(k[1])
  mlpaks = 0
  mlmods = 0
  for p in paks:
    mlpaks = max(mlpaks, len(p))
  for m in mods:
    mlmods = max(mlmods, len(m))
  
  lines = []
  mlline = 0
  allk = sorted(list(_moduleImportStatus.keys()))
  for k in allk:
    if mode == "fail":
      if _moduleImportStatus[k]:
        continue
    elif mode == "success":
      if not _moduleImportStatus[k]:
        continue
    lines.append(("%" + str(mlpaks) + "s  " + "%" + str(mlmods) + "s  %5s") % (k[0], k[1], str(_moduleImportStatus[k])))
    mlline = max(mlline, len(lines[-1]))
  
  sep = "-"*mlline
  head = ("{0: ^" + str(mlline) + "}").format("Import status of 'pyasl'")
  print(sep)
  print(head)
  print(sep)
  for l in lines:
    print(l)
  print(sep)
