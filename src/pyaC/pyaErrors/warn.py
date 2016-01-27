from __future__ import print_function
from __future__ import absolute_import
from .pyaErrTemplate import PyaErrTemplate

def warn(w):
  """
    Parameters:
      - `w` - Something which can be converted into a string.
              This may especially be a PyA exception, which is treated as a warning here.
    
    Warnings are printed to stdout.
  """
  print("--------------------")
  print("| PyA User warning |")
  print("--------------------")
  if isinstance(w, PyaErrTemplate):
    print(w.__str__(head=False))
    return
  print(str(w))