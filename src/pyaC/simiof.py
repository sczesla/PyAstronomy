import datetime
import re
import pyaErrors as PE

class SimIOF(file):
  """
    Simple Input/Output file.
    
    If a file is opened for writing, the origin and
    date will be added at the top of the file. If it
    is opened for reading, the properties will be read,
    converted into float of possible, and stored in the
    `args` attribute. If the float-conversion fails,
    the value is kept as a string.
    
    Parameters
    ----------
    origin : str
        Identifies the script/program writing the file.
    args : tuple
        Passed to the constructor of a file object.
    kwargs : dictionary
        Passed to the constructor of a file object.
  """ 
  
  def __init__(self, origin, *args, **kwargs):
    self.origin = origin
    file.__init__(self, *args, **kwargs)
    if self.mode == 'w':
      self.write("# File cerated by " + str(origin) + "\n")
      self.write("# Date: " + str(datetime.datetime.now()).split('.')[0] + "\n")
      self.write("#\n")
    elif self.mode == 'r':
      self._readProps()

  def __del__(self):
    self.close()

  def _readProps(self):
    """
      Reads the properties from input file.
    """
    self.args = {}
    for l in self:
      r = re.match("^#\s+([^:]+):\s*(.*)$", l)
      if r is not None:
        p = r.group(2)
        try:
          p = float(p)
        except ValueError:
          pass
        self.args[r.group(1)] = p
    self.seek(0)

  def _addOneProp(self, n, v, fmt=None):
    """
      Adds a property to an output file.
      
      Parameters
      ----------
      n : string
          Name of the property
      v : any
          Value of that property
      fmt : string, optional
          The format-string used to represent the property.
    """
    if fmt is None:
      self.write("# " + str(n) + ": " + str(v) + "\n")
    else:
      if fmt.startswith("%"):
        self.write("# " + str(n) + ((": " + fmt) % v) + "\n")
      else:
        self.write("# " + str(n) + ((": " + "%" + fmt) % v) + "\n")

  def addProp(self, name, value, fmt=None):
    """
      Add a property to the file.
      
      Parameters
      ----------
      name : string or list of strings
          Name of the property
      value : any, or list of any type
          Value of that property
      fmt : string, list, optional
          The format-string used to represent the property.
          If a single string is given, it is applied to all
          properties. Otherwise, one format string for every
          property should be given. 
    """
    if isinstance(name, basestring):
      # There is only one value to be added
      self._addOneProp(name, value, fmt=fmt)
    else:
      # There is more than one
      for i, (n, v) in enumerate(zip(name, value)):
        if (fmt is None) or isinstance(fmt, basestring):
          # There is only no or one format string
          self._addOneProp(n, v, fmt=fmt)
        else:
          # There is a list of formats
          self._addOneProp(n, v, fmt=fmt[i])