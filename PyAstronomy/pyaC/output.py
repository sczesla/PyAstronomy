from __future__ import print_function, division
import six.moves as smo
from . import pyaErrors as PE
import numpy as np

def matrix2doutput(m, oformat="% 6.2e", colsep=" | ", rowNames=None, colNames=None, transpose=False, toScreen=True):
  """
    Format a matrix in readable form and write it to screen.
    
    The column is specified by the second index, e.g., the first
    entry in the second column is given by m[0,1]. The first
    entry in the third row is, consequently, given by m[2,0].
    
    Parameters
    ----------
    m : 2-dimensional array
        The data to be formatted.
    oformat : string or list of strings, optional
        The output format. If string, the same format string will
        be used for all columns. If a list of strings is given,
        the associated specifier will be used for each individual
        column.
    colsep : string, optional
        The separator used between columns.
    rowNames : list of strings, optional
        The names of the rows.
    colNames : list of strings optional
        The names of the columns.
    transpose : boolean, optional
        If True, the input matrix will be transposed. In effect,
        this exchanges the roles of columns and rows. Note, however,
        that the role of `colNames` and `rowNames` remains
        unaltered. The default is False.
    toScreen : boolean, optional
        If True (default), the result will be written to screen. 
    
    Returns
    -------
    Formatted matrix : list of strings
        The formatted output is a list of strings, which
        might be written to screen.
  """
  if transpose:
    m = np.transpose(m)
  
  ncol = len(m[0,::])
  nrow = len(m[::,0])
  if rowNames is not None:
    if nrow != len(rowNames):
      raise(PE.PyAValError("Number of row names does not match the size of the matrix.", \
                           where="matrix2doutput", \
                           solution="You may try to use `transpose=True'."))
  if colNames is not None:
    if ncol != len(colNames):
      raise(PE.PyAValError("Number of column names does not match the size of the matrix.", \
                           where="matrix2doutput",\
                           solution="You may try to use `transpose=True'."))
  
  # Check whether format is a list (one format for each column)
  if isinstance(oformat, list):
    pass
  else:
    oformat = [oformat] * ncol
  
  if len(oformat) != ncol:
    raise(PE.PyAValError("The number of format specifiers does not match the number of columns.", \
                         where="matrix2doutput",\
                         solution=["Adjust the number of specifiers.", \
                                   "You may try to use `transpose=True'."]))                          
  
  # Find width of columns
  if colNames is not None:
    colWidth = [len(x) for x in colNames]
  else:
    colWidth = [0] * ncol
  if rowNames is not None:
    # Width of "zeros" column with row names
    colWidth.insert(0, max([len(x) for x in rowNames]))
  else:
    # Use -1, if there are no row names
    colWidth.insert(0, -1) 
  # Determine width of rest of columns
  for i in smo.range(ncol):
    l = max([len(oformat[i] % x) for x in m[::,i]])
    colWidth[i+1] = max(colWidth[i+1], l)
  
  lines = []
  
  if colNames is not None:
    # Create header
    if colWidth[0] != -1:
      h = " " * colWidth[0] + colsep
    else:
      h = ""
    for i, n in enumerate(colNames):
      h += ("%" + str(colWidth[i+1]) + "s") % n
      if i != (len(colNames) - 1):
        h += colsep
    lines.append(h)
    lines.append("-" * len(h))
  # Add the data
  for i in smo.range(nrow):
    l = ""
    if rowNames is not None:
      # First, add the name of the row
      l += ("%" + str(colWidth[0]) + "s") % rowNames[i]
      l += colsep
    # Loop over the column indices
    for j in smo.range(ncol - 1):
      l += ("%" + str(colWidth[j+1]) + "s") % (oformat[j] % m[i, j]) + colsep
    # Add last one
    l += ("%" + str(colWidth[-1]) + "s") % ((oformat[-1] % m[i, -1]))
    lines.append(l)
  
  if colNames is not None:
    lines.append(lines[1])
  
  if toScreen:
    for l in lines:
      print(l)
  return lines

  