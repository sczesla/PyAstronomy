from __future__ import print_function, division
import numpy as np
import six.moves as smo

def coordinateGrid(*args):
  """
    Construct a coordinate array for n-dimensional fitting.
    
    For n-dimensional fitting, `funcFit` requires a mapping from
    array index (i,j,k,...) to coordinate (x,y,z,...).
    This mapping must be given in the form of an array with
    dimension (n1, n2, ..., nd), where ni is the length of
    the i-th coordinate axis and nd is the number of coordinate
    axes. Then, e.g. in the 2d case, [0,0,0] gives the x-coordinate
    of index [0,0] and [0,0,1] given the associated y-coordinate.
    
    Parameters
    ----------
    args : arrays
        An arbitrary number of coordinate arrays.
    
    Returns
    -------
    Coordinate grid : array
        An array of dimension (n1, n2, ..., nd) where
        n1, n2, etc. is the length of the n-th coordinate
        axis and nd is the total number of dimensions.
  """
  # Number of axis
  na = len(args)
  # Length of individual axes
  la = []
  # Number of points in cubes
  npo = 1
  # Shape of output array
  shape = []
  for a in args:
    shape.append(len(a))
    la.append(len(a))
    npo *= la[-1]
  shape.append(na)
  # The resulting output grid
  g = np.zeros( shape )
  # Contains the number of data points in "sub-cubes"
  # Holds: nx, nx*ny, nx*ny*nz, etc.
  npcube = [1]
  for i in smo.range(len(la)):
    npcube.append(npcube[-1]*la[i])
  npcube.pop(0)
  # Loop over all points
  # ai = Array Index
  ai = [0]*na
  for i in smo.range(npo):
    itemp = i
    for j in smo.range(na-1,0,-1):
      ai[j] = itemp // npcube[j-1]
      itemp -= ai[j]*npcube[j-1]
    ai[0] = itemp
    # Assign the coordinates
    for k in smo.range(na):
      g[tuple(ai)][k] = args[k][ai[k]]
  return g