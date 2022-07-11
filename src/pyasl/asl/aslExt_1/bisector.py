import numpy as np
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy.pyaC import zerocross1d

def bisector(x, y, x0=None, nx=None, ylevels=None, sampling="redblue"):
    """
    """
    # Red and blue (left/right) part
    
    if x0 is None:
        x0 = x[np.argmin(y)]
    
    def bisy(x, y, y0):
        """ Individual bisector point, returns array (bisector, left, right), NaN if not applicable """
        zcs = zerocross1d(x, y-y0) - x0
        zb = zcs[zcs < 0]
        zr = zcs[zcs > 0]
        r = np.array([np.NaN, np.NaN, np.NaN], dtype=float)
        if len(zb) > 0:
            r[1] = max(zb)+x0
        if len(zr) > 0:
            r[2] = min(zr)+x0
        r[0] = (r[1]+r[2])/2.
        return r
    
    def ys(x, y, rb, nx):
        """
        nx Y levels on red or blue side; returns list with y values
        """
        if rb == "blue":
            yls = list(y[x < x0][-nx:])
        elif rb == "red":
            yls = list(y[x > x0][0:nx+1])
        else:
            raise(PE.PyAValError("Unknown 'rb' value.", \
                                 solution="Use 'red' or 'blue'."))
        return yls
    
    if not ylevels is None:
        yls = ylevels
    elif not nx is None:
        yls = []
        if sampling == "red":
            yls.extend(ys(x,y,"red",nx))
        elif sampling == "blue":
            yls.extend(ys(x,y,"blue",nx))
        elif sampling == "redblue":
            yls.extend(ys(x,y,"blue",nx))
            yls.extend(ys(x,y,"red",nx))
        else:
            raise(PE.PyAValError("Unknown string parameter for 'sampling'.", \
                                 solution="Use 'red', 'blue', or 'redblue'"))
    else:
        raise(PE.PyAValError("Could not construct ylevels for calculating bisector.", \
                             solution="Specify 'ylevels' manually or choose 'nx' and sampling parameter."))
        
    yls = np.sort(np.array(yls))
    
    r = np.zeros( (len(yls), 3), dtype=float )
    for i, yl in enumerate(yls):
        r[i,::] = bisy(x, y, yl)
    
    return x0, yls, r
        
        
        