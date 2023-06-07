import numpy as np
import matplotlib.pylab as plt
from PyAstronomy.pyaC import pyaErrors as PE

def sysrem_iter_c(rij, sigij, a, sigij2=None):
    """
    Estimate 'c'
    
    Parameters
    ----------
    rij : 2d array
        Residuals (rij[::,i] = ith of nds observation with nobs data points)
    sigij : 2d array
        Uncertainties
    a : 1d array
        Current estimate of 'a' (length is nobs)
    sigij2 : 2d array, optional
        Square of uncertainties. If given, sigij will be ignores.
    
    Returns
    -------
    c : 1d array
        Estimate of 'c'
    """
    nobs, nds = rij.shape
    
    a2 = a**2
    ss2 = sigij**2 if sigij2 is None else sigij2
    
    # Short form below
    #c = np.zeros(nds)
    #for i in range(nds):
    #    ss = ss2[::,i]
    #    c[i] = np.sum( rij[::,i]*a/ss ) / np.sum( a2/ss )
        
    c = np.sum( ((rij/ss2).T * a).T, axis=0) / np.sum( (a2 * (1/ss2).T).T, axis=0 )
        
    return c

def sysrem_iter_a(rij, sigij, c, sigij2=None):
    """
    Estimate 'a'
    
    Parameters
    ----------
    rij : 2d array
        Residuals (rij[::,i] = ith of nds observation with nobs data points)
    sigij : 2d array
        Uncertainties
    c : 1d array
        Current estimate of 'a' (length is nds)
    sigij2 : 2d array, optional
        Square of uncertainties. If given, sigij will be ignores.
    
    Returns
    -------
    a : 1d array
        Estimate of 'a'
    """
    nobs, nds = rij.shape
    
    c2 = c**2
    ss2 = sigij**2 if sigij2 is None else sigij2
    
    # Short form below
    #a = np.zeros(nobs)
    #for j in range(nobs):
    #    ss = ss2[j,::]
    #    a[j] = np.sum( rij[j,::]*c / ss ) / np.sum( c2/ss )
        
    a = np.sum(rij / ss2 * c, axis=1) / np.sum( c2 * (1/ss2), axis=1 )
    
    return a
    
def sysrem_data_prepare(obs, sigs, ms_obs=True, ms_feat=False):
    """
    Construct data sets for SysRem
    
    Parameters
    ----------
    obs : list of arrays
        Observations
    sigs : list of arrays, optional
        Uncertainties of the observations
    
    Returns
    -------
    rij : 2d array
        The residuals (mean-subtracted observations) so that
        rij[::,i] gives the i-th observation
    sm : 2d array
        The uncertainties in the same format as rij
    a0 : 1d array
        A dummy guess for the 'a' vector (linear between 0 and 1)
    ms_obs : boolean, optional
        Subtract mean from observations (columns of data matrix). Default is True.
    ms_feat : boolean, optional
        Subtract mean from features (e.g., spectral bins). Default is False.
    """
    nds = len(obs)
    nobs = len(obs[0])
    # Residuals and STDs
    rij = np.zeros( (nobs, nds) )
    sm = np.zeros( (nobs, nds) )
    
    if len(obs) != len(sigs):
        raise(ValueError("Obs and sigs must have the same length."))
    
    for n in range(nds):
        if np.isnan(np.min(obs[n])):
            raise(PE.PyAValError(f"NaN value found in observations no. {n}. No NaNs admissible."))
        if np.isnan(np.min(sigs[n])):
            raise(PE.PyAValError(f"NaN value found in error array no. {n}. No NaNs admissible."))
    
    for n in range(nds):
        rij[::,n] = obs[n]
        sm[::,n] = sigs[n]
    
    if ms_obs:
        # Subtract mean for observations/experiments
        for n in range(nds):
            rij[::,n] -= np.mean(rij[::,n])
    if ms_feat:
        # Subtract mean from features/bins
        for n in range(rij.shape[0]):
            rij[n,::] -= np.mean(rij[n,::])

    return rij, sm, np.linspace(0,1,nobs)
    
    
def sysrem_update_rij(rij, a, c):
    """
    Subtract current model from residuals
    
    In matrix notation Tamuz et al. write (Eq. 6) R' = R - c a^T. The operation here
    implemented corresponds to (R' = R - c^T)^T -> R'^T = R^T - a c^T (where ^T is the transpose,
    R is a matrix and a and c are vectors).  
    
    Returns
    -------
    New rij : 2d array
        Updated residuals
    """
    return rij - np.outer(a,c)


class SysRem:
    
    def __init__(self, obs, sigs, ms_obs=True, ms_feat=False, a0=None, ms_warn=True):
        """
        Implementation of the SysRem algorithm.
        
        SysRem was described by
        `Tamuz et al. 2005 (MNRAS 356, 1466) <https://ui.adsabs.harvard.edu/abs/2005MNRAS.356.1466T/abstract>`_
        in the context of correcting systematic effects in samples of
        light curves, but has been applied in other areas such as planetary atmospheric
        spectroscopy.
        
        .. note::
        
            The data (or residuals) are stored such that X[::,0] gives the first observation, i.e.,
            the observation is stored in the first column. The data matrix is, therefore, stored here as
            the transpose of that adopted for the PCA etc..
        
        .. note::
            
            Compared to the presentation by Tamuz et al., the roles of `a` and `c` are exchanged.
        
        Parameters
        ----------
        obs : list of 1d arrays or 2d array
            List of observations (e.g., light curves) or 2d array such that
            obs[::,0] is the first observation. 
        sigs : list of 1d arrays or 2d array
            Uncertainties (same format as obs)
        a0 : 1d array, optional
            Starting values for 'a' parameter (same length as
            the observations). If not given, values linearly
            increasing from 0 to 1 are assumed.
        ms_obs : boolean, optional
            Subtract mean from observations in initial data matrix (columns of data matrix). Default is True.
        ms_feat : boolean, optional
            Subtract mean from features in initial data matrix (e.g., spectral bins). Default is False.
        ms_warn : boolean, optional
            If True (default), a warning will be printed if either ms_obs or ms_feat is True so that the
            data will be manipulated by subtracting mean(s). Set False to suppress warning.
        
        Attributes
        ----------
        rijs : list of 2d arrays
            List of residual arrays. Each call of the `iterate`
            method adds an updated residual array to the list.
        ac : list of tuples of arrays
            The `a` and `c` parameters used to obtain the updated
            residual array from the previous one.
            
        """
        # Convert into list of arrays if 2d array
        if isinstance(obs, np.ndarray) and (obs.ndim == 2):
            obs = [obs[::,i] for i in range(obs.shape[1])]
        if isinstance(sigs, np.ndarray) and (sigs.ndim == 2):
            sigs = [sigs[::,i] for i in range(sigs.shape[1])]
        
        if (ms_obs or ms_feat) and ms_warn:
            PE.warn(PE.PyAValError(f"Be aware: Mean will be subtracted from the data (ms_feat={ms_feat}, ms_obs={ms_obs}).", \
                                   where="SYSREM", \
                                   solution=["Use ms_warn=False to suppress this warning or change ms_feat/obs.", \
                                   "Initial data can be viewed using the attribute rijs[0]."]))
        
        self.rij, self.sm, self.a0 = sysrem_data_prepare(obs, sigs, ms_obs=ms_obs, ms_feat=ms_feat)
        self.sm2 = self.sm**2
        if a0 is not None:
            self.a0 = a0
        
        self.nds = len(obs)
        self.nobs = len(obs[0])
    
        self.rijs = [self.rij]
        self.ac = [(None,None)]
    
    def get_model(self, last=True, a=None, c=None):
        """
        Get model corresponding to a,c vector
        
        Parameters
        ----------
        last : boolean, optional
            If True (default), the model for the last iteration will be
            returned (a,c must not be given then).
        a, c : arrays, optional
            If given, these will be used to calculate the model (`last` must be
            set False)
        """
        
        if (a is not None) and (c is not None) and (not last):
            return np.outer(a, c)
        elif last and (a is None) and (c is None):
            if len(self.rijs) == 0:
                raise(PE.PyAValError("Last model requested, but no model exists yet.", \
                                     where="SysRem", \
                                     solution="Call 'iterate' first."))
            a, c = self.ac[-1][0], self.ac[-1][1]
            return np.outer(a, c)
        else:
            raise(PE.PyAValError("Unknown combination of parameters given.", \
                                 where="SysRem.get_model", \
                                 solution="Specify either a and c or use 'last=True'."))
    
    def get_cumulative_model(self):
        """
        Get the cumulative model (sum of all available models)
        """
        cm = None
        for ac in self.ac:
            if ac[0] is None:
                continue
            m = self.get_model(last=False, a=ac[0], c=ac[1])
            cm = m if cm is None else cm+m
        return cm
    
    def iterate(self, atol=1e-3, rtol=0, imax=1001, a0=None):
        """
        A single SysRem iteration: Remove linear systematic effect
        
        Parameters
        ----------
        atol, rtol : float, optional
            Absolute and relative tolerance for a-c iteration
            required to stop. The pertinent quantity tested is the
            difference between consecutive models divided by the
            uncertainties.
            By default, relative tolerance is ignored.
        imax : int
            Maximum of iterations. Throws exception if convergence
            is not reached earlier.
        a0 : array, optional
            Starting value for SYSREM model PCA-like component. If not
            specified, the value defined during instance initialization
            will be used.
        
        Returns
        -------
        rij : 2d array
            Updated residual matrix with best-fit model subtracted.
        a : array
            Best-fit 'a' values
        c : array
            Best-fit 'c' values
        """
        if a0 is None:
            a = self.a0
        else:
            if a0.shape != self.a0.shape:
                raise(PE.PyAValError("a0 has wrong shape", \
                                     where="SYSREM iterate", \
                                     solution="Needs to be a 1d array with same length as 'observation'."))
            a = a0
        # First ac iteration
        c = sysrem_iter_c(self.rijs[-1], None, a, sigij2=self.sm2)
        a = sysrem_iter_a(self.rijs[-1], None, c, sigij2=self.sm2)
        m = np.outer(a,c)

        converge = False
        for i in range(imax):
            m0 = m.copy()
            c = sysrem_iter_c(self.rijs[-1], None, a, sigij2=self.sm2)
            a = sysrem_iter_a(self.rijs[-1], None, c, sigij2=self.sm2)
            m = np.outer(a,c)
            
            dm = (m-m0)/self.sm
            
            if np.allclose(0,dm,atol=atol,rtol=rtol):
                converge = True
                break

        self.last_ac_iterations = i

        if not converge:
            PE.warn(PE.PyAAlgorithmFailure(f"Iteration (a-c) did not converge with absolute tolerance {atol} and relative tolerance {rtol}. Current value is {np.max(np.abs(dm))}.", \
                                           where="SysRem (iterate)", \
                                           solution="Increase `imax` or adapt atol and/or rtol."))
            
        self.ac.append((a,c))
        self.rijs.append(sysrem_update_rij(self.rijs[-1], a, c))
        return self.rijs[-1], a, c
