# -*- coding: utf-8 -*-
from PyAstronomy.funcFit import ic
import os.path
from PyAstronomy.pyaC import pyaErrors as PE

def updateFitsHeader(model, hdr, clobber=False, conf=0.9):
    """
      Update the delivered fits-header with the parameters of the model.

      Parameters:
        - `hdr` - Pyfits header which should be updated with the model parameters.
        - `conf' - Confidence level for MCMC errors.
        - `clobber` - Allows to overwrite parameters which might be already present in the header.
    """
    if not ic.check["pyfits"]:
      raise(PE.PyARequiredImport("pyfits required to use fits file.", where="Params::updateFitsHeader"))
      return
    try:
       # @FIXME The next line can NEVER work
#       x=hdr[p]
       raise(PE.pyaOtherErrors(" Keyword PA_model already present in fits-header, aborting...", where="Params::updateFitsHeader"))
       return
    except:
        pass
    hdr.update('PA_model',model.naming.getRoot(),'PyAstronomy model type')
    for p in model.parameters():
        if clobber==False:
            try:
                x=hdr[p]
                raise(PE.pyaOtherErrors(" Parameter "+str(p)+" already present in fits-header", where="Params::updateFitsHeader"))
                return
            except:
                pass
        hdr.update(p, model[p])
    if ic.check["pymc"]:
        from pymc.utils import hpd
        hdr.update('Conf',conf,"Error Confidence Level")
        for p in model.parameters():
            try:
                v_err=hpd(model.MCMC.trace(p)[:], 1.0-conf)
                p0,p1=str(p+'_e0'),str(p+'_e1')
                if len(p0) > 8:
                    raise(PE.pyaOtherErrors(" Cannot save Error for parameter "+str(p)+" because len(" + p0 + ")>8", where="Params::updateFitsHeader"))
                    return
                hdr.update(p0,v_err[0],"Lower confidence boundary")
                hdr.update(p1,v_err[1],"Upper confidence boundary")
            except:
                pass

    return hdr 

def restoreFromFits(filename,ext=1, verbose=1):
    """
    Restores a fits file.

    Args:
        filename: (str): write your description
        ext: (str): write your description
        verbose: (bool): write your description
    """
    if not ic.check["pyfits"]:
      raise(PE.PyARequiredImport("pyfits required to use fits file.", where="Params::updateFitsHeader"))
      return
    import pyfits
    if not os.path.isfile(filename):
        raise(PE.PyAValError("No such file: "+str(filename), where="restoreFromFits", solution="Give correct filename."))
        return
    ff=pyfits.open(filename)[ext]
    try:
        modelname=ff.header['PA_model']
    except:
        raise(PE.pyaOtherErrors("File ("+str(filename)+") contains no model.", where="restoreFromFits", solution="Use another file."))
        return
    if modelname=='convLine':
        import convLine as cl
        nLines=ff.header['model_p0']
        model=cl.ConvLine(nLines)
    elif modelname=='multiGauss1d':
        from PyAstronomy.funcFit import MultiGauss1d
        nLines=ff.header['model_p0']
        model=MultiGauss1d(nLines)
    for p in model.parameters():
        model[p]=ff.header[p]
    return model
    
        
     