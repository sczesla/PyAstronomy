from __future__ import print_function, division
from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy.pyaC import pyaErrors as PE
import os
from .. import readFitsSpec
import numpy as np
import astropy.io.fits as pyfits


class SpectralLib:
    
    def refreshDescription(self):
        """
        (Re-)Download file description.txt
        
        This file contains an identifier and a brief description of the available spectral libraries.
        """
        self.pfs.downloadToFile(self.baseurl+"/description.txt", os.path.join(self.dpa, "description.txt"), True, self.dlv)
    
    def getREADME(self):
        """
        Get content of README file
        
        Returns
        -------
        README : list of strings
            Content of README file
        """
        rurl = self.baseurl + self.lib + "/README"
        rfn = os.path.join("pyasl", "resBased", "spectrallib", self.lib, "README")
        if not self.pfs.fileExists(rfn):
            self.pfs.downloadToFile(rurl, rfn, True, self.dlv)
        f = self.pfs.requestFile(rfn, 'rt')
        return f.readlines()
    
    def _readDescription(self):
        """
        Read description
        """
        self.description = {}
        ll = self.pfs.requestFile(os.path.join(self.dpa, "description.txt"), 'rt').readlines()
        for l in ll:
            l = l.rstrip("\n")
            s = l.split(":")
            if len(s) == 2:
                self.description[s[0].strip()] = s[1].strip()
    
    def listDescription(self):
        """
        List descriptions of available libraries
        """
        print("-"*80)
        print("--- Description of libraries")
        print("-"*80)
        for k, v in self.description.items():
            print("Key         : ", k)
            print("Description : ", v)
        print("-"*80)

    def refreshInventory(self):
        """
        Re-download inventory.txt file
        """
        self.pfs.downloadToFile(self.invurl, self.invfn, True, self.dlv)        

    def _initLib(self, lib):
        """
        Initialize spectral library
        """
        # Get inventory
        self.invfn = os.path.join("pyasl", "resBased", "spectrallib", lib, "inventory.txt")
        self.invurl = self.baseurl + lib + "/inventory.txt"
        if not self.pfs.fileExists(self.invfn):
            self.pfs.downloadToFile(self.invurl, self.invfn, True, self.dlv)
        # Read inventory
        f = self.pfs.requestFile(self.invfn, 'rt')
        ll = [l.rstrip("\n") for l in f.readlines()]
        self.inv = {}
        for l in ll:
            if l.startswith("#"):
                continue
            s = l.split()
            code = []
            # First is number, last is filename
            for i in range(1,len(s)-1):
                try:
                    x = float(s[i])
                except ValueError:
                    x = s[i]
                code.append(x)
            self.inv[tuple(code)] = s[-1]        
    
    def listInventory(self):
        """
        Print content of inventory
        
        The inventory is stored in the attribute `inv` (dictionary)
        """
        print("-"*80)
        print("--- Inventory")
        print("-"*80)
        for k,v in self.inv.items():
            print(k, v)
        print("-"*80)
    
    def read1dFitsSpec(self, fn):
        """
        Read 1d spectral model
        
        Returns
        -------
        w, f : arrays
            Wavelength [A] and flux
        """
        return readFitsSpec.read1dFitsSpec(fn)
        
    def readMuFits(self, fn):
        """
        Read specific intensity file
        
        Returns
        -------
        
        """
        hl = pyfits.open(fn)
        hkeys = hl[1].header
        w = ((np.arange(hkeys["NAXIS2"]) + 1.0) - hkeys["CRPIX1"]) * hkeys["CDELT1"] + hkeys["CRVAL1"]
        mus = np.array(hl["MUS"].data)
        specs = np.array(hl["I(mu)"].data)
        return w, mus, specs
        
    def requestModel(self, teff, logg, mh=0.0, nex="diskint"):
        """
        Get filename of model spectrum in the library
        
        Model files are downloaded on demand.
        
        Parameters
        ----------
        teff : float
            Effective temperature [K]
        logg : float
            Surface gravity
        mh : float
            Metallicity
        nex : string
            Name extension (e.g., diskint or muspecs)
        
        Returns
        -------
        Spectrum : numpy npz file
            The spectrum in the form of the context of a numpy npz file
        """
        if self.lib == "A":
            p = (float(teff), float(logg), float(mh), nex)
        if not p in self.inv:
            raise(PE.PyAValError("The key '" + str(p) + "' is not in the inventory.", \
                                 where="requestSpectrum", \
                                 solution="Use 'listInventory' to see available combinations."))
        fn = os.path.join(self.dpa, self.lib, self.inv[p])
        if not self.pfs.fileExists(fn):
            surl = self.baseurl + self.lib + "/" + self.inv[p]
            self.pfs.downloadToFile(surl, fn, True, self.dlv)
        fp = self.pfs.getFullPath(fn)
        return fp
    
    def __init__(self, lib="A", refreshDescr=False, dlv=True):
        """
        Access spectral library/ies
        
        Parameters
        ----------
        lib : string, optional
            Library key
        refreshDescr : boolean, optional
            Re-download library description file (default is False)
        dlv : boolean, optional
            Verbose download, default is True
        """
        
        # Download verbose
        self.dlv = dlv
        self.baseurl = "ftp://ftp.hs.uni-hamburg.de/pub/outgoing/czesla/spectrallib/"
        
        self.pfs = pp.pyaFS.PyAFS()
        # Data path
        self.dpa = os.path.join("pyasl", "resBased", "spectrallib")
        self.lib = lib
        
        self.fndescr = os.path.join(self.dpa, "description.txt")
        if (not self.pfs.fileExists(self.fndescr)) or refreshDescr:
            # Get description
            self.refreshDescription()
        self._readDescription()
        
        if not self.lib in self.description:
            raise(PE.PyAValError("No such library key: " + str(self.lib), \
                                 where="SpectralLib", \
                                 solution="Chose one of: " + ','.join(self.description.keys())))
        
        self._initLib(self.lib)
        
        