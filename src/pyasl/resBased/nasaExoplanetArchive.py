from __future__ import print_function, division
from PyAstronomy.pyaC import pyaPermanent as pp
from PyAstronomy.pyaC import pyaErrors as PE
import PyAstronomy.pyaC as pyaC
import os
import gzip
import csv
import numpy as np
import six


class NasaExoplanetArchive(pp.PyAUpdateCycle):
    """
    Easy access to NASA's exoplanet archive.

    This class downloads a table of exoplanet data from
    the "NASA Exoplanet Archive"
    (http://exoplanetarchive.ipac.caltech.edu/index.html)
    and provides access to these data. By default, the
    data will be re-downloaded every seven days.

    The following data are provided:

    ===========   ===================================  ======
    Column        Description                          Unit 
    -----------   -----------------------------------  ------
    pl_hostname   Name of host star                         
    pl_name       Name of the planet                        
    pl_letter     Planet letter (e.g., b, c, d, etc.)       
    ra            Right ascension                      deg  
    dec           Declination                          deg  
    pl_orbper     Planetary orbital period             d    
    pl_massj      Planetary mass                       MJ   
    pl_radj       Planetary radius                     RJ   
    pl_trandep    Central depth of transit             %    
    pl_trandur    Transit duration                     d    
    pl_tranmid    Transit midpoint                     BJD  
    pl_orbsmax    Semi-major-axis                      AU   
    pl_orbincl    Orbital inclination of planet        deg  
    st_rad        Stellar radii                        Solar
    st_dist       Distance to star                     pc
    st_mass       Stellar mass                         Solar
    st_teff       Effective temperature of star        K    
    st_vsini      Stellar vsin(i)                      km/s 
    st_logg       Stellar surface gravity              cm/s**2
    st_acts       Stellar S-Index                           
    st_vj         Stellar V-band brightness            mag
    ===========   ===================================  ======
    """

    def _downloadData(self):
        """
        Download data and store it to file in PyA's data directory.
        """        
        urlRoot = "http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?"
        table = "&table=PS"
        select = "&select="
        for v in six.itervalues(self._columns):
            select = ''.join([select, ',', v[0]])
        outformat = "&format=csv"

        url = ''.join((urlRoot, table + select, outformat))
        self._fs.downloadToFile(url, self.dataFileName, clobber=True, verbose=False,
                                openMethod=gzip.open)

    def __init__(self):
        self.data = None
        self.dataFileName = os.path.join("pyasl", "resBased", "NEXA.csv.gz")
        configFileName = os.path.join("pyasl", "resBased", "NEXA.cfg")
        pp.PyAUpdateCycle.__init__(self, configFileName, "NEXA")
        # Define columns to select
        # Column name, Description, Unit
        self._columns = {}
        self._columns[1] = ["pl_name", "Planet Name", "", "U15"]
        self._columns[2] = ["hostname", "Host Name", "", "U15"]
        self._columns[3] = ["default_flag", "Default Parameter Set", "", np.int]
        self._columns[4] = ["sy_snum", "Number of Stars", "", np.int]
        self._columns[5] = ["sy_pnum", "Number of Planets", "", np.int]
        self._columns[6] = ["discoverymethod", "Discovery Method", "", "U50"]
        self._columns[7] = ["disc_year", "Discovery Year", "", np.float]
        self._columns[8] = ["disc_facility", "Discovery Facility", "", "U100"]
        self._columns[9] = ["soltype", "Solution Type", "", "U50"]
        self._columns[10] = ["pl_controv_flag", "Controversial Flag", "", np.int]
        self._columns[11] = ["pl_refname", "Planetary Parameter Reference", "", "U200"]
        self._columns[12] = ["pl_orbper", "Orbital Period", "days", np.float]
        self._columns[13] = ["pl_orbpererr1", "Orbital Period Upper Unc.", "days", np.float]
        self._columns[14] = ["pl_orbpererr2", "Orbital Period Lower Unc.", "days", np.float]
        self._columns[15] = ["pl_orbperlim", "Orbital Period Limit Flag", "", np.float]
        self._columns[16] = ["pl_orbsmax", "Orbit Semi-Major Axis", "au])", np.float]
        self._columns[17] = ["pl_orbsmaxerr1", "Orbit Semi-Major Axis Upper Unc.", "au", np.float]
        self._columns[18] = ["pl_orbsmaxerr2", "Orbit Semi-Major Axis Lower Unc.", "au", np.float]
        self._columns[19] = ["pl_orbsmaxlim", "Orbit Semi-Major Axis Limit Flag", "", np.float]
        self._columns[20] = ["pl_rade", "Planet Radius", "Earth Radius", np.float]
        self._columns[21] = ["pl_radeerr1", "Planet Radius Upper Unc.", "Earth Radius", np.float]
        self._columns[22] = ["pl_radeerr2", "Planet Radius Lower Unc.", "Earth Radius", np.float]
        self._columns[23] = ["pl_radelim", "Planet Radius Limit Flag", "", np.float]
        self._columns[24] = ["pl_radj", "Planet Radius", "Jupiter Radius", np.float]
        self._columns[25] = ["pl_radjerr1", "Planet Radius Upper Unc.", "Jupiter Radius", np.float]
        self._columns[26] = ["pl_radjerr2", "Planet Radius Lower Unc.", "Jupiter Radius", np.float]
        self._columns[27] = ["pl_radjlim", "Planet Radius Limit Flag", "", np.float]
        self._columns[28] = ["pl_bmasse", "Planet Mass or Mass*sin(i)", "Earth Mass", np.float]
        self._columns[29] = ["pl_bmasseerr1", "Planet Mass or Mass*sin(i) Upper Unc.", "Earth Mass", np.float]
        self._columns[30] = ["pl_bmasseerr2", "Planet Mass or Mass*sin(i) Lower Unc.", "Earth Mass", np.float]
        self._columns[31] = ["pl_bmasselim", "Planet Mass or Mass*sin(i) Limit Flag", "Earth Mass", np.float]
        self._columns[32] = ["pl_bmassj", "Planet Mass or Mass*sin(i)", "Jupiter Mass", np.float]
        self._columns[33] = ["pl_bmassjerr1", "Planet Mass or Mass*sin(i) Upper Unc.", "Jupiter Mass", np.float]
        self._columns[34] = ["pl_bmassjerr2", "Planet Mass or Mass*sin(i) Lower Unc.", "Jupiter Mass", np.float]
        self._columns[35] = ["pl_bmassjlim", "Planet Mass or Mass*sin(i) Limit Flag", "Jupiter Mass", np.float]
        self._columns[36] = ["pl_bmassprov", "Planet Mass or Mass*sin(i) Provenance", "", "U30"]
        self._columns[37] = ["pl_orbeccen", "Eccentricity", "", np.float]
        self._columns[38] = ["pl_orbeccenerr1", "Eccentricity Upper Unc.", "", np.float]
        self._columns[39] = ["pl_orbeccenerr2", "Eccentricity Lower Unc.", "", np.float]
        self._columns[40] = ["pl_orbeccenlim", "Eccentricity Limit Flag", "", np.float]
        self._columns[41] = ["pl_insol", "Insolation Flux", "Earth Flux", np.float]
        self._columns[42] = ["pl_insolerr1", "Insolation Flux Upper Unc.", "Earth Flux", np.float]
        self._columns[43] = ["pl_insolerr2", "Insolation Flux Lower Unc.", "Earth Flux", np.float]
        self._columns[44] = ["pl_insollim", "Insolation Flux Limit Flag", "", np.float]
        self._columns[45] = ["pl_eqt", "Equilibrium Temperature", "K", np.float]
        self._columns[46] = ["pl_eqterr1", "Equilibrium Temperature Upper Unc.", "K", np.float]
        self._columns[47] = ["pl_eqterr2", "Equilibrium Temperature Lower Unc.", "K", np.float]
        self._columns[48] = ["pl_eqtlim", "Equilibrium Temperature Limit Flag", "", np.float]
        self._columns[49] = ["ttv_flag", "Data show Transit Timing Variations", "", np.float]
        self._columns[50] = ["st_refname", "Stellar Parameter Reference", "", "U200"]
        self._columns[51] = ["st_spectype", "Spectral Type", "", "U10"]
        self._columns[52] = ["st_teff", "Stellar Effective Temperature", "K", np.float]
        self._columns[53] = ["st_tefferr1", "Stellar Effective Temperature Upper Unc.", "K", np.float]
        self._columns[54] = ["st_tefferr2", "Stellar Effective Temperature Lower Unc.", "K", np.float]
        self._columns[55] = ["st_tefflim", "Stellar Effective Temperature Limit Flag", "", np.float]
        self._columns[56] = ["st_rad", "Stellar Radius", "Solar Radius", np.float]
        self._columns[57] = ["st_raderr1", "Stellar Radius Upper Unc.", "Solar Radius", np.float]
        self._columns[58] = ["st_raderr2", "Stellar Radius Lower Unc.", "Solar Radius", np.float]
        self._columns[59] = ["st_radlim", "Stellar Radius Limit Flag", "", np.float]
        self._columns[60] = ["st_mass", "Stellar Mass", "Solar mass", np.float]
        self._columns[61] = ["st_masserr1", "Stellar Mass Upper Unc.", "Solar mass", np.float]
        self._columns[62] = ["st_masserr2", "Stellar Mass Lower Unc.", "Solar mass", np.float]
        self._columns[63] = ["st_masslim", "Stellar Mass Limit Flag", "", np.float]
        self._columns[64] = ["st_met", "Stellar Metallicity", "dex", np.float]
        self._columns[65] = ["st_meterr1", "Stellar Metallicity Upper Unc.", "dex", np.float]
        self._columns[66] = ["st_meterr2", "Stellar Metallicity Lower Unc.", "dex", np.float]
        self._columns[67] = ["st_metlim", "Stellar Metallicity Limit Flag", "", np.float]
        self._columns[68] = ["st_metratio", "Stellar Metallicity Ratio", "", "U10"]
        self._columns[69] = ["st_logg", "Stellar Surface Gravity", "log10(cm/s**2)", np.float]
        self._columns[70] = ["st_loggerr1", "Stellar Surface Gravity Upper Unc.", "log10(cm/s**2)", np.float]
        self._columns[71] = ["st_loggerr2", "Stellar Surface Gravity Lower Unc.", "log10(cm/s**2)", np.float]
        self._columns[72] = ["st_logglim", "Stellar Surface Gravity Limit Flag", "", np.float]
        self._columns[73] = ["sy_refname", "System Parameter Reference", "", "U200"]
        self._columns[74] = ["rastr", "RA", "sexagesimal", "U13"]
        self._columns[75] = ["ra", "RA", "deg", np.float]
        self._columns[76] = ["decstr", "Dec", "sexagesimal", "U13"]
        self._columns[77] = ["dec", "Dec", "deg", np.float]
        self._columns[78] = ["sy_dist", "Distance", "pc", np.float]
        self._columns[79] = ["sy_disterr1", "Distance Upper Unc", "pc", np.float]
        self._columns[80] = ["sy_disterr2", "Distance Lower Unc", "pc", np.float]
        self._columns[81] = ["sy_vmag", "V (Johnson) Magnitude", "mag", np.float]
        self._columns[82] = ["sy_vmagerr1", "V (Johnson) Magnitude Upper Unc", "mag", np.float]
        self._columns[83] = ["sy_vmagerr2", "V (Johnson) Magnitude Lower Unc", "mag", np.float]
        self._columns[84] = ["sy_kmag", "Ks (2MASS) Magnitude", "mag", np.float]
        self._columns[85] = ["sy_kmagerr1", "Ks (2MASS) Magnitude Upper Unc", "mag", np.float]
        self._columns[86] = ["sy_kmagerr2", "Ks (2MASS) Magnitude Lower Unc", "mag", np.float]
        self._columns[87] = ["sy_gaiamag", "Gaia Magnitude", "mag", np.float]
        self._columns[88] = ["sy_gaiamagerr1", "Gaia Magnitude Upper Unc", "mag", np.float]
        self._columns[89] = ["sy_gaiamagerr2", "Gaia Magnitude Lower Unc", "mag", np.float]
        self._columns[90] = ["rowupdate", "Date of Last Update", "", "U10"]
        self._columns[91] = ["pl_pubdate", "Planetary Parameter Reference Publication Date", "", "U10"]
        self._columns[92] = ["releasedate", "Release Date", "", "U10"]
        
        # Check whether data file exists
        self._fs = pp.PyAFS()

        if self.needsUpdate():
            # Data needs update
            print("Downloading exoplanet data from NASA exoplanet archive")
            self._update(self._downloadData)
            print("Saved data to file: ",
                  self.dataFileName, " in data directory,")
            print("  which has been configured as: ", self._fs.dpath)
            print("By default, the data will be downloaded anew every 7 days.")
            print("You can use the `changeDownloadCycle` to change this behavior.")
        self._readData()

    def downloadData(self):
        """
        Trigger download of data.
        """
        self._update(self._downloadData)

    def _readData(self):
        """
        Read the data from local file into numpy recarray.
        """
        r = list(csv.DictReader(self._fs.requestFile(
            self.dataFileName, 'rt', gzip.open), delimiter=','))
        dtype = [(self._columns[x+1][0], self._columns[x+1][3])
                for x in range(len(self._columns))]

        self.data = np.recarray((len(r) + 1,), dtype=dtype)
        for i, x in enumerate(r):
            for k, v in six.iteritems(x):
                self.data[k][i] = v if len(v) > 0 else None

    def availableColumns(self, verbose=True):
        """
        Shows a list of available data columns.

        Parameters
        ----------
        verbose : boolean, optional
            If True (default), prints information
            to screen.

        Returns
        -------
        columns : list of strings
            The names of the available data columns.
        """
        if self.data is None:
            return None
        if verbose:
            print("{0:12s}  {1:35s}  {2:5s}".format(
                "Column", "Description", "Unit"))
            print("-" * 56)
            for v in six.itervalues(self._columns):
                print("{0:12s}  {1:35s}  {2:5s}".format(*v))

        return [self._columns[x+1][0] for x in range(len(self._columns))]

    def selectByPlanetName(self, planetName, caseSensitive=False):
        """
        Get entry by planet name.

        Parameters
        ----------
        planetName : string
            The name of the planet (includes planet letter,
            e.g., "corot-2 b"
        caseSensitive : boolean, optional
            If False (default), the search will be case-insensitive.

        Returns
        -------
        Data entry : dictionary
            A dictionary with a key for every data column holding
            the associated value from the data table.
        """
        if hasattr(self.data.pl_name[0], "decode"):
            pnames = [name.decode("utf8") for name in self.data.pl_name]
        else:
            pnames = self.data.pl_name
        r = pyaC.fuzzyMatch(planetName, pnames,
                            caseSensitive=caseSensitive, raises=True)
        result = {}
        for c in six.itervalues(self._columns):
            result[c[0]] = self.data[c[0]][r["index"]]
        return result

    def getAllData(self):
        """
        Get all available data.

        Returns
        -------
        Data : numpy recarray
            All data stored in the table as a numpy recarray.
        """
        return self.data.copy()
