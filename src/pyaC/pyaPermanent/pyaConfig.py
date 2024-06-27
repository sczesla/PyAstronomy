from __future__ import print_function
import os
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves.configparser as ConfigParser
import six.moves as smo

class PyARC(object):
    """
    """
    
    def __init__(self):
        self.rc = ConfigParser.RawConfigParser()
        self.rc.add_section("pyaRC")
        # By default, keep online access on
        self.goOnline()
    
    def takeOffline(self):
        self.rc.set("pyaRC", "onlineKillSwitch", 1)
    def goOnline(self):
        self.rc.set("pyaRC", "onlineKillSwitch", 0)
    def supposedOnline(self):
        return (not self.rc["pyaRC"]["onlineKillSwitch"])
       
pyaRC = PyARC() 


class PyAConfig(object):
    """
    Provide access to permanent PyA configuration.

    Permanent configuration and data created by/for
    PyAstronomy are stored in a data directory appointed
    to PyAstronomy. The directory name is stored in a file named
    ".pyaConfigWhere", which is located in the user's home directory
    (defined by environment variable `home`). This file contains only
    a single line, which is exactly the name of the PyA data directory.

    In the PyA data directory, this class creates a stub configuration
    file named "pyaConfig.cfg", if the file does not already exist. 

    Attributes
    ----------
    dpath : string
        The name of PyA's data directory (defined by the content of
        the file ".pyaConfigWhere" in the home directory).
        None if no valid directory could be found or created.           .
    configWhere : string
        The full name of the ".pyaConfigWhere" file.
    """

    _rootConfig = None

    def __createConfigStub(self):
        """
        Creates a stub file for configuration if not yet existing.
        """
        if not os.path.isfile(os.path.join(self.dpath, "pyaConfig.cfg")):
            config = ConfigParser.RawConfigParser()
            with open(os.path.join(self.dpath, "pyaConfig.cfg"), 'wt') as configfile:
                config.write(configfile)
        else:
            return

    def __checkDataDir(self):
        """
        Check whether a valid data directory is present.

        This method throws an exception if no valid data directory
        could be found. The check is performed simply by checking
        the `dpath` attribute.
        """
        if self.dpath is None:
            raise(PE.PyAValError("No valid data and configuration directory for PyA could be found.",
                                 solution="Check/Repeat initialization (maybe remove '.pyaConfigWhere' in home directory)."))

    def getDataRoot(self):
        """
        Access the root path to PyA's permanent data and configuration.

        Returns
        -------
        Path : string
            The root path as defined in the file ".pyaConfigWhere" in the
            home directory.
        """
        self.__checkDataDir()
        return self.dpath

    def accessRootConfig(self):
        """
        Provides access to PyA's root configuration file.

        Returns
        -------
        Configuration : ConfigParser object
            The content of the root configuration file "pyaConfig.cfg".
        """
        return PyAConfig._rootConfig

    def _ynQuestion(self, message):
        """
        Ask yes/no question only accepting definite answer

        Parameters
        ----------
        message : string
            The question being asked.

        Returns
        -------
        Answer : boolean
            True if the answer was yes and False if it was no.
        """
        while True:
            yn = smo.input(message)
            if (yn.lower() == "y" or yn.lower() == "yes"):
                return True
            elif (yn.lower() == "n" or yn.lower() == "no"):
                return False
            else:
                print()
                print("  What do you by '" + str(yn) +
                      "'? Please provide a valid input (y/n).")
                print()

    def _locateHomeDirectory(self):
        """ Try to locate HOME directory and return if possible (else return None) """
        hd = None
        hd = os.getenv("HOME")
        if hd is None:
            # Try alternative solution
            hd = os.path.expanduser("~")
        return hd
        
    def __init__(self):
        self.dpath = None
        # Try to locate home directory via environment variable
        self.homeDir = self._locateHomeDirectory()
        if self.homeDir is None:
            PE.warn(PE.PyAValError("Could not find a home directory. Data directory cannot be set up.",
                                   solution="Set 'HOME' environment variable."))
            return

        # Hard-coded name of file, which will only contain path to "real" configuration.
        # This is done, because at one point, information must be stored, which can be
        # found without any extra information.
        self.configWhere = os.path.join(self.homeDir, ".pyaConfigWhere")
        if not os.path.isfile(self.configWhere):
            # Ask user whether she/he wants to configure data path
            print(" -------------- Configure PyA's data path ------------------------")
            print("  Why do you get this message?")
            print("    Most probably PyA tries to save permanent data on your")
            print("    system for the first time. This may be the case when, e.g.,")
            print("    a table is downloaded and saved. These data are stored")
            print("    under PyA's 'data path', i.e., a regular directory on")
            print("    your disk, which PyA can access.")
            print("  Is anything else saved to the disk?")
            print("    The location of the data path")
            print("    directory will be written to the file '.pyaConfigWhere' in")
            print("    your home directory, so that you need not provide it")
            print("    again.")
            print("  Can I delete it?")
            print("    Yes, you are free to delete anything of this at any time.")
            print("  What is a valid data path?")
            print("    You may provide any existing or non-existing path, although")
            print("    it is strongly encouraged to use a fresh directory to avoid any")
            print("    confusion. The given path needs to be absolute, to uniquely")
            print("    identify it.")
            print(" -----------------------------------------------------------------")
            print("")
            if self._ynQuestion("Configure PyA's data path now (y/n)? "):
                print()
                print("  Configure data path now.")
                print()
            else:
                print()
                print("  Configure data path later.")
                print()
                return

            try:
                while True:
                    print(
                        "Please provide a directory where PyA can store data (may already exist):")
                    print("Press enter to accept default; use 'exit' to abort.")
                    suggestion = os.path.join(self.homeDir, "PyAData")
                    dpath = smo.input(
                        "  Path (default = " + suggestion + "): ")
                    if dpath.lower().strip() == "exit":
                        print("  Process aborted. No data path configured.")
                        dpath = None
                        break
                    if dpath == "":
                        dpath = suggestion
                    # Check whether data path is an absolute path
                    if not os.path.isabs(dpath):
                        print("The path you specified (" + dpath +
                              ") is not absolute, but it needs to be.")
                        dpatha = os.path.abspath(dpath)
                        print("Did you intend to use this path:")
                        print("  ", dpatha)
                        if self._ynQuestion("(y/n) ?"):
                            dpath = dpatha
                            break
                        else:
                            # Retry to ask for path
                            continue
                    if not os.path.isdir(dpath):
                        # Try to create the directory
                        os.makedirs(dpath)
                        break
                    else:
                        print("Directory '" + dpath + "' exists.")
                        if self._ynQuestion("Use as PyA data directory (you may want to use, e.g., a subdirectory instead) (y/n) "):
                            if not os.access(dpath, os.W_OK):
                                # Exists and shall be used. Is it writable?
                                print("")
                                print(
                                    "The directory '" + dpath + "' is not writable! Therefore, it cannot be used.")
                                print("")
                                continue
                        break

                if dpath is None:
                    # Process has been aborted
                    print("No valid data path configured.")
                    return

                self.dpath = dpath

                # There is not yet a file '.pyaConfigWhere', which stores the place to look
                # for the real configure-file and data files.
                with open(self.configWhere, 'wt') as f:
                    f.write(self.dpath)

                # All could be created appropriately
                print("PyA data path configured successfully. Using path: ")
                print("  " + self.dpath)

            except Exception as e:
                PE.warn(PE.PyAValError("The directory: '" + str(dpath) + "' " +
                                       " could not be created. Data directory cannot be set up.\n" +
                                       "Error message: " + str(e)))
                os.remove(self.configWhere)
                return

        else:
            # There is a .pyaConfigWhere file.
            try:
                with open(self.configWhere) as cfi:
                    self.dpath = cfi.readline()
            except Exception as e:
                PE.warn(PE.PyAValError("The file " + self.configWhere +
                                       " exists, but could not be opened for reading.",
                                       solution="Check permissions of the file.",
                                       addInfo="Error message: " + str(e)))
                self.dpath = None
            try:
                self.dpath = os.path.realpath(self.dpath)
            except Exception as e:
                PE.warn(PE.PyAValError("Obtained the path '" + self.dpath + "' from .pyaConfigWhere, but" +
                                       "could not expand soft-links etc. (using os.path.realpath).",
                                       solution="Check the path written to .pyaConfigWhere in home directory.",
                                       addInfo="Error message: " + str(e)))
                self.dpath = None
            if not os.path.isdir(self.dpath):
                PE.warn(PE.PyAValError("The directory " + self.dpath + "' " +
                                       "was specified as data path in file '" + self.configWhere + "', " +
                                       "but it does no longer exist! No data can be stored.",
                                       solution=["Delete file '" + self.configWhere + "' to allow reconfiguration.",
                                                 "Modify the file content by giving the (modified) location " +
                                                 "of the data directory."]))
                self.dpath = None
                return
        # Create cfg file unless it already exists
        self.__createConfigStub()
        # Open the "root" configuration file for later access
        if PyAConfig._rootConfig is None:
            PyAConfig._rootConfig = ConfigParser.RawConfigParser()
            PyAConfig._rootConfig.read(
                os.path.join(self.dpath, "pyaConfig.cfg"))

    def saveConfigToFile(self):
        """
        Save current state of configuration to file 'pyaConfig.cfg'.
        """
        with open(os.path.join(self.dpath, "pyaConfig.cfg"), 'wt') as configfile:
            PyAConfig._rootConfig.write(configfile)

    def set(self, section, option, value):
        """
        Add an entry to the global configuration file.

        If the specified section does not already exist, it is created.
        """
        if not PyAConfig._rootConfig.has_section(section):
            PyAConfig._rootConfig.add_section(section)
        PyAConfig._rootConfig.set(section, option, value)
        self.saveConfigToFile()

    def remove_option(self, section, option):
        """
        Remove option
        """
        PyAConfig._rootConfig.remove_option(section, option)
        self.saveConfigToFile()

    def remove_section(self, section):
        """
        Remove section
        """
        PyAConfig._rootConfig.remove_section(section)
        self.saveConfigToFile()



