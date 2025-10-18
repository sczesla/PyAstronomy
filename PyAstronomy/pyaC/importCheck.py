from . import pyaErrors as PE
import importlib
import traceback


def pyaimportallfrom(mn, pak, globs, cra=False, rf=False, re=False):
    """
    Mimics the behavior of 'from X import *'.

    Get some help here.
    # https://stackoverflow.com/questions/21221358/python-how-to-import-all-methods-and-attributes-from-a-module-dynamically

    Parameters
    ----------
    mn : string
        Module name
    pak : string
        Name of the package to which the imported
        names are added.
    globs : dictionary
        The `globals()` dictionary of that package.
    cra : boolean, optional
        If True, symbol reassignments will be checked
        and some info printed.
    rf : boolean, optional
        If True, a warning will be given on import failure.
        Default is False.
    re : boolean, optional
        If False (default), no exception will be raised on
        import failure (note, also `rf` must be set True to
        make this happen).
        
    Returns
    -------
    status : boolean
        True if import was successful
    exception : string
        If an exception occurred, the text describing the exception. Otherwise
        an empty string
    traceback : list of strings
        If an exception occurred, the traceback. Otherwise an empty string.
    """
    try:
        n = importlib.import_module("." + mn, pak)
        try:
            # If __all__ is defined, use it!
            to_import = n.__all__
        except AttributeError:
            # Otherwise, import all names not starting with an underscore
            to_import = [
                name for name in n.__dict__ if not name.startswith('_')]
            if cra:
                # Check reassignment
                for name in to_import:
                    if name in globs:
                        print("Reassigning: ", name, " from ", n)

        globs.update({name: n.__dict__[name] for name in to_import})
    except Exception as e:
        if rf:
            f = PE.PyAImportFailure("Could not import module " + str(mn) + " to package " + str(pak) + ". " +
                                    "Received message: " + str(e))
            if not re:
                # Report failure without raising exception
                PE.warn(f)
            else:
                # Raise exception
                raise(f)
        return False, str(e), [l.rstrip("\n") for l in traceback.format_exc().splitlines(True)]
    return True, "", ""


class ImportCheck:

    def __init__(self, modules, required=None):
        """
        Checks whether individual modules can be imported.

        Parameters
        ----------
        modules : list of strings
            The names of the modules to be checked.
        required : list of strings, optional
            List of modules which are required. If not found,
            an import error is raised.             
        """
        if required is None:
            required = []
        # List of required modules that could not be imported,
        # list of all modules that could not be imported
        self.requiredFail = []
        self.fail = []

        self.check = {}
        self.versions = {}
        for module in modules:
            self.check[module] = True
            try:
                mi = __import__(module)
            except ImportError:
                self.check[module] = False
                self.fail.append(module)
            except PE.PyARequiredImport:
                self.check[module] = False
                self.fail.append(module)
            except Exception as e:
                # In fact, any error prevents import
                self.check[module] = False
                self.fail.append(module)

            self.versions[module] = "undefined"
            if self.check[module]:
                try:
                    self.versions[module] = mi.__version__
                except:
                    pass

            if (not self.check[module]) and (module in required):
                self.requiredFail.append(module)

        if len(self.requiredFail) > 0:
            ms = ', '.join(self.requiredFail)
            raise(PE.PyARequiredImport(
                "Could not import required module(s): " + ms, solution="Please install " + ms))
