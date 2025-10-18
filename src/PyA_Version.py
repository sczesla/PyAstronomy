import sys
from pathlib import Path
    
try:
    if sys.version_info >= (3, 11):
        import tomllib
    else:
        import tomli as tomllib
except ImportError:
    tomllib = None


def _extract_version(ll):
    """ Finds and returns version tag using simplistic parsing """
    for l in ll:
        s = l.split("=")
        if s[0].strip() == "version":
            return s[1].strip()
    return "0.0.0"

def _get_version_from_toml():
    """
    Finds and robustly parses pyproject.toml to get the version.
    Returns "0.0.0-unknown" on failure.
    """
    pyproject_path = Path(__file__).resolve().parent.parent / "pyproject.toml"
    if not tomllib:
        # Try simple parse
        with open(pyproject_path, 'rt') as to:
            return _extract_version(to.readlines())
    try:
        with open(pyproject_path, "rb") as f:
            pyproject_data = tomllib.load(f)

        return pyproject_data["project"]["version"]

    except (FileNotFoundError, KeyError, tomllib.TOMLDecodeError) as e:
        print(e)
        return "0.0.0-toml-error"
    
    
# Try to get the metadata library
try:
    if sys.version_info >= (3, 8):
        from importlib import metadata
    else:
        # For older Python, you need the 'importlib-metadata' backport
        import importlib_metadata as metadata
except ImportError:
    # If the backport is not installed, we cannot use this method.
    metadata = None

__version__ = "0.0.0-unknown"

if metadata is not None:
    try:
        __version__ = metadata.version("PyAstronomy")
    except metadata.PackageNotFoundError:
        __version__ = _get_version_from_toml()


def PyA_Version():
    return __version__

