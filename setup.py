# -*- coding: utf-8 -*-
from __future__ import print_function
try:
  from numpy.distutils.core import setup, Command
except(ImportError):
  print("Please install 'numpy' first.")
  
import glob
import sys
import re
sys.path.append("src")
from PyA_Version import PyA_Version

from distutils.extension import Extension as old_Extension

import re
cxx_ext_re = re.compile(r'.*[.](cpp|cxx|cc)\Z',re.I).match
fortran_pyf_ext_re = re.compile(r'.*[.](f90|f95|f77|for|ftn|f|pyf)\Z',re.I).match

class Extension(old_Extension):
    def __init__ (self, name, sources,
                  include_dirs=None,
                  define_macros=None,
                  undef_macros=None,
                  library_dirs=None,
                  libraries=None,
                  runtime_library_dirs=None,
                  extra_objects=None,
                  extra_compile_args=None,
                  extra_link_args=None,
                  export_symbols=None,
                  swig_opts=None,
                  depends=None,
                  language=None,
                  f2py_options=None,
                  module_dirs=None,
                  optional=False
                 ):
        old_Extension.__init__(self,name, [],
                               include_dirs,
                               define_macros,
                               undef_macros,
                               library_dirs,
                               libraries,
                               runtime_library_dirs,
                               extra_objects,
                               extra_compile_args,
                               extra_link_args,
                               export_symbols)
        # Avoid assert statements checking that sources contains strings:
        self.sources = sources

        # Python 2.4 distutils new features
        self.swig_opts = swig_opts or []

        # Python 2.3 distutils new features
        self.depends = depends or []
        self.language = language

        # numpy_distutils features
        self.f2py_options = f2py_options or []
        self.module_dirs = module_dirs or []

        return

    def has_cxx_sources(self):
        for source in self.sources:
            if cxx_ext_re(str(source)):
                return True
        return False

    def has_f2py_sources(self):
        for source in self.sources:
            if fortran_pyf_ext_re(source):
                return True
        return False










# By default build public distribution
sdist = False
withExt = False
versionAdd = ""

for s in sys.argv:
  if s == "sdist":
    sdist = True
  if s == "--with-ext":
    withExt = True

# Has to be removed. Otherwise distutils complain...
if withExt:
  sys.argv.pop(sys.argv.index("--with-ext"))

# The list of packages
packages = ['PyAstronomy', \
            'PyAstronomy.funcFit', \
            'PyAstronomy.funcFit.utils', \
            'PyAstronomy.pyasl', \
            'PyAstronomy.pyasl.asl', \
            'PyAstronomy.pyasl.asl.aslExt_1', \
            'PyAstronomy.pyasl.resBased', \
            'PyAstronomy.pyasl.phoenixUtils', \
            'PyAstronomy.modelSuite', \
            'PyAstronomy.modelSuite.XTran', \
            'PyAstronomy.modelSuite.XTran.palTrans', \
            'PyAstronomy.modelSuite.XTran.forTrans', \
            'PyAstronomy.pyTiming', \
            'PyAstronomy.pyTiming.pyPDM', \
            'PyAstronomy.pyTiming.pyPeriod', \
            'PyAstronomy.pyaC', \
            'PyAstronomy.pyaC.pyaErrors', \
            'PyAstronomy.pyaC.pyaPermanent', \
            'PyAstronomy.pyaC.mtools', \
            'PyAstronomy.constants', \
            'PyAstronomy.pyaGui']

# "doc/*/*.rst", "doc/*/*.png"],
package_data = {"PyAstronomy":["setup.cfg"],
                "PyAstronomy.modelSuite.XTran.palTrans":["ellint/makefile_template", "ellint/ell.cpp"], \
                "PyAstronomy.modelSuite.XTran.forTrans":["*.f", "*.pyf"], \
                "PyAstronomy.pyasl.asl":["testPro/*.pro", "*.cfg", "*.dat"], \
                "PyAstronomy.pyasl.resBased":["*.dat"], \
                "PyAstronomy.pyasl.asl.aslExt_1":["*.dat"], \
                "PyAstronomy.constants":["*.dat"], \
                "PyAstronomy.pyaC":["*.dat"], \
                "PyAstronomy.pyaC.mtools":["*.dat"]
                }

extOccultnl = Extension('PyAstronomy.modelSuite.XTran.forTrans.occultnl',
                    sources = ['src/modelSuite/XTran/forTrans/occultnl.pyf', \
                               'src/modelSuite/XTran/forTrans/occultnl.f'], optional=True)
extOccultquad = Extension('PyAstronomy.modelSuite.XTran.forTrans.occultquad',
                    sources = ['src/modelSuite/XTran/forTrans/occultquad.pyf', \
                               'src/modelSuite/XTran/forTrans/occultquad.f'], optional=True)

# If the --with-ext is not specified, these modules will not be built!
_ext_modules = [extOccultnl, extOccultquad]

if not withExt:
  ext_modules = []
else:
  ext_modules = _ext_modules


class WithExtCommand(Command):
    description = "Dummy command to allow with-ext option."
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        print("Version: ", PyA_Version())



if sdist:  
  
  manin = open("MANIFEST.in_template").readlines()
  # Add package data to Manifest.in
  for p in packages:
    if not p in package_data: continue
    manin.append("# Package data for: "+p+"\n")
    # Find the path from where to include
    path = p.replace('.', '/').replace("PyAstronomy", "src")
    for d in package_data[p]:
      manin.append("include "+path+"/"+d+"\n")
  
  # Add documentation in Manifest.in
  for p in packages:
    # Regular expression ignoring sub-packages.
    # These are automatically included in the
    # documentation.
    r = re.match("[^\.]+\.([^\.]*)$", p)
    if r is not None:
      manin.append("# Documentation for: "+p+"\n")
      manin.append("recursive-include src/doc/"+r.group(1)+"Doc *.rst *.png\n")
  
  open("MANIFEST.in", 'w').writelines(manin)



setup(cmdclass = {"with-ext":WithExtCommand},
      name='PyAstronomy',
      url="http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/index.html",
      description='A collection of astronomy related tools for Python.',
      version=PyA_Version()+versionAdd,
      packages = packages,
      ext_modules = ext_modules,
      package_dir = {'PyAstronomy':'src'},
      # Do not forget to give the ``correct'' name for the module! (here, e.g., PyAstronomy.funcFit)
      package_data = package_data,
      author='PyA group',
      author_email='stefan.czesla@hs.uni-hamburg.de',
      )

if not withExt:
  print("")
  print("")
  print("  USER INFO: External modules have not been built!")
  print("    The following modules have not been compiled:")
  for e in _ext_modules:
    print("    \""+e.name+"\", Sources: ", e.sources)
  print("  USE 'python setup.py --with-ext install' to build external modules")
  print("") 
  