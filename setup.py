#!/usr/bin/env python
# To use:
#       python setup.py install
#
import sys
import os
import os.path
import string
import site
from Forthon.compilers import FCompiler
import getopt
import logging

version='8.0.0'

try:
    os.environ['PATH'] += os.pathsep + site.USER_BASE + '/bin'
    import setuptools
    import distutils
    from distutils.core import setup
    from distutils.core import Extension
    from distutils.dist import Distribution
    from distutils.command.build import build
    from distutils.command.install import install
    from subprocess import call
    import numpy
except:
    raise SystemExit("Distutils problem")

optlist, args = getopt.getopt(sys.argv[1:], 'gt:F:', ['parallel', 'petsc', 'omp'])
machine = sys.platform
debug = 0
fcomp = None
parallel = 0
petsc = 0

for o in optlist:
    if o[0] == '-g':
        debug = 1
    elif o[0] == '-t':
        machine = o[1]
    elif o[0] == '-F':
        fcomp = o[1]
    elif o[0] == '--parallel':
        parallel = 1
    elif o[0] == '--petsc':
        petsc = 1
    elif o[0] == '--omp':
        os.putenv("OMP","1")
        


if petsc == 1 and os.getenv('PETSC_DIR') == None:
    raise SystemExit("PETSc requested but PETSC_DIR not set")
if os.getenv('PETSC_DIR') != None:
    petsc = 1
if petsc == 1 and os.getenv('PETSC_ARCH') == None:
    raise SystemExit("PETSc requested but PETSC_ARCH not set")


sys.argv = ['setup2.py']+args
fcompiler = FCompiler(machine=machine,
                      debug=debug,
                      fcompname=fcomp)


class uedgeInstall(build):
    def run(self):
       install.run(self)
       logging.basicConfig(stream=sys.stderr,level=logging.INFO)
       log = logging.getLogger()
       log.info("test")
class uedgeBuild(build):
    def run(self):
        # with python2 everything is put into a single uedgeC.so file
        if sys.hexversion < 0x03000000:
            raise SystemExit("Python versions < 3 not supported")
        else:
            if petsc == 0:
                call(['make', '-f','Makefile.Forthon'])
            else:
                call(['make', '-f', 'Makefile.PETSc'])
            build.run(self)


class uedgeClean(build):
    def run(self):
        if sys.hexversion < 0x03000000:
            raise SystemExit("Python versions < 3 not supported")
        else:
            if petsc == 0:
                call(['make', '-f', 'Makefile.Forthon', 'clean'])
            else:
                call(['make', '-f', 'Makefile.PETSc', 'clean'])

uedgepkgs = ['aph', 'api', 'bbb', 'com', 'flx', 'grd', 'svr', 'wdf', 'ncl']


def makeobjects(pkg):
    return [pkg+'_p.o', pkg+'pymodule.o']


uedgeobjects = []

# add here any extra dot o files other than pkg.o, pkg_p.o


if sys.hexversion < 0x03000000:
    raise SystemExit("Python versions < 3 not supported")
else:
    dummydist = Distribution()
    dummydist.parse_command_line()
    dummybuild = dummydist.get_command_obj('build')
    dummybuild.finalize_options()
    builddir = dummybuild.build_temp

uedgeobjects = map(lambda p: os.path.join(builddir, p), uedgeobjects)

if os.getenv('PACT_DIR') != None:
    library_dirs = fcompiler.libdirs + [
        os.path.join(os.getenv('PACT_DIR'), 'lib')]
    libraries = ['pdb', 'pml', 'score', 'blas', 'm'] + fcompiler.libs
else:
    library_dirs = fcompiler.libdirs
    libraries = fcompiler.libs

if petsc:
    # PETSC_DIR = '/homes/mccomic/petsc-uedge'
    # PETSC_ARCH = 'linux-uedge'
    PETSC_DIR = os.getenv('PETSC_DIR')
    PETSC_ARCH = os.getenv('PETSC_ARCH')
    library_dirs = fcompiler.libdirs + \
        [os.path.join(PETSC_DIR, PETSC_ARCH, 'lib')]
    libraries = ['petscts', 'petscsnes', 'petscksp', 'petscdm', 'petscmat',
                 'petscvec', 'petsc', 'HYPRE', 'mpich', 'lapack', 'blas', 'X11',
                 'pthread', 'rt', 'stdc++', 'm'] + fcompiler.libs
    libraries = ['petsc'] + fcompiler.libs

if parallel:
    library_dirs = fcompiler.libdirs + ['/usr/lpp/ppe.poe/lib']
    libraries = fcompiler.libs + ['mpi']
    # uedgeobjects = uedgeobjects + ['/usr/local/mpi/ifc_farg.o']

with open('pyscripts/__version__.py','w') as ff:
    ff.write("__version__ = '%s'\n"%version)
with open('pyscripts/__src__.py','w') as ff:
    ff.write("__src__ = '%s'\n"%os.getcwd())

define_macros=[("WITH_NUMERIC", "0"),
               ("FORTHON_PKGNAME", '\"uedgeC\"'),
               ("FORTHON","1")]

# check for readline
rlncom = "echo \"int main(){}\" | gcc -x c -lreadline - "
rln = os.system(rlncom)
if rln == 0: 
   define_macros = define_macros + [("HAS_READLINE","1")]
   os.environ["READLINE"] = "-l readline"
   libraries = ['readline'] + libraries


setup(name="uedge",
      version=version,
      author='Tom Rognlien',
      author_email="trognlien@llnl.gov",
      maintainer='Bill Meyer',
      maintainer_email='meyer8@llnl.gov',
      description="2D Fluid simulation of plasma and neutrals in magnetic fusion devices",
      platforms="Unix, Windows (cygwin), Mac OSX",
      packages=['uedge'],
      package_dir={'uedge': 'pyscripts'},
      # include_package_data=True,
      scripts=['pyscripts/pdb2hdf5', 'pyscripts/bas2py', 'pyscripts/hdf52pdb'],
      ext_modules=[Extension('uedge.uedgeC',
                             ['uedgeC_Forthon.c',
                              os.path.join(builddir, 'Forthon.c'),
                              'com/handlers.c', 'com/vector.c','bbb/exmain.c'],
                             include_dirs=[builddir, numpy.get_include()],
                             library_dirs=library_dirs,
                             libraries=libraries,
                             define_macros=define_macros,
                             extra_objects=uedgeobjects,
                             extra_link_args=['-g','-DFORTHON'] +
                             fcompiler.extra_link_args,
                             extra_compile_args=fcompiler.extra_compile_args
                             )],

      cmdclass={'build': uedgeBuild, 'clean': uedgeClean},
      test_suite="pytests",
      install_requires=['forthon', 'mppl','easygui'],
      # note that include_dirs may have to be expanded in the line above
      classifiers=['Programming Language :: Python',
                   'Programming Language :: Python :: 3']

      )
