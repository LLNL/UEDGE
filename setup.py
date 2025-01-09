#!/usr/bin/env python
# To use:
#       python setup.py install
#
import os
import os.path
import string
import site
from Forthon.compilers import FCompiler
import getopt
import logging
from subprocess import call
from sys import hexversion, argv, platform
import numpy
try:
    from setuptools import Extension, setup, Distribution
    from setuptools.command.build import build
except:
    from distutils.core import Extension, setup
    from distutils.dist import Distribution
    from distutils.command.build import build
# Check Python version
if hexversion < 0x03000000:
    raise SystemExit("Python versions < 3 not supported")

# UEDGE packages to be included in the build
uedgepkgs = ['aph', 'bbb', 'api', 'com', 'flx', 'grd', 'svr', 'wdf', 'ppp', 'ncl']
# *.o-files to be included in the build (in addition to pkg.o and pkg_p.o)
uedgeobjects = []


os.environ['PATH'] += os.pathsep + site.USER_BASE + '/bin'


arglist = {
        'debug': False,
        'serial': False,
        'checkbounds': False,
        'noclean': False,
#        'petsc': False,
}



optlist, args = getopt.getopt(argv[1:], 'gt:F:', list(arglist.keys()))
machine = platform
fcomp = None
parallel = 0

for o in optlist:
    if o[0] == '-g':
        arglist['debug'] = True
    elif o[0] == '-t':
        machine = o[1]
    elif o[0] == '-F':
        fcomp = o[1]
    elif o[0] == '--parallel':
        parallel = 1
    elif o[0].replace('--','') in arglist.keys():
        arglist[o[0].replace('--','')] = True

# Capture and extract PIP options passed using -C
for arg, _ in arglist.items():
    if "--{}".format(arg) in argv:
        argv.remove("--{}".format(arg))
        arglist[arg] = True

        
# Flags for makefile. Flags are easier to handle from setup.py and it prevents dealing with the makefile.)

FLAGS  = {
    # Specify FORTRAN compiler arguments
    'FARGS': [ 
        '-g -fmax-errors=15', 
        '-DFORTHON',
        '-cpp',
        '-Wconversion',
        '-fimplicit-none',
        '-fopenmp'
    ],
    # Specify C compiler arguments
    'CARGS': ['-fopenmp'],
    # Specify debugging flags
    'FARGSDEBUG': [
        '-fbacktrace',
        '-ffree-line-length-0', 
        '-fcheck=all',
        '-fbounds-check',
        '-ffpe-trap=invalid,overflow,underflow -finit-real=snan',
        '-Og'
    ],
    # Specify optional flags used if not debugging
    'FARGSOPT': ['-O3'],
    # Specify OMP args
    'OMPargs': ['--omp'],
}

fcompiler = FCompiler(machine=machine,
                      debug=int(arglist['debug']),
                      fcompname=fcomp)


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

 
if parallel:
    library_dirs = fcompiler.libdirs + ['/usr/lpp/ppe.poe/lib']
    libraries = fcompiler.libs + ['mpi']
    # uedgeobjects = uedgeobjects + ['/usr/local/mpi/ifc_farg.o']

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

# Modify Build to inject call to building UEDGE source
class uedgeBuild(build):
    def run(self):
        '''
        # Perform PETSc build [deprecated]
        if arglist['petsc']:
            if arglist['noclean'] is False:
                if call(['make', '-f', 'Makefile.PETSc', 'clean']) != 0:
                        raise SystemExit("Pre-PETSc UEDGE build clean failure")
            if os.getenv('PETSC_DIR') == None:
                raise SystemExit("PETSc requested but PETSC_DIR not set")
            if os.getenv('PETSC_ARCH') == None:
                raise SystemExit("PETSc requested but PETSC_ARCH not set")
            PETSC_DIR = os.getenv('PETSC_DIR')
            PETSC_ARCH = os.getenv('PETSC_ARCH')
            library_dirs = fcompiler.libdirs + [
                    os.path.join(PETSC_DIR, PETSC_ARCH, 'lib')
            ]
            libraries = [
                    'petscts', 'petscsnes', 'petscksp', 'petscdm', 'petscmat',
                    'petscvec', 'petsc', 'HYPRE', 'mpich', 'lapack', 'blas', 
                    'X11', 'pthread', 'rt', 'stdc++', 'm'] + fcompiler.libs
            libraries = ['petsc'] + fcompiler.libs
            if call(['make', '-f', 'Makefile.PETSc']) !=0:
                raise SystemExit("PETSc build failure")
            build.run(self)
            return
        '''
        # If serial buid requested, remove parallel flags
        if arglist['serial']:
            FLAGS['OMPargs'].remove('--omp')
            FLAGS['FARGS'].remove('-fopenmp')
            FLAGS['CARGS'].remove('-fopenmp')
        # Add boundary checking to non-debug build
        if arglist['checkbounds']:
            FLAGS['FARGSOPT'] = FLAGS['FARGSOPT']  + [
                '-fbounds-check',
                '-fcheck=all',
            ]
        # Check whether to add debug flags
        if arglist['debug']:
            FLAGS['FARGS'] = FLAGS['FARGS'] + FLAGS['FARGSDEBUG']
        else:
            FLAGS['FARGS'] = FLAGS['FARGS'] + FLAGS['FARGSOPT']
        # Clean build directory unless flagged not to
        if arglist['noclean'] is False:
            status = call(['make', '-f', 'Makefile.Forthon', 'clean'])
            if status != 0: 
                raise SystemExit("Pre-UEDGE build clean failure")

        # Compile arguments into compiler flags    
        OMPFLAGS='OMPFLAGS = {}'.format(' '.join(FLAGS['OMPargs']))
        COMPILEFLAGS ='DEBUG = -v --fargs "{}"'.format(' '.join(FLAGS['FARGS']))
        if FLAGS['CARGS']!=[]:
            COMPILEFLAGS = COMPILEFLAGS+' --cargs="{}"'.format(' '.join(FLAGS['CARGS']))
        # Call compiler
        if call(['make',COMPILEFLAGS,OMPFLAGS,'-f','Makefile.Forthon']) != 0:
            raise SystemExit("UEDGE build failure")
        build.run(self)

# Modify clean to use make clean
class uedgeClean(build):
    def run(self):
#        if call(['make', '-f', 'Makefile.'+ "PETSc"*(arglist['petsc']) + "Forthon"*(not arglist['petsc']), 'clean']) != 0:
        if call(['make', '-f', 'Makefile.Forthon', 'clean']) != 0:
            raise SystemExit("Clean failure")


setup(name="uedge",
      ext_modules=[Extension('uedge.uedgeC',
                             ['uedgeC_Forthon.c',
                              os.path.join(builddir, 'Forthon.c'),
                              'com/handlers.c', 'com/vector.c','bbb/exmain.c'],
                             include_dirs=[builddir, numpy.get_include()],
                             library_dirs=library_dirs,
                             libraries=libraries,
                             define_macros=define_macros,
                             extra_objects=uedgeobjects,
                             extra_link_args=FLAGS['CARGS']+['-g','-DFORTHON'] +
                             fcompiler.extra_link_args,
                             extra_compile_args=fcompiler.extra_compile_args
                             )],

      cmdclass={
            'build': uedgeBuild,
            'clean': uedgeClean},
)
