#!/usr/bin/env python
# To use:
#       python setup.py install
#
import sys,os,os.path,string
from Forthon.compilers import FCompiler
import getopt

try:
    import distutils
    from distutils.core import setup, Extension
    from distutils.dist import Distribution
    from distutils.command.build import build
except:
    raise SystemExit, "Distutils problem"

optlist,args = getopt.getopt(sys.argv[1:],'gt:F:',['parallel','petsc'])
machine = sys.platform
debug   = 0
fcomp   = None
parallel = 0
petsc = 0

for o in optlist:
  if   o[0] == '-g': debug = 1
  elif o[0] == '-t': machine = o[1]
  elif o[0] == '-F': fcomp = o[1]
  elif o[0] == '--parallel': parallel = 1
  elif o[0] == '--petsc': petsc = 1  

if petsc: fcomp = 'g95'

sys.argv = ['setup.py']+args
fcompiler = FCompiler(machine=machine,
                      debug=debug,
                      fcompname=fcomp)

#dummydist = Distribution()
#dummybuild = build(dummydist)
#dummybuild.finalize_options()
#builddir = dummybuild.build_temp
builddir = 'build'

uedgepkgs = ['aph','api','bbb','com','flx','grd','svr','wdf']

def makeobjects(pkg):
  return [pkg+'_p.o',pkg+'pymodule.o']

uedgeobjects = []
for pkg in uedgepkgs:
  uedgeobjects = uedgeobjects + makeobjects(pkg)

# add here any extra dot o files other than pkg.o, pkg_p.o


uedgeobjects = uedgeobjects + ['aphrates.o','aphread.o',
                               'apifcn.o','apip93.o','apisorc.o',
                               'fimp.o','fmombal.o','inelrates.o',
                               'sputt.o','boundary.o','convert.o',
                               'geometry.o','griddubl.o','oderhs.o',
                               'odesetup.o','odesolve.o','potencur.o',
                               'turbulence.o','blasext.o','brent.o',
                               'comutil.o','misc.o','mnbrak.o',
                               'flxcomp.o','flxdriv.o','flxread.o',
                               'flxwrit.o','grdcomp.o','grddriv.o',
                               'grdinit.o','grdread.o','grdwrit.o',
                               'daspk.o','nksol.o','svrut1.o','svrut2.o',
                               'svrut3.o','svrut4.o','vodpk.o','uoa.o',
                               'dsum.o','dummy_py.o','error.o','getmsg.o',
                               'dmath.o','ssum.o','daux1.o','wdf.o']

if petsc:
  uedgeobjects = uedgeobjects +['petsc-uedge.o']

if parallel:
  # add extra dot o's needed if we're parallel
  uedgeobjects = uedgeobjects + []

uedgeobjects = map(lambda p:os.path.join(builddir,p),uedgeobjects)

library_dirs = fcompiler.libdirs
libraries = fcompiler.libs

if petsc:
  PETSC_DIR = '/homes/mccomic/petsc-uedge'
  PETSC_ARCH = 'linux-uedge'
  library_dirs = fcompiler.libdirs + [os.path.join(PETSC_DIR,PETSC_ARCH,'lib')]
  libraries =  ['petscts','petscsnes','petscksp','petscdm','petscmat',
                'petscvec','petsc','HYPRE','mpich','lapack','blas','X11',
                'pthread','rt','stdc++','m'] + fcompiler.libs

if parallel:
  library_dirs = fcompiler.libdirs + ['/usr/lpp/ppe.poe/lib']
  libraries = fcompiler.libs + ['mpi']
  #uedgeobjects = uedgeobjects + ['/usr/local/mpi/ifc_farg.o']

setup (name = "uedgeC",
       version = '3.0',
       author = 'Tom Rognlien',
       author_email = "trognlien@llnl.gov",
       description = "Combines uedge's packages into one",
       platforms = "Unix, Windows (cygwin), Mac OSX",
       ext_modules = [Extension('uedgeC',
                                ['uedgeC_Forthon.c',
                                 os.path.join(builddir,'Forthon.c'),
                                 'com/handlers.c','com/vector.c'],
                                include_dirs=[builddir],
                                library_dirs=library_dirs,
                                libraries=libraries,
                                define_macros=[("WITH_NUMERIC","1")],
                                extra_objects=uedgeobjects,
                                extra_link_args=['-g']+
                                             fcompiler.extra_link_args,
                                extra_compile_args=fcompiler.extra_compile_args
                               )]
       ## note that include_dirs may have to be expanded in the line above

       )
