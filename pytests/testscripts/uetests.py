import unittest
import os
import numpy as np
import sys,getopt
import uuid
from importlib import import_module as im
from multiprocessing import Process



def saverefs(filename):
    import uedge as ue
    import uedge.hdf5 as h5
    fnrm = ue.bbb.get_fnrm(ue.bbb.dtreal)
    nfe = ue.bbb.nfe[0][0]
    sclyl = ue.bbb.yldot * ue.bbb.sfscal
    nodeid = uuid.getnode()
    h5.hdf5_dump(filename,vars=['nodeid','fnrm','nfe','sclyl'],globals=locals())

   
    

def prhelp():
    print("""
   Usage: python test_slab.py [-r|--ref] [-h|--help]       or
          pytest --forked test_slab.py
   
               -r|--ref  to produce reference files for future test runs

          """)

def startup(case,ftol=1.e-9):
    """
       startup("<case>",ftol=1.e-9)
       Reconverge specified case. Done for every test.
    """
    import uedge as ue
    im(case)
    ue.bbb.ftol = ftol
    ue.bbb.exmain()

def perturb(case,ftol=1.e-9):
    """
       perturn("<case>",ftol=1.e-9)
       Perturb solution and Reconverge specified case.
    """
    import uedge as ue
    ue.bbb.ngs = ue.bbb.ngs*1.02
    ue.bbb.nis = ue.bbb.nis*1.02
    ue.bbb.phis = ue.bbb.phis*1.02
    ue.bbb.tes = ue.bbb.tes*1.02
    ue.bbb.tis = ue.bbb.tis*1.02
    ue.bbb.ups = ue.bbb.ups*1.02
    ue.bbb.ftol = ftol


def steadystate():
    """
       steadystate()
       Evolve current case to steady state with rdcontdt.
    """
    import uedge as ue
    import rdinitdt
    import rdcontdt

def makeperturb(filename,case,ftol=1.e-9):
    import uedge as ue
    startup(case,ftol=ftol)
    perturb(case,ftol=ftol)
    steadystate()
    saverefs(filename)
    

def makeref(filename,case,ftol=1.e-9):
    """
    Produce and save reference data for comparison in future test runs.
    """
    import uedge as ue
    startup(case,ftol=ftol)
    saverefs(filename)
    

identical=1.e-10
close=0.5


def check_fnorm(name,filename,case,doassert=True):
    import uedge as ue
    import uedge.hdf5 as h5
    ref={}
    h5.hdf5_restore_dump(filename,scope=ref)
    fnrm = ue.bbb.get_fnrm(ue.bbb.dtreal)
    nfe = ue.bbb.nfe[0][0]


    if np.isclose(fnrm,ref['fnrm'],atol=0.0,rtol=identical):
       #
       # If it makes it here then fnorm is basically identical.
       # Most likely running on the same os/compiler versions 
       # as when the reference files were produced.
       #
       if doassert:
          print()
          print(name,' fnorm identical.')
          assert True
       else:
          return True
    elif fnrm > ue.bbb.ftol:
       if doassert:
          print('Relative change in Fnorm too large. Threshold is ',close)
          print('fnrm: ',fnrm,'  ref: ',ref['fnrm'])
          print('   rel change: ',np.abs(fnrm - ref['fnrm'])/np.abs(ref['fnrm']))
          print('   abs change: ',np.abs(fnrm - ref['fnrm']))
          print()
       sclyl = ue.bbb.yldot * ue.bbb.sfscal
       rsclyl = ref['sclyl']
       iva = np.abs(rsclyl - sclyl) / (np.abs(rsclyl) + np.abs(sclyl) + 1e-20)
       ind = np.where(iva == np.max(iva))
       iv = ind[0][0]
       (ix,iy) = ue.bbb.igyl[iv,0:2]
       loc_troub_eqn = np.mod(iv,ue.bbb.numvar)+1
       numvar = ue.bbb.numvar
       if doassert:
          print("** Number of variables is:")
          print("numvar = ", numvar)
          print(" ")
          print("** Troublemaker equation is:")
          print("loc_troub_eqn = ",loc_troub_eqn)
          print(" ")
          print("** Troublemaker cell (ix,iy) is:")
          print(ue.bbb.igyl[iv,:])
          print(" ")
          print("** Timestep for troublemaker equation:")
          print(ue.bbb.dtuse[iv])
          print(" ")
          print("** yl for troublemaker equation:")
          print(ue.bbb.yl[iv])
          print(" ")
          assert False
       else:
          return False
    else:
       if doassert:
          assert True
       else:
          return True
       

       

def check_nfe(name,filename,case,doassert=None):
    import uedge as ue
    import uedge.hdf5 as h5
    ref={}
    h5.hdf5_restore_dump(filename,scope=ref)
    fnrm = ue.bbb.get_fnrm(ue.bbb.dtreal)
    nfe = ue.bbb.nfe[0][0]

    if np.isclose(fnrm,ref['fnrm'],atol=0.0,rtol=identical):
       #
       # If it makes it here then fnorm is basically identical.
       # Most likely running on the same os/compiler versions 
       # as when the reference files were produced.
       #
       if doassert:
          print()
          print(name,' fnorm identical.')
          assert True
       else:
          return True
    elif not np.isclose(nfe,ref['nfe'],atol=0.0,rtol=0.02):
       if doassert:
          print("The number of Krylov iterations for a 2% perturbation is  ",nfe/ref['nfe']," times the ref case")
          assert False
       else:
          return False
    else:
       if doassert:
          assert True
       else:
          return True
       
       

       

    

        
    
       

    

        
    
