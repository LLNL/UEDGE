import unittest
import os
import numpy as np
import sys,getopt
from multiprocessing import Process
thisfile=os.path.realpath(__file__)
thispath=os.path.dirname(thisfile)
sys.path.insert(0,thispath)
sys.path.insert(0,os.path.dirname(os.path.dirname(thispath))+'/testscripts')
import uetests as uet

ftol = 1.e-9

name = 'd3d_snull_D+C_orthog'
restart = 'rd_d3dHmsCog_strahl_Rp95'

def prhelp():
    print("""
   Usage: python test_slab.py [-r|--ref] [-h|--help]       or
          pytest --forked test_slab.py
   
               -r|--ref  to produce reference files for future test runs

          """)

#
# Note the numbering of the tests. Purists will say that 
# tests should be order independent. Because Uedge variables
# are in shared objects it is naturally stateful. Because of
# this stateful feature the order needs to be controlled.
#



class TestRun(unittest.TestCase):
    def setUp(self):
        """
        This is run pre-test. Reads in restart data and 
        re-converges solution.
        """
        import uedge as ue
        os.chdir(thispath)

    def test_reconv(self):
        """
        Test that initial re-converged solution fnrm is low.
        """
        global ftol
        import uedge as ue
        uet.startup(restart,ftol=ftol) 
        uet.check_fnorm(name+' Reconverge','ref_reconv.h5',restart,doassert=True)

    def test_perturb(self):
        """
        Test that initial re-converged solution fnrm is low.
        """
        global ftol
        import uedge as ue
        uet.startup(restart,ftol=ftol) 
        uet.perturb(restart,ftol=ftol) 
        uet.steadystate()
        uet.check_nfe(name+' Perturb','ref_perturb.h5',restart,doassert=True)

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hr',['help,ref'])
    except getopt.GetoptError:
        prhelp()
        sys.exit(2)

    for opt,arg in opts:
        if opt in ('-h','--help'):
           prhelp()
           sys.exit(2)
        elif opt in ('-r','--ref'):
           kargs = {'ftol':ftol}
           p1 = Process(target=uet.makeref,args=('ref_reconv.h5',restart),kwargs=kargs)
           p1.start()
           p2 = Process(target=uet.makeperturb,args=('ref_perturb.h5',restart),kwargs=kargs)
           p2.start()
           p1.join()
           p2.join()
           sys.exit(2)


    unittest.main()


