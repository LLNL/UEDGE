import unittest
import os
import numpy as np
import sys,getopt



def prhelp():
    print("""
   Usage: python test_slab.py [-r|--ref] [-h]       or
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
        This is run pre-test.
        """

    def test_one(self):
        """
        Test that uedge will import
        """
        global ftol
        try:
            import uedge as ue
            assert True
        except:
            assert False



if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hr',['ref'])
    except getopt.GetoptError:
        prhelp()
        sys.exit(2)

    for opt,arg in opts:
        if opt in ('-h'):
           prhelp()
           sys.exit(2)
        elif opt in ('-r','--ref'):
           """
           Do something that will store reference data
           """
           sys.exit(2)


    unittest.main()


