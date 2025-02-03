

#
# Bill Meyer - 7/24/2019
# meyer8@llnl.gov
#
#
import sys


def sources():
   """ 
       This routine simply dumps all the modules as a sorted list along
       with the "__file__" attribute. Prints "unknown" for those that don't
       have this attribute. This is just for debug. When a user reports a
       problem it may be useful to have them run this to make sure they are
       getting the modules from expected places.
   """
   for m in sorted(sys.modules.keys()):
       t = sys.modules[m]
       try:
          f = t.__file__
       except:
          f = 'unknown'
       finally:
          print(m,'\t--\t',f)

