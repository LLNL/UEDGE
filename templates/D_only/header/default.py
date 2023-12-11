# This is the test header file
# It is copied into the start of every test case
# created from this subdirectory. It is the 
# best place to include lines common to
# all test cases, such as rate paths etc.

from uedge import *
from uedge.hdf5 import hdf5_restore

try:
    bbb.oldseec = 0
    bbb.jhswitch = 0
except:
    pass

