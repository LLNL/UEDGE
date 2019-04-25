# Show the range of a NymPy array
#
# Usage example:
# >>> ShowRange('Te [ev]',  bbb.te/ev)
#
# First coding: MVU, 26-jul-17
#=========================================================#

def ShowRange(name,v):
    print("%s: [%8.2e,%8.2e]" % (name, v.min(),v.max()))
#=========================================================#
