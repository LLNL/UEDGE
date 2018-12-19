# Import all required libraries

import numpy
import pylab
import string

#
# The purpose of this class is to read a gridue file and present the
# data in a format friendly for use elsewhere
#
class UeReadGrid:

#
# Constructor
#
    def __init__(self):
        pass

    def readData(self, fname):
        fh = open(fname, 'r')

        lns = fh.readlines()

        # Read the header information including metadata and grid shape
        ln1 = lns.pop(0).split()

        self.xxr = ln1[0]
        self.yyr = ln1[1]

        self.xsp1 = ln1[2]
        self.xsp2 = ln1[3]
        self.ysp = ln1[4]

        lns.pop(0)

        # Reshape the grid data to be linear
        data = numpy.zeros( (lns.__len__() * 3), numpy.float )
        print data.shape
        for i in range(0, lns.__len__()-1):
            ll = lns[i].split()
            print ll
            data[3*i  ] = string.atof( ll[0] )
            data[3*i+1] = string.atof( ll[1] )
            data[3*i+2] = string.atof( ll[2] )

        rml = 0
        rmh = (self.xxr+2) * (self.yyr+2) * 5
        self.rm = data[rml:rmh].reshape( (self.xxr+2, self.yyr+2, 5) )

        zml = rmh
        zmh = rmh + (self.xxr+2) * (self.yyr+2) * 5
        self.zm = data[zml:zmh].reshape( (self.xxr+2, self.yyr+2, 5) )

uegrid = UeReadGrid()
uegrid.readData('gridue')
pylab.plot(uegrid.rm[:,:,0], uegrid.zm[:,:,0], 'b+')
