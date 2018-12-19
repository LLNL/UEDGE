# Import all required libraries

import numpy
import pylab
import string

#
# This class will generate a lineout of the given data across the
# midplane starting at the magnetic axis and going all the way to the
# edge.  Eventually it should be capable of generating a lineout along
# any given ray, but currently only midplane is implemented
#
class UeGenLineout:

#
# Constructor
#
    def __init__(self):
        pass

#
# Use this to set the core region data and grids
#
# @param xx the XX array of the mapped grid
# @param yy the YY array of the mapped grid
# @param da the data at each x,y point described by xx and yy
#
    def setCoreData(self, xx, yy, da):
        self.eqXX = xx
        self.eqYY = yy
        self.eqDA = da

#
# Use this to set the interior edge region data and grids
#
# @param xx the XX array of the mapped grid centers
# @param yy the YY array of the mapped grid centers
# @param da the data at each x,y point described by xx and yy
#
    def setInteriorEdgeData(self, xx, yy, da):
        self.coXX = xx
        self.coYY = yy
        self.coDA = da

#
# Use this to set the scrape-off layer region data and grids
#
# @param xx the XX array of the mapped grid
# @param yy the YY array of the mapped grid
# @param da the data at each x,y point described by xx and yy
#
    def setSOLEdgeData(self, xx, yy, da):
        self.scXX = xx
        self.scYY = yy
        self.scDA = da

#
# Use this to get a lineout of the data input.  Call all the set
# functions before calling this
#
# @param fname file from which to read processed uedge grid points
#
# @return a tuple of numpy arrays representing the lineout data
#     1: The array of X coordinates of the data
#     2: The array of data representing the lineout
#
    def getLineout(self, fname):
        eqRetX = self.eqXX[0,:]
        eqRetD = self.eqDA[0,:]

        fh = open(fname, 'r')

        lns = fh.readlines()

        da = numpy.zeros((lns.__len__(), 5), numpy.float)
        for i in range(0, lns.__len__()):
            ln = lns[i].split()
            for j in range(0, ln.__len__()):
                da[i, j] = string.atof( ln[j] )

        n = 0
        for i in range(0, self.coXX.shape[1]):
            arrX = self.coXX[:,i]
            arrY = self.coYY[:,i]
            arrD = self.coDA[:,i]

            q = 0;
            for j in range (1, arrY.shape[0]):
                if ( ( arrX[j-1] > eqRetX[0] ) and
                     ( ( arrY[j-1] < da[n, 1] and
                         arrY[j]   > da[n, 1] ) or
                       ( arrY[j]   < da[n, 1] and
                         arrY[j-1] > da[n, 1] ) ) ):
                    q = j
                    break

            rat = (da[n, 1] - arrY[q-1]) / (arrY[q] - arrY[q-1])
            da[n, 4] = rat * arrD[q] + (1-rat)*arrD[q-1]

            n = n+1

        for i in range(0, self.scXX.shape[1]):
            arrX = self.scXX[:,i]
            arrY = self.scYY[:,i]
            arrD = self.scDA[:,i]

            q = 0;
            for j in range (1, arrY.shape[0]):
                if ( ( arrX[j-1] > eqRetX[0] ) and
                     ( ( arrY[j-1] < da[n, 1] and
                         arrY[j]   > da[n, 1] ) or
                       ( arrY[j]   < da[n, 1] and
                         arrY[j-1] > da[n, 1] ) ) ):
                    q = j
                    break

            rat = (da[n, 1] - arrY[q-1]) / (arrY[q] - arrY[q-1])
            da[n, 4] = rat * arrD[q] + (1-rat)*arrD[q-1]

            n = n+1

        retX = numpy.zeros( (eqRetX.shape[0] + da.shape[0] - 1), numpy.float )
        retD = numpy.zeros( (eqRetX.shape[0] + da.shape[0] - 1), numpy.float )

        retX[ :eqRetX.shape[0]-1 ] = eqRetX[:-1]
        retX[ eqRetX.shape[0]-1: ] = da[:,0]

        retD[ :eqRetX.shape[0]-1 ] = eqRetD[:-1]
        retD[ eqRetX.shape[0]-1: ] = da[:,4]
        
        return retX, retD
