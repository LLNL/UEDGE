# Import all required libraries

import numpy
import tables

#
# The purpose of this class is to read in a UEDGE HDF5 file and then
# make it available in a friendly fashion for plotting
#
class UePlotData:

#
# Constructor
#
# @param fname the file from which this PlotData object should read
#
    def __init__(self, fname):
        self.myName = fname
        self.readData(fname)

#
# This function is called to read data from a given HDF5 file
#
# @param fname the name of the file from which to read
#
    def readData(self, fname, block='tis'):
        fh = tables.openFile(fname)

        # Read the Uedge grid information
        uen = fh.getNode('/bbb/uegrid')
        self.ueg = uen.uegrid.read()
        self.psi = uen.psiVals.read()

        # Read the ion temperature data
        udn = fh.getNode('/bbb')
        self.ued = udn.__getattr__(block).read()

        # Read the additional metadata associated with the branch
        # cut grid
        ifg = fh.getNode('/bbb/uegrid/blockInfo')

        self.ysp  = ifg._v_attrs.__getattr__('ySeparatrix')
        self.xsp1 = ifg._v_attrs.__getattr__('xLeftSeparatrix')
        self.xsp2 = ifg._v_attrs.__getattr__('xRightSeparatrix')

        self.yoff = ifg._v_attrs.__getattr__('yMagneticAxis')
        self.xoff = ifg._v_attrs.__getattr__('xMagneticAxis')

        fh.close()

#
# Returns a tuple of numpy arrays representing the mapped grid of the
# private flux region at the cell nodes.  Grid is +1 on each side since
# all cells are completed
#
# @return a tuple of numpy arrays representing the mapped grid
#     1: The XX array of the mapped grid
#     2: The YY array of the mapped grid
#
    def getNodalPFGrid(self):
        xx = numpy.zeros( (self.ueg.shape[0] - self.xsp2 + self.xsp1 - 1,
                           self.ysp + 1), numpy.float )
        yy = numpy.zeros( (self.ueg.shape[0] - self.xsp2 + self.xsp1 - 1,
                           self.ysp + 1), numpy.float )

        xx[:self.xsp1, :] = self.ueg[1:self.xsp1+1, 1:self.ysp+2, 1, 0]
        xx[self.xsp1:, :] = self.ueg[self.xsp2+1:,  1:self.ysp+2, 1, 0]

        yy[:self.xsp1, :] = self.ueg[1:self.xsp1+1, 1:self.ysp+2, 1, 1] - self.yoff
        yy[self.xsp1:, :] = self.ueg[self.xsp2+1:,  1:self.ysp+2, 1, 1] - self.yoff

        return xx, yy

#
# Returns a tuple of numpy arrays representing the mapped grid of the
# scrape-off region at the cell nodes.  Grid is +1 on each side since
# all cells are completed
#
# @return a tuple of numpy arrays representing the mapped grid
#     1: The XX array of the mapped grid
#     2: The YY array of the mapped grid
#
    def getNodalSOLGrid(self):
        xx = numpy.zeros( (self.ueg.shape[0] - 1,
                           self.ueg.shape[1] - self.ysp - 1), numpy.float )
        yy = numpy.zeros( (self.ueg.shape[0] - 1,
                           self.ueg.shape[1] - self.ysp - 1), numpy.float )

        xx[:, :] = self.ueg[1:, self.ysp+1:, 1, 0]

        yy[:, :] = self.ueg[1:, self.ysp+1:, 1, 1] - self.yoff

        return xx, yy

#
# Returns a tuple of numpy arrays representing the mapped grid of the
# core region at the cell nodes.  Grid is +1 on each side since all
# cells are completed
#
# @return a tuple of numpy arrays representing the mapped grid
#     1: The XX array of the mapped grid
#     2: The YY array of the mapped grid
#
    def getNodalCoreGrid(self):
        xx = numpy.zeros( (self.xsp2 - self.xsp1 + 1,
                           self.ysp + 1), numpy.float )
        yy = numpy.zeros( (self.xsp2 - self.xsp1 + 1,
                           self.ysp + 1), numpy.float )

        xx[:-1, :] = self.ueg[self.xsp1+1:self.xsp2+1, 1:self.ysp+2, 1, 0]
        xx[-1:, :] = self.ueg[self.xsp1+1,             1:self.ysp+2, 1, 0]

        yy[:-1, :] = self.ueg[self.xsp1+1:self.xsp2+1, 1:self.ysp+2, 1, 1] - self.yoff
        yy[-1:, :] = self.ueg[self.xsp1+1,             1:self.ysp+2, 1, 1] - self.yoff

        return xx, yy

#
# Returns a tuple of numpy arrays representing the mapped grid of the
# private flux region at the cell centers.
#
# @return a tuple of numpy arrays representing the mapped grid
#     1: The XX array of the mapped grid
#     2: The YY array of the mapped grid
#
    def getCenteredPFGrid(self):
        xx = numpy.zeros( (self.ueg.shape[0] - self.xsp2 + self.xsp1 - 2,
                           self.ysp), numpy.float )
        yy = numpy.zeros( (self.ueg.shape[0] - self.xsp2 + self.xsp1 - 2,
                           self.ysp), numpy.float )

        xx[:self.xsp1, :] = self.ueg[1:self.xsp1+1,  1:self.ysp+1, 0, 0]
        xx[self.xsp1:, :] = self.ueg[self.xsp2+1:-1, 1:self.ysp+1, 0, 0]

        yy[:self.xsp1, :] = self.ueg[1:self.xsp1+1,  1:self.ysp+1, 0, 1] - self.yoff
        yy[self.xsp1:, :] = self.ueg[self.xsp2+1:-1, 1:self.ysp+1, 0, 1] - self.yoff

        return xx, yy

#
# Returns a tuple of numpy arrays representing the mapped grid of the
# scrape-off region at the cell centers.
#
# @return a tuple of numpy arrays representing the mapped grid
#     1: The XX array of the mapped grid
#     2: The YY array of the mapped grid
#
    def getCenteredSOLGrid(self):
        xx = numpy.zeros( (self.ueg.shape[0] - 2,
                           self.ueg.shape[1] - self.ysp - 2), numpy.float )
        yy = numpy.zeros( (self.ueg.shape[0] - 2,
                           self.ueg.shape[1] - self.ysp - 2), numpy.float )

        xx[:, :] = self.ueg[1:-1, self.ysp+1:-1, 0, 0]

        yy[:, :] = self.ueg[1:-1, self.ysp+1:-1, 0, 1] - self.yoff

        return xx, yy

#
# Returns a tuple of numpy arrays representing the mapped grid of the
# core region at the cell centers.
#
# @return a tuple of numpy arrays representing the mapped grid
#     1: The XX array of the mapped grid
#     2: The YY array of the mapped grid
#
    def getCenteredCoreGrid(self):
        xx = numpy.zeros( (self.xsp2 - self.xsp1 + 1,
                           self.ysp), numpy.float )
        yy = numpy.zeros( (self.xsp2 - self.xsp1 + 1,
                           self.ysp), numpy.float )

        xx[:-1, :] = self.ueg[self.xsp1+1:self.xsp2+1, 1:self.ysp+1, 0, 0]
        xx[-1:, :] = self.ueg[self.xsp1+1,             1:self.ysp+1, 0, 0]

        yy[:-1, :] = self.ueg[self.xsp1+1:self.xsp2+1, 1:self.ysp+1, 0, 1] - self.yoff
        yy[-1:, :] = self.ueg[self.xsp1+1,             1:self.ysp+1, 0, 1] - self.yoff

        return xx, yy

#
# Returns a numpy array representing the uedge psi values at the centers
# of every cell in the prival flux region.
#
# @return a numpy array of the psi values
#
    def getCenteredPFPsi(self):
        ps = numpy.zeros( (self.ueg.shape[0] - self.xsp2 + self.xsp1 - 2,
                           self.ysp), numpy.float )

        ps[:self.xsp1, :] = self.psi[1:self.xsp1+1,  1:self.ysp+1, 0]
        ps[self.xsp1:, :] = self.psi[self.xsp2+1:-1, 1:self.ysp+1, 0]

        return ps

#
# Returns a numpy array representing the uedge psi values at the centers
# of every cell in the scrape off layer.
#
# @return a numpy array of the psi values
#
    def getCenteredSOLPsi(self):
        ps = numpy.zeros( (self.ueg.shape[0] - 2,
                           self.ueg.shape[1] - self.ysp - 2), numpy.float )

        ps[:, :] = self.psi[1:-1, self.ysp+1:-1, 0]

        return ps

#
# Returns a numpy array representing the uedge psi values at the centers
# of every cell in the core region.
#
# @return a numpy array of the psi values
#
    def getCenteredCorePsi(self):
        ps = numpy.zeros( (self.xsp2 - self.xsp1 + 1,
                           self.ysp), numpy.float )

        ps[:-1, :] = self.psi[self.xsp1+1:self.xsp2+1, 1:self.ysp+1, 0]
        ps[-1:, :] = self.psi[self.xsp1+1,             1:self.ysp+1, 0]

        return ps

#
# Returns a numpy array represeing the uedge data at the cell nodes
# of each grid cell in the private flux region
#
# @return a numpy array of the output data
#
    def getNodalPFData(self):
        da = numpy.zeros( (self.ueg.shape[0] - self.xsp2 + self.xsp1 - 1,
                           self.ysp + 1), numpy.float )

        for i in range(1, self.xsp1+1):
            for j in range(1, self.ysp+2):
                xx = numpy.array( [self.ueg[i-1, j-1, 0, 0],
                                   self.ueg[i  , j-1, 0, 0],
                                   self.ueg[i  , j  , 0, 0],
                                   self.ueg[i-1, j  , 0, 0]] )
                yy = numpy.array( [self.ueg[i-1, j-1, 0, 1],
                                   self.ueg[i  , j-1, 0, 1],
                                   self.ueg[i  , j  , 0, 1],
                                   self.ueg[i-1, j  , 0, 1]] )
                zz = numpy.array( [self.ued[i-1, j-1],
                                   self.ued[i  , j-1],
                                   self.ued[i  , j  ],
                                   self.ued[i-1, j  ]] )

                fn = self.bilinearInterpolate(xx, yy, zz)

                da[i-1, j-1] = fn( self.ueg[i-1, j-1, 1, 0],
                                   self.ueg[i-1, j-1, 1, 1] )

        for i in range(self.xsp2+1, ueg.shape[0]-1):
            for j in range(1, self.ysp+2):
                xx = numpy.array( [self.ueg[i-1, j-1, 0, 0],
                                   self.ueg[i  , j-1, 0, 0],
                                   self.ueg[i  , j  , 0, 0],
                                   self.ueg[i-1, j  , 0, 0]] )
                yy = numpy.array( [self.ueg[i-1, j-1, 0, 1],
                                   self.ueg[i  , j-1, 0, 1],
                                   self.ueg[i  , j  , 0, 1],
                                   self.ueg[i-1, j  , 0, 1]] )
                zz = numpy.array( [self.ued[i-1, j-1],
                                   self.ued[i  , j-1],
                                   self.ued[i  , j  ],
                                   self.ued[i-1, j  ]] )

                fn = self.bilinearInterpolate(xx, yy, zz)

                da[i-1, j-1] = fn( self.ueg[i-1, j-1, 1, 0],
                                   self.ueg[i-1, j-1, 1, 1] )

        return da

#
# Returns a numpy array represeing the uedge data at the cell nodes
# of each grid cell in the scrape-off region
#
# @return a numpy array of the output data
#
    def getNodalSOLData(self):
        da = numpy.zeros( (self.ueg.shape[0] - 2,
                           self.ueg.shape[1] - self.ysp - 2), numpy.float)

        for i in range(1, self.ueg.shape[0]-1):
            for j in range(self.ysp+1, self.ueg.shape[1]-1):
                xx = numpy.array( [self.ueg[i-1, j-1, 0, 0],
                                   self.ueg[i  , j-1, 0, 0],
                                   self.ueg[i  , j  , 0, 0],
                                   self.ueg[i-1, j  , 0, 0]] )
                yy = numpy.array( [self.ueg[i-1, j-1, 0, 1],
                                   self.ueg[i  , j-1, 0, 1],
                                   self.ueg[i  , j  , 0, 1],
                                   self.ueg[i-1, j  , 0, 1]] )
                zz = numpy.array( [self.ued[i-1, j-1],
                                   self.ued[i  , j-1],
                                   self.ued[i  , j  ],
                                   self.ued[i-1, j  ]] )

                fn = self.bilinearInterpolate(xx, yy, zz)

                da[i-1, j-self.ysp-1] = fn( self.ueg[i-1, j-1, 1, 0],
                                            self.ueg[j-1, j-1, 1, 1] )

        return da

#
# Returns a numpy array represeing the uedge data at the cell nodes
# of each grid cell in the core region
#
# @return a numpy array of the output data
#
    def getNodalCoreData(self):
        da = numpy.zeros( (self.xsp2 - self.xsp1 + 1,
                           self.ysp + 1), numpy.float )

        for i in range(self.xsp1+1, self.xsp2+1):
            for j in range(1, self.ysp+2):
                xx = numpy.array( [self.ueg[i-1, j-1, 0, 0],
                                   self.ueg[i  , j-1, 0, 0],
                                   self.ueg[i  , j  , 0, 0],
                                   self.ueg[i-1, j  , 0, 0]] )
                yy = numpy.array( [self.ueg[i-1, j-1, 0, 1],
                                   self.ueg[i  , j-1, 0, 1],
                                   self.ueg[i  , j  , 0, 1],
                                   self.ueg[i-1, j  , 0, 1]] )
                zz = numpy.array( [self.ued[i-1, j-1],
                                   self.ued[i  , j-1],
                                   self.ued[i  , j  ],
                                   self.ued[i-1, j  ]] )

                fn = self.bilinearInterpolate(xx, yy, zz)

                da[i-self.xsp1-1, j-1] = fn( self.ueg[i-1, j-1, 1, 0],
                                             self.ueg[i-1, j-1, 1, 1] )

        return da

#
# Returns a numpy array representing the uedge data at the cell centers
# of each grid cell in the private flux region.
#
# @return a numpy array of the output data
#
    def getCenteredPFData(self, block='tis'):
        self.readData(self.myName, block)

        if ( self.ued.shape.__len__() == 3 ):
            cmpts = self.ued.shape[2]
            deep = True
        else:
            cmpts = 1
            deep = False
        da = numpy.zeros( (self.ueg.shape[0] - self.xsp2 + self.xsp1 - 2,
                           self.ysp, cmpts), numpy.float )

        if ( deep ):
            da[:self.xsp1, :] = self.ued[1:self.xsp1+1,  1:self.ysp+1]
            da[self.xsp1:, :] = self.ued[self.xsp2+1:-1, 1:self.ysp+1]
        else:
            da[:self.xsp1, :, 0] = self.ued[1:self.xsp1+1,  1:self.ysp+1]
            da[self.xsp1:, :, 0] = self.ued[self.xsp2+1:-1, 1:self.ysp+1]

        return da

#
# Returns a numpy array representing the uedge data at the cell centers
# of each grid cell in the scrape-off region
#
# @return a numpy array of the output data
#
    def getCenteredSOLData(self, block='tis'):
        self.readData(self.myName, block)

        if ( self.ued.shape.__len__() == 3 ):
            cmpts = self.ued.shape[2]
            deep = True
        else:
            cmpts = 1
            deep = False
        da = numpy.zeros( (self.ueg.shape[0] - 2,
                           self.ueg.shape[1] - self.ysp - 2, cmpts), numpy.float )

        if ( deep ):
            da[:, :] = self.ued[1:-1, self.ysp+2:]
        else:
            da[:, :, 0] = self.ued[1:-1, self.ysp+2:]

        return da

#
# Returns a numpy array representing the uedge data at the cell centers
# of each grid cell in the core region
#
# @return a numpy array of the output data
#
    def getCenteredCoreData(self, block='tis'):
        self.readData(self.myName, block)

        if ( self.ued.shape.__len__() == 3 ):
            cmpts = self.ued.shape[2]
            deep = True
        else:
            cmpts = 1
            deep = False
        da = numpy.zeros( (self.xsp2 - self.xsp1 + 1,
                           self.ysp, cmpts), numpy.float )

        if ( deep ):
            da[:-1, :] = self.ued[self.xsp1+1:self.xsp2+1, 1:self.ysp+1]
            da[-1:, :] = self.ued[self.xsp1+1,             1:self.ysp+1]
        else:
            da[:-1, :, 0] = self.ued[self.xsp1+1:self.xsp2+1, 1:self.ysp+1]
            da[-1:, :, 0] = self.ued[self.xsp1+1,             1:self.ysp+1]

        return da

#
# Returns a lambda function that can be used to query the value of a
# bilinear interpolation given 4 coordinates and the value at each
# coordinate
#
# @param xx x coordinates of the points to be interpolated
# @param yy y coordinates of the points to be interpolated
# @param zz z values of the points to be interpolated
#
# @return a lambda function representing the interpolation
#
    def bilinearInterpolate(self, xx, yy, zz):

        # These coefficients were auto-generated by Mathematica
        a = (-(xx[0]*xx[3]*(yy[0] - yy[3])*(yy[2]*zz[1] - yy[1]*zz[2])) + xx[2]*(-(xx[3]*(yy[2] - yy[3])*(yy[1]*zz[0]    \
            - yy[0]*zz[1])) + xx[0]*(yy[0] - yy[2])*(yy[3]*zz[1] - yy[1]*zz[3])) + xx[1]*(xx[3]*(yy[1] - yy[3])          \
            * (yy[2]*zz[0] - yy[0]*zz[2]) - xx[2]*(yy[1] - yy[2])*(yy[3]*zz[0] - yy[0]*zz[3]) - xx[0]*(yy[0] - yy[1])    \
            * (yy[3]*zz[2] - yy[2]*zz[3])))/(xx[1]*(xx[2]*(yy[1] - yy[2])*(yy[0] - yy[3]) - xx[3]*(yy[0] - yy[2])*(yy[1] \
            - yy[3])) + xx[0]*(xx[3]*(yy[1] - yy[2])*(yy[0] - yy[3]) - xx[2]*(yy[0] - yy[2])*(yy[1] - yy[3])             \
            + xx[1]*(yy[0] - yy[1])*(yy[2] - yy[3])) + xx[2]*xx[3]*(yy[0] - yy[1])*(yy[2] - yy[3]))
        b = (xx[3]*yy[3]*(yy[2]*(zz[0] - zz[1]) + yy[0]*(zz[1] - zz[2]) + yy[1]*(-zz[0] + zz[2])) + xx[1]*yy[1]*(yy[3]   \
            * (zz[0] - zz[2]) + yy[0]*(zz[2] - zz[3]) + yy[2]*(-zz[0] + zz[3])) + xx[2]*yy[2]*(yy[3]*(-zz[0] + zz[1])    \
            + yy[1]*(zz[0] - zz[3]) + yy[0]*(-zz[1] + zz[3])) + xx[0]*yy[0]*(yy[3]*(-zz[1] + zz[2]) + yy[2]*(zz[1]       \
            - zz[3]) + yy[1]*(-zz[2] + zz[3])))/(xx[1]*(xx[2]*(yy[1] - yy[2])*(yy[0] - yy[3]) - xx[3]*(yy[0] - yy[2])    \
            * (yy[1] - yy[3])) + xx[0]*(xx[3]*(yy[1] - yy[2])*(yy[0] - yy[3]) - xx[2]*(yy[0] - yy[2])*(yy[1] - yy[3])    \
            + xx[1]*(yy[0] - yy[1])*(yy[2] - yy[3])) + xx[2]*xx[3]*(yy[0] - yy[1])*(yy[2] - yy[3]))
        c = (xx[0]*xx[3]*(yy[0] - yy[3])*(zz[1] - zz[2]) + xx[2]*(xx[3]*(yy[2] - yy[3])*(zz[0] - zz[1]) - xx[0]*(yy[0]   \
            - yy[2])*(zz[1] - zz[3])) + xx[1]*(-(xx[3]*(yy[1] - yy[3])*(zz[0] - zz[2])) + xx[2]*(yy[1] - yy[2])          \
            * (zz[0] - zz[3]) + xx[0]*(yy[0] - yy[1])*(zz[2] - zz[3])))/(xx[1]*(xx[2]*(yy[1] - yy[2])*(yy[0] - yy[3])    \
            - xx[3]*(yy[0] - yy[2])*(yy[1] - yy[3])) + xx[0]*(xx[3]*(yy[1] - yy[2])*(yy[0] - yy[3]) - xx[2]*(yy[0]       \
            - yy[2])*(yy[1] - yy[3]) + xx[1]*(yy[0] - yy[1])*(yy[2] - yy[3])) + xx[2]*xx[3]*(yy[0] - yy[1])*(yy[2]       \
            - yy[3]))
        d = (xx[3]*(yy[2]*(-zz[0] + zz[1]) + yy[1]*(zz[0] - zz[2]) + yy[0]*(-zz[1] + zz[2])) + xx[2]*(yy[3]*(zz[0]       \
            - zz[1]) + yy[0]*(zz[1] - zz[3]) + yy[1]*(-zz[0] + zz[3])) + xx[0]*(yy[3]*(zz[1] - zz[2]) + yy[1]*(zz[2]     \
            - zz[3]) + yy[2]*(-zz[1] + zz[3])) + xx[1]*(yy[3]*(-zz[0] + zz[2]) + yy[2]*(zz[0] - zz[3]) + yy[0]           \
            * (-zz[2] + zz[3])))/(xx[1]*(xx[2]*(yy[1] - yy[2])*(yy[0] - yy[3]) - xx[3]*(yy[0] - yy[2])*(yy[1]            \
            - yy[3])) + xx[0]*(xx[3]*(yy[1] - yy[2])*(yy[0] - yy[3]) - xx[2]*(yy[0] - yy[2])*(yy[1] - yy[3]) + xx[1]     \
            * (yy[0] - yy[1])*(yy[2] - yy[3])) + xx[2]*xx[3]*(yy[0] - yy[1])*(yy[2] - yy[3]))

        fn = lambda x, y: a + b*x + c*y + d*x*y

        return fn
