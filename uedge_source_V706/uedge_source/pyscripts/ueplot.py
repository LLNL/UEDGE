#!/usr/bin/env python

# Import all required libraries

import optparse
import sys

import numpy
import pylab
import string
import tables

import ueplotdata

#
# This function tests to see if a point lies inside a given polygon
#
# @param polyX array of x point values for the polygon
# @param polyY array of y point values for the polygon
# @param x value of the x coordinate of the point to test
# @param y value of the y coordinate of the point to test
#
# @return bool specifying if point is inside the polygon
#
def polyContains(polyX, polyY, x, y):
    retval = False
    
    j = polyX.shape[0] - 1
    for i in range(0, polyX.shape[0]):
        if ( polyY[i] < y and polyY[j] >= y or
             polyY[j] < y and polyY[i] >= y ):
            if ( polyX[i] + (y-polyY[i]) / (polyY[j]-polyY[i])
                 * (polyX[j]-polyX[i]) < x ):
                retval = not retval
        j = i

    return retval

#
# Construct a 2D mapped grid interpolator
#
# @param xxarry 2d array of x coordinate values of the mapped grid
# @param yyarry 2d array of y coordinate values of the mapped grid
# @param zzvals 2d array of z values at each grid point
# @param nxarry N array of x coordinates at which to interpolate
# @param nyarry N array of y coordinates at which to interpolate
#
# @return N array of values at the interpolation points
#
def mappedInterpolate(xxarry, yyarry, zzvals, nxarry, nyarry):
    xx = nxarry
    yy = nyarry
    zz = numpy.zeros(xx.shape, numpy.float)

    oset = [(0,0), (1, 0), (1, 1), (0, 1)]

    for i in range(0, nxarry.shape[0]):

        doBreak = False
        mval = nval = 0
        for m in range(0, xxarry.shape[0]-1):
            for n in range(0, xxarry.shape[1]-1):

                polyX = numpy.array([ xxarry[m + oset[0][0], n + oset[0][1]],
                                      xxarry[m + oset[1][0], n + oset[1][1]],
                                      xxarry[m + oset[2][0], n + oset[2][1]],
                                      xxarry[m + oset[3][0], n + oset[3][1]] ])
                polyY = numpy.array([ yyarry[m + oset[0][0], n + oset[0][1]],
                                      yyarry[m + oset[1][0], n + oset[1][1]],
                                      yyarry[m + oset[2][0], n + oset[2][1]],
                                      yyarry[m + oset[3][0], n + oset[3][1]] ])

                contVal = polyContains(polyX, polyY, xx[i], yy[i])
                if ( contVal ):
                    print "Box found at ", m, n
                    mval = m
                    nval = n
                    doBreak = True
                    break

            if ( doBreak ):
                break

        if ( xxarry[mval + oset[0][0], nval + oset[0][1]] ==
             xxarry[mval + oset[1][0], nval + oset[1][1]] and
             yyarry[mval + oset[0][0], nval + oset[0][1]] ==
             yyarry[mval + oset[1][0], nval + oset[1][1]] ):

            x1 = xxarry[mval + oset[0][0], nval + oset[0][1]]
            x2 = xxarry[mval + oset[2][0], nval + oset[2][1]]
            x3 = xxarry[mval + oset[3][0], nval + oset[3][1]]

            y1 = yyarry[mval + oset[0][0], nval + oset[0][1]]
            y2 = yyarry[mval + oset[2][0], nval + oset[2][1]]
            y3 = yyarry[mval + oset[3][0], nval + oset[3][1]]

            z1 = zzvals[mval + oset[0][0], nval + oset[0][1]]
            z2 = zzvals[mval + oset[2][0], nval + oset[2][1]]
            z3 = zzvals[mval + oset[3][0], nval + oset[3][1]]

            a = -((-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) \
                / (-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3))
            b = -((-(x2*z1) + x3*z1 + x1*z2 - x3*z2 - x1*z3 + x2*z3) \
                / (x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3))
            c = -((-(x3*y2*z1) + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2 \
                * y1*z3 + x1*y2*z3)/(x2*y1 - x3*y1 - x1*y2 + x3*y2   \
                + x1*y3 - x2*y3))

            fn = lambda x, y: a*x + b*y + c

        elif ( xxarry[mval + oset[1][0], nval + oset[1][1]] ==
               xxarry[mval + oset[2][0], nval + oset[2][1]] and
               yyarry[mval + oset[1][0], nval + oset[1][1]] ==
               yyarry[mval + oset[2][0], nval + oset[2][1]] ):

            x1 = xxarry[mval + oset[0][0], nval + oset[0][1]]
            x2 = xxarry[mval + oset[1][0], nval + oset[1][1]]
            x3 = xxarry[mval + oset[3][0], nval + oset[3][1]]

            y1 = yyarry[mval + oset[0][0], nval + oset[0][1]]
            y2 = yyarry[mval + oset[1][0], nval + oset[1][1]]
            y3 = yyarry[mval + oset[3][0], nval + oset[3][1]]

            z1 = zzvals[mval + oset[0][0], nval + oset[0][1]]
            z2 = zzvals[mval + oset[1][0], nval + oset[1][1]]
            z3 = zzvals[mval + oset[3][0], nval + oset[3][1]]

            a = -((-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) \
                / (-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3))
            b = -((-(x2*z1) + x3*z1 + x1*z2 - x3*z2 - x1*z3 + x2*z3) \
                / (x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3))
            c = -((-(x3*y2*z1) + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2 \
                * y1*z3 + x1*y2*z3)/(x2*y1 - x3*y1 - x1*y2 + x3*y2   \
                + x1*y3 - x2*y3))

            fn = lambda x, y: a*x + b*y + c

        elif ( xxarry[mval + oset[2][0], nval + oset[2][1]] ==
               xxarry[mval + oset[3][0], nval + oset[3][1]] and
               yyarry[mval + oset[2][0], nval + oset[2][1]] ==
               yyarry[mval + oset[3][0], nval + oset[3][1]] ):

            x1 = xxarry[mval + oset[0][0], nval + oset[0][1]]
            x2 = xxarry[mval + oset[1][0], nval + oset[1][1]]
            x3 = xxarry[mval + oset[2][0], nval + oset[2][1]]

            y1 = yyarry[mval + oset[0][0], nval + oset[0][1]]
            y2 = yyarry[mval + oset[1][0], nval + oset[1][1]]
            y3 = yyarry[mval + oset[2][0], nval + oset[2][1]]

            z1 = zzvals[mval + oset[0][0], nval + oset[0][1]]
            z2 = zzvals[mval + oset[1][0], nval + oset[1][1]]
            z3 = zzvals[mval + oset[2][0], nval + oset[2][1]]

            a = -((-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) \
                / (-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3))
            b = -((-(x2*z1) + x3*z1 + x1*z2 - x3*z2 - x1*z3 + x2*z3) \
                / (x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3))
            c = -((-(x3*y2*z1) + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2 \
                * y1*z3 + x1*y2*z3)/(x2*y1 - x3*y1 - x1*y2 + x3*y2   \
                + x1*y3 - x2*y3))

            fn = lambda x, y: a*x + b*y + c

        elif ( xxarry[mval + oset[3][0], nval + oset[3][1]] ==
               xxarry[mval + oset[0][0], nval + oset[0][1]] and
               yyarry[mval + oset[3][0], nval + oset[3][1]] ==
               yyarry[mval + oset[0][0], nval + oset[0][1]] ):

            x1 = xxarry[mval + oset[1][0], nval + oset[1][1]]
            x2 = xxarry[mval + oset[2][0], nval + oset[2][1]]
            x3 = xxarry[mval + oset[3][0], nval + oset[3][1]]

            y1 = yyarry[mval + oset[1][0], nval + oset[1][1]]
            y2 = yyarry[mval + oset[2][0], nval + oset[2][1]]
            y3 = yyarry[mval + oset[3][0], nval + oset[3][1]]

            z1 = zzvals[mval + oset[1][0], nval + oset[1][1]]
            z2 = zzvals[mval + oset[2][0], nval + oset[2][1]]
            z3 = zzvals[mval + oset[3][0], nval + oset[3][1]]

            a = -((-(y2*z1) + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3) \
                / (-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3))
            b = -((-(x2*z1) + x3*z1 + x1*z2 - x3*z2 - x1*z3 + x2*z3) \
                / (x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3))
            c = -((-(x3*y2*z1) + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2 \
                * y1*z3 + x1*y2*z3)/(x2*y1 - x3*y1 - x1*y2 + x3*y2   \
                + x1*y3 - x2*y3))

            fn = lambda x, y: a*x + b*y + c

        else:

            x1 = xxarry[mval + oset[0][0], nval + oset[0][1]]
            x2 = xxarry[mval + oset[1][0], nval + oset[1][1]]
            x3 = xxarry[mval + oset[2][0], nval + oset[2][1]]
            x4 = xxarry[mval + oset[3][0], nval + oset[3][1]]

            y1 = yyarry[mval + oset[0][0], nval + oset[0][1]]
            y2 = yyarry[mval + oset[1][0], nval + oset[1][1]]
            y3 = yyarry[mval + oset[2][0], nval + oset[2][1]]
            y4 = yyarry[mval + oset[3][0], nval + oset[3][1]]

            z1 = zzvals[mval + oset[0][0], nval + oset[0][1]]
            z2 = zzvals[mval + oset[1][0], nval + oset[1][1]]
            z3 = zzvals[mval + oset[2][0], nval + oset[2][1]]
            z4 = zzvals[mval + oset[3][0], nval + oset[3][1]]

            a = (-(x1*x4*(y1 - y4)*(y3*z2 - y2*z3)) + x3*(-(x4*(y3 - y4)*(y2*z1    \
                - y1*z2)) + x1*(y1 - y3)*(y4*z2 - y2*z4)) + x2*(x4*(y2 - y4)       \
                * (y3*z1 - y1*z3) - x3*(y2 - y3)*(y4*z1 - y1*z4) - x1*(y1 - y2)    \
                * (y4*z3 - y3*z4)))/(x2*(x3*(y2 - y3)*(y1 - y4) - x4*(y1 - y3)*(y2 \
                - y4)) + x1*(x4*(y2 - y3)*(y1 - y4) - x3*(y1 - y3)*(y2 - y4)       \
                + x2*(y1 - y2)*(y3 - y4)) + x3*x4*(y1 - y2)*(y3 - y4))
            b = (x4*y4*(y3*(z1 - z2) + y1*(z2 - z3) + y2*(-z1 + z3)) + x2*y2*(y4   \
                * (z1 - z3) + y1*(z3 - z4) + y3*(-z1 + z4)) + x3*y3*(y4*(-z1 + z2) \
                + y2*(z1 - z4) + y1*(-z2 + z4)) + x1*y1*(y4*(-z2 + z3) + y3*(z2    \
                - z4) + y2*(-z3 + z4)))/(x2*(x3*(y2 - y3)*(y1 - y4) - x4*(y1 - y3) \
                * (y2 - y4)) + x1*(x4*(y2 - y3)*(y1 - y4) - x3*(y1 - y3)*(y2 - y4) \
                + x2*(y1 - y2)*(y3 - y4)) + x3*x4*(y1 - y2)*(y3 - y4))
            c = (x1*x4*(y1 - y4)*(z2 - z3) + x3*(x4*(y3 - y4)*(z1 - z2) - x1*(y1   \
                - y3)*(z2 - z4)) + x2*(-(x4*(y2 - y4)*(z1 - z3)) + x3*(y2 - y3)    \
                * (z1 - z4) + x1*(y1 - y2)*(z3 - z4)))/(x2*(x3*(y2 - y3)*(y1 - y4) \
                - x4*(y1 - y3)*(y2 - y4)) + x1*(x4*(y2 - y3)*(y1 - y4) - x3*(y1    \
                - y3)*(y2 - y4) + x2*(y1 - y2)*(y3 - y4)) + x3*x4*(y1 - y2)*(y3    \
                - y4))
            d = (x4*(y3*(-z1 + z2) + y2*(z1 - z3) + y1*(-z2 + z3)) + x3*(y4*(z1    \
                - z2) + y1*(z2 - z4) + y2*(-z1 + z4)) + x1*(y4*(z2 - z3) + y2*(z3  \
                - z4) + y3*(-z2 + z4)) + x2*(y4*(-z1 + z3) + y3*(z1 - z4) + y1     \
                * (-z3 + z4)))/(x2*(x3*(y2 - y3)*(y1 - y4) - x4*(y1 - y3)*(y2      \
                - y4)) + x1*(x4*(y2 - y3)*(y1 - y4) - x3*(y1 - y3)*(y2 - y4) + x2  \
                * (y1 - y2)*(y3 - y4)) + x3*x4*(y1 - y2)*(y3 - y4))

            fn = lambda x, y: a + b*x + c*y + d*x*y

        zz[i] = fn(xx[i], yy[i])
        if ( not doBreak ):
            zz[i] = numpy.nan
            print "Box NOT found!"

    return zz

# get options from the command line

cpsr = optparse.OptionParser()

cpsr.add_option('-e', '--edge-file', dest='efile', action='store',
                help='Specify the edge data file to plot', default='')
cpsr.add_option('-s', '--edge-species', dest='species', action='store',
                help='Specify the edge species temperature to plot (ion, elec)',
                default='ion')

cpsr.add_option('-c', '--core-file', dest='cfile', action='store',
                help='Specify the core data file to plot', default='')

cpsr.add_option('-g', '--grid-file', dest='gfile', action='store',
                help='Specify the grid data file for core plots', default='')

cpsr.add_option('-l', '--log-plot', dest='logbool', action='store_true',
                help='Pass this flag to generate logarithmic scale plots',
                default=False)

cpsr.add_option('-f', '--file-suffix', dest='sufx', action='store',
                help='Use this as a suffix to uniquely name the outfiles',
                default='')

opts, args = cpsr.parse_args(sys.argv)

# Get all the data from the HDF5 file

fh = tables.openFile(opts.efile)

uen = fh.getNode('/bbb/uegrid')
ueg = uen.uegrid.read()

udn = fh.getNode('/bbb')
if opts.species == 'ion':
    ued = udn.tis.read() * 6.21e18
elif opts.species == 'elec':
    ued = udn.tes.read() * 6.21e18
else:
    print 'Invalid species plot requested'
    exit(1)

ifg = fh.getNode('/bbb/uegrid/blockInfo')

ysp = ifg._v_attrs.__getattr__('ySeparatrix')
xsp1 = ifg._v_attrs.__getattr__('xLeftSeparatrix')
xsp2 = ifg._v_attrs.__getattr__('xRightSeparatrix')

yoff = ifg._v_attrs.__getattr__('yMagneticAxis')
xoff = ifg._v_attrs.__getattr__('xMagneticAxis')

fh.close()

udata = ueplotdata.UePlotData(opts.efile)

# Process the private flux region

pfXX, pfYY = udata.getNodalPFGrid()
pfDA = udata.getCenteredPFData()

# Process the core layer

coXX, coYY = udata.getNodalCoreGrid()
coDA = udata.getCenteredCoreData()

# Process the scrape off layer

scXX, scYY = udata.getNodalSOLGrid()
scDA = udata.getCenteredSOLData()

# Read in all the data from the equilibrium data file

fh = open(opts.gfile, 'r')

lns = fh.readlines()[1:]
l1 = lns[0].replace(',', ' ').split()
ii = string.atoi( l1[ l1.index('i=')+1 ] )
jj = string.atoi( l1[ l1.index('j=')+1 ] ) - 5

fh.close()

# Process the core region grid

eqXX = numpy.zeros((ii, jj), numpy.float)
eqYY = numpy.zeros((ii, jj), numpy.float)
eqDA = numpy.zeros((ii, jj), numpy.float)

lns = lns[1:]
k = 0

for j in range(0, jj):
    for i in range(0, ii):
        l = lns[k].split()
        eqXX[i,j] = string.atof( l[0] )
        eqYY[i,j] = string.atof( l[1] )
        eqDA[i,j] = string.atof( l[2] )
        k = k+1

# Read in the data from the core grid fild

fh = tables.openFile(opts.cfile)

con = fh.getNode('/facets/core')
cog = con.qOld.read()
temp = cog[:,1]/(1.5*cog[:,0])*1e3

fh.close()

# Process the core region data

dx = eqDA.max() / (cog.shape[0])

eqDAt = numpy.interp(eqDA.reshape(ii*jj),
                     numpy.arange(dx / 2, eqDA.max(), dx),
                     temp)

eqDA = eqDAt.reshape(eqDA.shape)

eqDA = eqDA

# Construct the slice arrays for midplane plots

pfiXX = numpy.zeros((ueg.shape[0] - xsp2 + xsp1 - 2, ysp), numpy.float)
pfiYY = numpy.zeros((ueg.shape[0] - xsp2 + xsp1 - 2, ysp), numpy.float)

pfmXX = numpy.zeros((ueg.shape[0] - xsp2 + xsp1 - 2, ysp, 5), numpy.float)
pfmYY = numpy.zeros((ueg.shape[0] - xsp2 + xsp1 - 2, ysp, 5), numpy.float)

pfiXX[:xsp1, :] = ueg[1:xsp1+1,  1:ysp+1, 0, 0]
pfiXX[xsp1:, :] = ueg[xsp2+1:-1, 1:ysp+1, 0, 0]

pfmXX[:xsp1, :] = ueg[1:xsp1+1,  1:ysp+1, :, 0]
pfmXX[xsp1:, :] = ueg[xsp2+1:-1, 1:ysp+1, :, 0]

pfiYY[:xsp1, :] = ueg[1:xsp1+1,  1:ysp+1, 0, 1]
pfiYY[xsp1:, :] = ueg[xsp2+1:-1, 1:ysp+1, 0, 1]

pfmYY[:xsp1, :] = ueg[1:xsp1+1,  1:ysp+1, :, 1]
pfmYY[xsp1:, :] = ueg[xsp2+1:-1, 1:ysp+1, :, 1]

coiXX = numpy.zeros((xsp2 - xsp1 + 1, ysp + 1), numpy.float)
coiYY = numpy.zeros((xsp2 - xsp1 + 1, ysp + 1), numpy.float)

cocXX = numpy.zeros((xsp2 - xsp1 + 1, ysp + 1), numpy.float)
cocYY = numpy.zeros((xsp2 - xsp1 + 1, ysp + 1), numpy.float)

cocXX = numpy.zeros((xsp2 - xsp1 + 1, ysp + 1), numpy.float)
cocYY = numpy.zeros((xsp2 - xsp1 + 1, ysp + 1), numpy.float)

comXX = numpy.zeros((xsp2 - xsp1 + 1, ysp + 1, 5), numpy.float)
comYY = numpy.zeros((xsp2 - xsp1 + 1, ysp + 1, 5), numpy.float)

coiDA = numpy.zeros((xsp2 - xsp1 + 1, ysp + 1), numpy.float)

cocXX[:-1, :] = ueg[xsp1+1:xsp2+1, 0:ysp+1, 0, 0]
cocXX[-1, :] = ueg[xsp1+1, 0:ysp+1, 0, 0]

coiXX[:-1, :] = ueg[xsp1+1:xsp2+1, 0:ysp+1, 1, 0]
coiXX[-1, :] = ueg[xsp1+1, 0:ysp+1, 1, 0]

comXX[:-1, :] = ueg[xsp1+1:xsp2+1, 0:ysp+1, :, 0]
comXX[-1, :] = ueg[xsp1+1, 0:ysp+1, :, 0]

cocYY[:-1, :] = ueg[xsp1+1:xsp2+1, 0:ysp+1, 0, 1] - yoff
cocYY[-1, :] = ueg[xsp1+1, 0:ysp+1, 0, 1] - yoff

coiYY[:-1, :] = ueg[xsp1+1:xsp2+1, 0:ysp+1, 1, 1] - yoff
coiYY[-1, :] = ueg[xsp1+1, 0:ysp+1, 1, 1] - yoff

comYY[:-1, :] = ueg[xsp1+1:xsp2+1, 0:ysp+1, :, 1] - yoff
comYY[-1, :] = ueg[xsp1+1, 0:ysp+1, :, 1] - yoff

coiDA[:-1, :] = ued[xsp1+1:xsp2+1, 1:ysp+2]
coiDA[-1, :] = ued[xsp1+1, 1:ysp+2]

sciXX = numpy.zeros((ueg.shape[0] - 2, ueg.shape[1] - ysp - 2), numpy.float)
sciYY = numpy.zeros((ueg.shape[0] - 2, ueg.shape[1] - ysp - 2), numpy.float)

scmXX = numpy.zeros((ueg.shape[0] - 2, ueg.shape[1] - ysp - 2, 5), numpy.float)
scmYY = numpy.zeros((ueg.shape[0] - 2, ueg.shape[1] - ysp - 2, 5), numpy.float)

sciXX[:, :] = ueg[1:-1, ysp+1:-1, 0, 0]
sciYY[:, :] = ueg[1:-1, ysp+1:-1, 0, 1] - yoff

scmXX[:, :] = ueg[1:-1, ysp+1:-1, :, 0]
scmYY[:, :] = ueg[1:-1, ysp+1:-1, :, 1] - yoff

fh = open('uedge-coregrid.dat', 'w')

fh.write( 'i = %d, j = %d\n' % (comXX.shape[0], comYY.shape[1]) )

for j in range(0, comYY.shape[1]):
    for i in range(0, comXX.shape[0]):
        fh.write( '%e %e %e\n' % (comXX[i, j, 1], comYY[i, j, 1], coiDA[i, j]) )

fh.close()

fh = open('uedge-solgrid.dat', 'w')

fh.write( 'i = %d, j = %d\n' % (scmXX.shape[0], scmYY.shape[1]) )

for j in range(0, scmYY.shape[1]):
    for i in range(0, scmXX.shape[0]):
        fh.write( '%e %e %e\n' % (scmXX[i, j, 1], scmYY[i, j, 1], scDA[i, j]) )

fh.close()

fh = open('uedge-pflgrid.dat', 'w')

fh.write( 'i = %d, j = %d\n' % (pfmXX.shape[0], pfmYY.shape[1]) )

for j in range(0, pfmYY.shape[1]):
    for i in range(0, pfmXX.shape[0]):
        fh.write( '%e %e %e\n' % (pfmXX[i, j, 1], pfmYY[i, j, 1], pfDA[i, j]) )

fh.close()

# Plot the data

# xsp = ueg[xsp1+1, ysp+2, 1, 0]
# ysp = ueg[xsp1+1, ysp+2, 1, 1]

# xxs1 = numpy.arange(xoff, xsp + 0.0001, (xsp - xoff) / 150)
# yys1 = numpy.arange(0, ysp - yoff, (ysp - yoff) / 150)

# sln1 = mappedInterpolate(eqXX, eqYY, eqDA, xxs1, yys1)
# sln2 = mappedInterpolate(coiXX, coiYY, coDA, xxs1, yys1)

# pylab.figure(0)
# pylab.plot(numpy.sqrt((xxs1 - xoff)**2 + (yys1)**2), sln1)
# pylab.plot(numpy.sqrt((xxs1 - xoff)**2 + (yys1)**2), sln2)
# pylab.show()

fig1 = pylab.figure(1)

if ( opts.logbool ):
    ziVal = numpy.log10(numpy.min([pfDA.min(), coDA.min(), scDA.min()]))
    zaVal = numpy.log10(numpy.max([pfDA.max(), coDA.max(), scDA.max()]))
else:
    ziVal = (numpy.min([pfDA.min(), coDA.min(), scDA.min()]))
    zaVal = (numpy.max([pfDA.max(), coDA.max(), scDA.max()]))

if ( opts.logbool ):
    pylab.pcolor(pfXX, pfYY-yoff, numpy.log10(pfDA), shading='flat')
    pylab.clim(ziVal, 2*zaVal-ziVal)

    pylab.pcolor(coXX, coYY-yoff, numpy.log10(coDA), shading='flat')
    pylab.clim(ziVal, 2*zaVal-ziVal)

    pylab.pcolor(scXX, scYY-yoff, numpy.log10(scDA), shading='flat')
    pylab.clim(ziVal, 2*zaVal-ziVal)
else:
    pylab.pcolor(pfXX, pfYY-yoff, (pfDA), shading='flat')
    pylab.clim(ziVal, 2*zaVal-ziVal)

    pylab.pcolor(coXX, coYY-yoff, (coDA), shading='flat')
    pylab.clim(ziVal, 2*zaVal-ziVal)

    pylab.pcolor(scXX, scYY-yoff, (scDA), shading='flat')
    pylab.clim(ziVal, 2*zaVal-ziVal)

if ( opts.logbool ):
    iVal = numpy.log10(eqDA.min())
    aVal = numpy.log10(eqDA.max())
else:
    iVal = (eqDA.min())
    aVal = (eqDA.max())

if ( opts.logbool ):
    pylab.pcolor(eqXX, eqYY, numpy.log10(eqDA), shading='flat')
    pylab.clim(2*iVal-aVal, aVal)
else:
    pylab.pcolor(eqXX, eqYY, (eqDA), shading='flat')
    pylab.clim(2*iVal-aVal, aVal)
    
pylab.xlabel('R')
pylab.ylabel('Z')

xiVal = scXX.min()
xaVal = scXX.max()
xrVal = xaVal - xiVal

yiVal = pfYY.min() - yoff
yaVal = scYY.max() - yoff
yrVal = yaVal - yiVal

fig1.gca().set_xlim([xiVal - 0.1*xrVal, xaVal + 0.1*xrVal])
fig1.gca().set_ylim([yiVal - 0.1*yrVal, yaVal + 0.1*yrVal])
fig1.gca().set_aspect('equal')

# if ( opts.logbool ):
#     pylab.colorbar(format='1.0 x %.0f')
# else:
#     pylab.colorbar()

fig1.savefig('coupled-overplot%s.png' % opts.sufx)

fig2 = pylab.figure(2)

iVal = (ued.min())
aVal = (ued.max())

pylab.pcolor(pfXX, pfYY-yoff, (pfDA), shading='flat')
pylab.clim(iVal, aVal)

pylab.pcolor(coXX, coYY-yoff, (coDA), shading='flat')
pylab.clim(iVal, aVal)

pylab.pcolor(scXX, scYY-yoff, (scDA), shading='flat')
pylab.clim(iVal, aVal)

pylab.xlabel('R')
pylab.ylabel('Z')
pylab.title('Maximum: %g' % aVal)

fig2.gca().set_xlim([xiVal - 0.1*xrVal, xaVal + 0.1*xrVal])
fig2.gca().set_ylim([yiVal - 0.1*yrVal, yaVal + 0.1*yrVal])
fig2.gca().set_aspect('equal')

pylab.colorbar()

fig2.savefig('edge-section%s.png' % opts.sufx)

fig3 = pylab.figure(3)

iVal = (eqDA.min())
aVal = (eqDA.max())

pylab.pcolor(eqXX, eqYY, (eqDA), shading='flat')
pylab.clim(iVal, aVal)

pylab.xlabel('R')
pylab.ylabel('Z')
pylab.title('Minimum: %g' % iVal)

fig3.gca().set_xlim([xiVal - 0.1*xrVal, xaVal + 0.1*xrVal])
fig3.gca().set_ylim([yiVal - 0.1*yrVal, yaVal + 0.1*yrVal])
fig3.gca().set_aspect('equal')

pylab.colorbar()

fig3.savefig('core-section%s.png' % opts.sufx)

print "********** Core Interpolation **********"
xx1 = numpy.arange(xoff, eqXX.max()+0.0001, (eqXX.max() - xoff) / 100)
yy1 = numpy.zeros(xx1.shape, numpy.float)
sln1 = mappedInterpolate(eqXX, eqYY, eqDA, xx1, yy1)

print "********** Edge Interpolation **********"
xx2 = numpy.arange(eqXX.max()-0.02, coXX.max() + 0.0001, (coXX.max() - eqXX.max()) / 25)
yy2 = numpy.zeros(xx2.shape, numpy.float)
sln2 = mappedInterpolate(coXX, coYY, udata.getNodalCoreData(), xx2, yy2)

print "********** SOL  Interpolation **********"
xx3 = numpy.arange(coiXX.max(), sciXX.max() + 0.0001, (sciXX.max() - coiXX.max()) / 25)
yy3 = numpy.zeros(xx3.shape, numpy.float)
sln3 = mappedInterpolate(scXX, scYY, udata.getNodalSOLData(), xx3, yy3)

fig4 = pylab.figure(4)

pylab.plot(xx1, sln1)
pylab.plot(xx2, sln2)
pylab.plot(xx3, sln3)

fig4.gca().set_xlim([2.2, 2.3])
fig4.gca().set_ylim([0, 1000 ])

pylab.xlabel('Arbitrary')
pylab.ylabel('Temperature')

fig4.savefig('midplane-slice%s.png' % opts.sufx)
#pylab.show()
