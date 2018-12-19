#!/usr/bin/env python

# @description This script will convert a pdb data file into an hdf5
#              data file
#
# @keywords    Py-UEDGE
#
# @contact     John Cary
#
# @version     $id: $

# import sys for system commands
import sys, os

# Import uefacets for interface convenience
import uefacets

# import getopts for cmd line parsing
import getopt

# import uedge and all its structures
from uedge import *

def usage(code):
    print "pdb2h5 [-h] -i <infile> -o <outfile> [-p <paramfile>]"
    print " -h: Print this help"
    sys.exit(code)

try:
    olst, _ = getopt.getopt(sys.argv[1:], "hi:o:p:")
except:
    usage(1)

for o in olst:
    if o[0] == "-h":
        usage(0)
    if o[0] == "-i":
        infile = o[1]
    if o[0] == "-o":
        outfile = o[1]
    if o[0] == "-p":
        prmfile = o[1]

# initialize the uedge object

uefacets.init()

ue = uefacets.Uedge()

try:
    ue.readParams(prmfile)
except NameError, _:
    print "No parameter file specified, continuing..."

ue.buildData()
ue.restore(infile)
ue.dump(outfile)

uefacets.final()

sys.exit(0)

