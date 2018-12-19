#!/usr/bin/env python

import convertvsh5

import sys
import optparse

cpsr = optparse.OptionParser()
cpsr.add_option( '-i', '--input-file', dest='file', action='store',
                 help='Specify the file to annotate' )
cpsr.add_option( '-o', '--output-file', dest='ofile', action='store',
                 help='Specify the file to write to' )

opts, args = cpsr.parse_args(sys.argv)

cvtr = convertvsh5.UeConvertVsH5()
cvtr.readData(opts.file)
cvtr.createUnstructGrid()
cvtr.annotateFile(opts.ofile)
cvtr.annotateStructured(opts.ofile)
