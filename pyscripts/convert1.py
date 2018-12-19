#!/usr/bin/env python
#
# $Id: convert1.py,v 7.0 2018/02/28 18:43:48 meyer8 Exp $
#
# To try solving linear critical gradient

import sys
import os
import getopt
import string
import convert
from convert import *

# define the mppl to f90 class
class M2F(generic):
  suffixin = "m"
  suffixout = "F"
  subrules = globalsubrules + M2Fsubrules

def usage():
  print "Usage: convert1.py -i <indir> -o <outdir> <infile>"

r"""
main(argv: array of strings)
"""
def main(argv):
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:", ["help",
	"indir=", "outdir="])
  except getopt.GetoptError:
# print help information and exit:
    usage()
    sys.exit(2)

  indir = "."
  outdir = "."
# Go through args
  for o, a in opts:
    if o in ("-h", "--help"):
      usage()
      sys.exit()
    if o in ("-i", "--indir"):
      indir = a
    elif o in ("-o", "--outdir"):
      outdir = a
  # print "args =", args

# Extract the file
  fn = args[0]
  # print "Converting " + fn

# Do the conversion
  convert.M2F = M2F
  m2f = M2F(indir, outdir)
  # m2f.outdir = outdir
  m2f.processfile(fn)

  return 0

if __name__ == '__main__':
  sys.exit(main(sys.argv))

