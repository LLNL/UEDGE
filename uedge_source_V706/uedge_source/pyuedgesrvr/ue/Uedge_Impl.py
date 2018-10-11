#
# File:          Uedge_Impl.py
# Symbol:        ue.Uedge-v0.2
# Symbol Type:   class
# Babel Version: 1.5.0 (Revision: 6860 trunk)
# Description:   Implementation of sidl class ue.Uedge in Python.
# 
# WARNING: Automatically generated; changes will be lost
# 
#


import uefacets
import ue.Uedge_Aux
import facets.Facets0dIfc
import facets.FacetsIfc
import sidl.BaseClass
import sidl.BaseInterface
import sidl.ClassInfo
import sidl.RuntimeException
import ue.Uedge
import sidl.NotImplementedException

class Uedge(ue.Uedge_Aux.Uedge):

  def __init__(self, IORself = None):
    ue.Uedge_Aux.Uedge.__init__(self, IORself)

  def preinitialize(self, usePetsc, localComm):
    uefacets.init(usePetsc, localComm)
    self.uedge = uefacets.Uedge()

  def readParams(self, inFile):
    self.uedge.readParams(inFile)
    return 0

  def setMpiComm(self, comm):
    self.uedge.setComm(comm)
    return 0

  def setLogFile(self, logFile):
    self.uedge.setLogFile(logFile)
    return 0
  
  def getRankOfInterface(self, name, rank):
    return 1

  def buildData(self):
    self.uedge.buildData()
    return 0

  def buildUpdaters(self):
    self.uedge.buildUpdaters()
    return 0

  def initialize(self):
    self.uedge.initialize()
    return 0

  def complete(self):
    self.uedge.finalize()
    return 0

  def update(self, tnew):
    return self.uedge.advance(tnew)

  def revert(self):
    self.uedge.reset()
    return 0

  def dumpToFile(self, fname):
    self.uedge.dump(fname)
    return 0 

  def restoreFromFile(self, fname):
    self.uedge.restore(fname)
    return 0

  def get0dDouble(self, varname):
    value = (0, self.uedge.getDouble(varname))
    return value 

  def set0dDouble(self, varname, value):
    self.uedge.setDouble(varname, value)
    return 0

  def get0dInt(self, varname):
    value= (0, self.uedge.getInt(varname))
    return value 

  def set0dInt(self, varname, value):
    self.uedge.setInt(varname, value)
    return 0

  def get0dString(self, varname):
    value = (0, self.uedge.getString(varname))
    return 0

  def set0dString(self, varname, value):
    self.uedge.setString(varname, value)
    return 0

  def get0dDoubleAtLocation(self, varname, location, dim):
    # NOT IMPLEMENTED
    return (1, 0)

  def set0dDoubleAtLocation(self, varname, location, dim, value):
    # NOT IMPLEMENTED
    return 1

  def get0dIntAtLocation(self, varname, location, dim):
    # NOT IMPLEMENTED
    return (1, 0)

  def set0dIntAtLocation(self, varname, location, dim, value):
    # NOT IMPLEMENTED
    return 1

