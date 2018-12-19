# Interface class for UEDGE quantities needed to communicate with wall model
from uedge import *

joule_per_ev = 1.6022e-19

class UEwall:
  def __init__(self):
# create empty set and dget function dictionaries; will register
# all the set functions into this dict with their proper
# string names.  This will get appended to the corresponding uefacets dictionaries
    self.setMap = {}
    self.getMap = {}

  def initialize(self):
    #as default, turn off switches is that allow setting the neutral flux on divertors and walls 
    #   (default is uewall is inactive)
    bbb.isextpltmod=0
    bbb.isextwallmod=0
    #
    # get limits of the single boundary index for the various surface
    # types (inner lower divertor, inner wall, etc)
    # Note most of these are not really settable but are set by sheath physics;
    #  invoking them will evoke a suitable error message
    bbb.getbdyindexlims()
    # self.setMap["foo"] = 12345.    #test
    self.setMap["density_EW_H2p1"]                 = [self.setNi_EW,1,2,1]
    self.setMap["density_EW_H2"]                   = [self.setNg_EW,1,2]
    self.setMap["rotationFrequency_EW_H2p1"]       = [self.setUp_EW,1,2,1]
    self.setMap["temperature_EW_H2p1"]             = [self.setTi_EW,1,2,1]
    self.setMap["temperature_EW_electron"]         = self.setTe_EW

    self.setMap["ptclFlux_EW_H2p1"]             = [self.setNiFlux_EW,1,2,1]
    self.setMap["ptclFlux_EW_H2"]          = [self.setNgFlux_EW,1,2]
    self.setMap["angularMomentumFlux_EW_H2p1"] = [self.setUpFlux_EW,1,2,1]
    self.setMap["energyFlux_EW_H2p1"]             = [self.setTiFlux_EW,1,2,1]
    self.setMap["energyFlux_EW_electron"]         = self.setTeFlux_EW
    self.setMap["energyFlux_EW_H2"]          = [self.setTgFlux_EW,1,2]

    self.getMap["density_EW_H2"]                   = [self.getNg_EW,1,2]
    self.getMap["density_EW_H2p1"]                 = [self.getNi_EW,1,2,1]
    self.getMap["rotationFrequency_EW_H2p1"]       = [self.getUp_EW,1,2,1]
    self.getMap["temperature_EW_H2"]                   = [self.getTg_EW,1,2]
    self.getMap["temperature_EW_H2p1"]             = [self.getTi_EW,1,2,1]
    self.getMap["temperature_EW_electron"]         = self.getTe_EW

    self.getMap["ptclFlux_EW_H2p1"]             = [self.getNiFlux_EW,1,2,1]
    self.getMap["ptclFlux_EW_H2"]          = [self.getNgFlux_EW,1,2]
    self.getMap["angularMomentumFlux_EW_H2p1"] = [self.getUpFlux_EW,1,2,1]
    self.getMap["energyFlux_EW_H2p1"]             = [self.getTiFlux_EW,1,2,1]
    self.getMap["energyFlux_EW_electron"]         = self.getTeFlux_EW
    self.getMap["energyFlux_EW_H2"]          = [self.getTgFlux_EW,1,2]

  # utility functions to turn on/off external (vs. internal) setting of neutral fluxes at divertors, wall
  def setextneutrals_at_div(self):
     bbb.isextpltmod=1
  def setextneutrals_at_wall(self):
     bbb.isextwallmod=1
  def unsetextneutrals_at_div(self):
     bbb.isextpltmod=0
  def unsetextneutrals_at_wall(self):
     bbb.isextwallmod=0
  def setextneutrals(self):
     setextneutrals_at_div(self)
     setextneutrals_at_wall(self)
  def unsetextneutrals(self):
     unsetextneutrals_at_div(self)
     unsetextneutrals_at_wall(self)

  def setNi_EW(self, val, loc_index,spec_index):
      raise "Ni_EW is computed, can't be set"

  def setNg_EW(self, val, loc_index,spec_index):
      raise "Ng_EW inappropriate variable, set gas flux instead"

  def setUp_EW(self, val, loc_index,spec_index):
    raise "Up_EW is computed, can't be set"
  
  def setTi_EW(self, val, loc_index,spec_index):
    raise "Ti_EW is computed, can't be set"
  
  def setTe_EW(self, val,loc_index):
    raise "Te_EW is computed, can't be set"

  def setTg_EW(self, val,loc_index):
    raise "Tg_EW inappropriate variable, set gas energy flux instead"

  def setNiFlux_EW(self, val, loc_index,spec_index):
    raise "NiFlux_EW is computed, can't be set"

  def setNgFlux_EW(self, val, loc_index,spec_index):
    bbb.get_ixiybdy(loc_index,ix,iy,wallname)
    # note uedge variable to be set is flux*area
    if wallname == "innerdiv":
      bbb.fngxextlb[iy,spec_index,0] = val*com.sy(ix,iy)
    if wallname == "outerdiv":
      bbb.fngxextrb[iy,spec_index,0] = val*com.sy(ix,iy)
    if wallname == "outerwall":
      bbb.fngyexto[ix,spec_index] = val*com.sx(ix,iy)
    if wallname == "privwall":
      bbb.fngyexti[ix,spec_index] = val*com.sx(ix,iy)

  def setUpFlux_EW(self, val, loc_index,spec_index):
    raise "UpFlux_EW is computed, can't be set"

  def setTiFlux_EW(self, val, loc_index,spec_index):
    raise "TiFlux_EW is computed, can't be set"

  def setTeFlux_EW(self, val,loc_index):
    raise "TeFlux_EW is computed, can't be set"

  def setTgFlux_EW(self, val, loc_index,spec_index):
    bbb.get_ixiybdy(loc_index,ix,iy,wallname)
    if wallname == "innerdiv":
      print "will set bbb.fegxextlb[iy,spec_index,0] = val*com.sy(ix,iy)"
    if wallname == "outerdiv":
      print "will set bbb.fegxextrb[iy,spec_index,1] = val*com.sy(ix,iy)"
    if wallname == "outerwall":
      print "will set bbb.fegyexto[ix,spec_index,1] = val*com.sx(ix,iy)"
    if wallname == "privwall":
      print "will set bbb.fegyexti[ix,spec_index,1] = val*com.sx(ix,iy)"
    
  def getNi_EW(self, loc_index,spec_index):
    return bbb.get1dpoint(bbb.nis[0:nx+1,0:ny+1,spec_index],loc_index)

  def getNg_EW(self, loc_index,spec_index):
    return bbb.get1dpoint(bbb.ngs[0:nx+1,0:ny+1,spec_index],loc_index)

  def getUp_EW(self, loc_index,spec_index):
    val = bbb.get1dpoint(bbb.ups[0:nx+1,0:ny+1,spec_index],loc_index)
    # convert from uedge's parallel velocity array to toroidal rotation
    # frequency.   
    # convert from uedge's parallel velocity array to toroidal rotation
    rbfbt = bbb.get1dpoint(bbb.rbfbt,loc_index)
    rm = bbb.get1dpoint(com.rm[0:nx+1,0:ny+1,0],loc_index)
    return val*rbfbt/rm

  def getTi_EW(self, loc_index,spec_index):
    return bbb.get1dpoint(bbb.tis[0:nx+1,0:ny+1,spec_index],loc_index)/bbb.qe

  def getTe_EW(self,loc_index):
    return bbb.get1dpoint(bbb.tes,loc_index)/bbb.qe

  def getTg_EW(self, loc_index,spec_index):
    return bbb.get1dpoint(bbb.tg[0:nx+1,0:ny+1,spec_index],loc_index)/bbb.qe

  def getNiFlux_EW(self, loc_index,spec_index):
    ix=0;iy=0;wallname=""
    bbb.get_ixiybdy(loc_index,ix,iy,wallname)
    if (wallname == "outerwall" or wallname == "privwall"):
      val = bbb.get1dpoint(bbb.fniy[0:nx+1,0:ny+1,spec_index],loc_index)
      return val/bbb.get1dpoint(com.sy,loc_index)
    if (wallname == "inner_div" or wallname == "outerdiv"):
      val = bbb.get1dpoint(bbb.fnix[0:nx+1,0:ny+1,spec_index],loc_index)
      return val/bbb.get1dpoint(com.sx,loc_index)

  def getNgFlux_EW(self, loc_index,spec_index):
    ix=0;iy=0;wallname=""
    bbb.get_ixiybdy(loc_index,ix,iy,wallname)
    if (wallname == "outerwall" or wallname == "privwall"):
      val = bbb.get1dpoint(bbb.fngy[0:nx+1,0:ny+1,spec_index],loc_index)
      return val/bbb.get1dpoint(com.sy,loc_index)
    if (wallname == "inner_div" or wallname == "outerdiv"):
      val = bbb.get1dpoint(bbb.fngx[0:nx+1,0:ny+1,spec_index],loc_index)
      return val/bbb.get1dpoint(com.sx,loc_index)

  def getUpFlux_EW(self, loc_index,spec_index):
    ix=0;iy=0;wallname=""
    bbb.get_ixiybdy(loc_index,ix,iy,wallname)
    if (wallname == "outerwall" or wallname == "privwall"):
      val = bbb.get1dpoint(bbb.fmiy[0:nx+1,0:ny+1,spec_index],loc_index)
      return val/bbb.get1dpoint(com.sy,loc_index)
    if (wallname == "inner_div" or wallname == "outerdiv"):
      val = bbb.get1dpoint(bbb.fmix[0:nx+1,0:ny+1,spec_index],loc_index)
      return val/bbb.get1dpoint(com.sx,loc_index)

  def getTiFlux_EW(self, loc_index,spec_index):
    ix=0;iy=0;wallname=""
    bbb.get_ixiybdy(loc_index,ix,iy,wallname)
    if (wallname == "outerwall" or wallname == "privwall"):
      val = bbb.get1dpoint(bbb.feiy[0:nx+1,0:ny+1,spec_index],loc_index)
      return val/bbb.get1dpoint(com.sy,loc_index)
    if (wallname == "inner_div" or wallname == "outerdiv"):
      val = bbb.get1dpoint(bbb.feix[0:nx+1,0:ny+1,spec_index],loc_index)
      return val/bbb.get1dpoint(com.sx,loc_index)

  def getTeFlux_EW(self,loc_index):
    ix=0;iy=0;wallname=""
    bbb.get_ixiybdy(loc_index,ix,iy,wallname)
    if (wallname == "outerwall" or wallname == "privwall"):
      val = bbb.get1dpoint(bbb.feey[0:nx+1,0:ny+1,spec_index],loc_index)
      return val/bbb.get1dpoint(com.sy,loc_index)
    if (wallname == "inner_div" or wallname == "outerdiv"):
      val = bbb.get1dpoint(bbb.feex[0:nx+1,0:ny+1,spec_index],loc_index)
      return val/bbb.get1dpoint(com.sx,loc_index)

  def getTgFlux_EW(self, loc_index,spec_index):
    # For now the neutral energy flux is just taken as the gas temperature times the
    # neutral particle flux
    return self.getNgFlux_EW(loc_index,spec_index)*self.getTg_EW(loc_index,spec_index)* \
           joule_per_ev
        
