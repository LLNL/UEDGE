#
# $Id: uefacets.py,v 7.0 2018/02/28 18:43:48 meyer8 Exp $
# Implementation of uedge methods for FACETS
#
# Substantial revisions R. Cohen 4/21/08
#
# Caveat 1 for now: We are not distinguishing between parallel and
# toroidal momentum fluxes at the interface.  (This is better than
# previously when the toroidal and parallel velocities themselves
# were treated as indistinguishable.)
# Doing better will require further work to ensure consistency between
# the core and edge models.
# Caveat 2 for now: no coupling of the potential equation.  Closely
# related to caveat 1.

import __main__
import string
import types
from copy import copy
from uedge import *
from uewall import *

joule_per_ev = 1.6022e-19

# Import pytables
try:
#SEK  import tables
  # bbb.writeToLog( "tables imported." )
  # from tables.nodes import FileNode
  # from tables import FileNode
  # bbb.writeToLog( "FileNode imported." )
  import numpy as numerix
  have_tables = 1
except ImportError, inst:
  # bbb.writeToLog( "ImportError:", inst )
  have_tables = 0
except:
  # bbb.writeToLog( "Exception:", sys.exc_info() )
  have_tables = 0

class UedgeError:

  def __init__(self, msg):
    self.msg = msg

  def __str__(self):
    return self.msg

iset_petsc_local=0
iset_mpi_local=0
def init(localPetsc=True, localComm=True):
# Initialize MPI and PETSc
    global iset_petsc_local
    global iset_mpi_local
    flx.init_hdf5()
    if (localComm):
       bbb.uedge_mpiInit()
       iset_mpi_local=1
       # if not bbb.mype:
         # bbb.writeToLog( "NOTE: localComm is True on rank %d.  Should not be a component."%bbb.mype )
    if (localPetsc):
       bbb.uedge_petscInit()
       iset_petsc_local=1
       # if not bbb.mype:
         # bbb.writeToLog( "PETSc initialized on rank %d."%bbb.mype )

def final():
# Finalize petsc and MPI
    global iset_petsc_local
    global iset_mpi_local
    if iset_petsc_local:
       if not bbb.mype:
         bbb.writeToLog( "uedge_petscFinal called on rank %d."%bbb.mype )
       bbb.uedge_petscFinal()
    if iset_mpi_local:
       if not bbb.mype:
         bbb.writeToLog( "uedge_mpiFinal called on rank %d."%bbb.mype )
       bbb.uedge_mpiFinal()

class Uedge:

  def __init__(self):
    """Create an instance of UEDGE that can be executed through Python
    This will register the dictionaries at the top level.
    """
# Initialize the time, step, ...
    self.t = 0.
    self.t_old = 0.
    self.initcall = 1

# create an empty setter function dictionary; will register
# all the set functions into this dict with their proper
# string names
    self.setMap = {}
    self.setMap["nspecies_CE"] = 0 # = self.setNSpecies

# create an empty getter function dictionary; will register
# all the get functions into this dict with their proper
# string names.
    self.getMap = {}

# Get all of the dictionaries registered
    __main__.__dict__['bbb'] = bbb
    __main__.__dict__['aph'] = aph
    __main__.__dict__['com'] = com
    __main__.__dict__['flx'] = flx
    __main__.__dict__['api'] = api
    __main__.__dict__['grd'] = grd
    __main__.__dict__['svr'] = svr
    __main__.__dict__['wdf'] = wdf

# For coupled runs, turn off uedge's core recycling.  Eventually
# we should have a way of telling UEDGE whether we are running
# coupled or not.
    bbb.recycc = 0.

    # Create an instance of the UEDGE wall interface class
    self.uewall = UEwall()

# Set up dictionary of methods
    # density setters and flux setters will be handled with different
    # string signatures.  density_* will be for density (below) and
    # ptclFlux_* will be for fluxes
    # Note RHS is a list [method, atomic number, atomic mass, charge state]
    self.setMap["density_CE_H2p1"]                 = [self.setNi_CE,1,2,1]
    self.setMap["density_CE_H2"]                   = [self.setNg_CE,1,2]
    self.setMap["rotationFrequency_CE_H2p1"]       = [self.setUp_CE,1,2,1]
    self.setMap["temperature_CE_H2p1"]             = [self.setTi_CE,1,2,1]
    self.setMap["temperature_CE_electron"]         = self.setTe_CE
    self.setMap["temperature_CE_H2"]               = [self.setTg_CE,1,2]

    self.setMap["ptclFlux_CE_H2p1"]                = [self.setNiFlux_CE,1,2,1]
    self.setMap["ptclFlux_CE_H2"]                  = [self.setNgFlux_CE,1,2]
    self.setMap["angularMomentumFlux_CE_H2p1"]     = [self.setUpFlux_CE,1,2,1]
    self.setMap["energyFlux_CE_H2p1"]              = [self.setTiFlux_CE,1,2,1]
    self.setMap["energyFlux_CE_electron"]          = self.setTeFlux_CE
    self.setMap["energyFlux_CE_H2"]                = [self.setTgFlux_CE,1,2]

    self.getMap["density_CE_H2"]                   = [self.getNg_CE,1,2]
    self.getMap["density_CE_H2p1"]                 = [self.getNi_CE,1,2,1]
    self.getMap["rotationFrequency_CE_H2p1"]       = [self.getUp_CE,1,2,1]
    self.getMap["temperature_CE_H2p1"]             = [self.getTi_CE,1,2,1]
    self.getMap["temperature_CE_electron"]         = self.getTe_CE
    self.getMap["temperature_CE_H2"]               = [self.getTg_CE,1,2]

    self.getMap["ptclFlux_CE_H2p1"]                = [self.getNiFlux_CE,1,2,1]
    self.getMap["ptclFlux_CE_H2"]                  = [self.getNgFlux_CE,1,2]
    self.getMap["angularMomentumFlux_CE_H2p1"]     = [self.getUpFlux_CE,1,2,1]
    self.getMap["energyFlux_CE_H2p1"]              = [self.getTiFlux_CE,1,2,1]
    self.getMap["energyFlux_CE_electron"]          = self.getTeFlux_CE
    self.getMap["energyFlux_CE_H2"]                = [self.getTgFlux_CE,1,2]

    # getters for various diffuson coefficients
    self.getMap["densityDiffusion_CE_H2p1"]          = [self.getDensDif_CE,1,2,1]
    self.getMap["densityConvection_CE_H2p1"]         = [self.getDensConv_CE,1,2,1]
    self.getMap["parmomentumDiffusion_CE_H2p1"]      = [self.getParmomDif_CE,1,2,1]
    self.getMap["parmomentumConvection_CE_H2p1"]     = [self.getParmomConv_CE,1,2,1]
    self.getMap["temperatureDiffusion_CE_H2p1"]      = [self.getTiDif_CE,1,2,1]
    self.getMap["temperatureConvection_CE_H2p1"]     = [self.getTiConv_CE,1,2,1]
    self.getMap["temperatureDiffusion_CE_electron"]  = self.getTeDif_CE
    self.getMap["temperatureConvection_CE_electron"] = self.getTeConv_CE

    # map entires for provenance data
    self.setMap["dataDirectory"]                     = self.setShareDir

    # add edge-wall interface quantities to get, setMap
    self.uewall.initialize()

  def readParams(self, sourcefile):
    """Read the parameters from the file, sourcefile.
    This file should set parameters only.  It should not call any methods.
    """
    if not bbb.mype:
      # bbb.writeToLog( "Uedge.readParams called with file, %s."%sourcefile )
      pass
    try:
      execfile(sourcefile)
    except:
      raise
    if not bbb.mype:
      bbb.writeToLog( "Uedge.readParams: %s sourced."%sourcefile )
    # if petsc enabled - readin petsc params [options] aswell.
    bbb.uedge_petscInsertOpts()

  def setLogFile(self, logFile):
    bbb.setLogFile(logFile)

  def setComm(self, mpiComm=None):
# Eventually: will set up MPI
    bbb.set_uedgeComm(mpiComm)

  def buildData(self):
    """buildData should now allocate all of the data.
    """
    if not bbb.mype:
      # bbb.writeToLog( "Uedge.buildData: entered." )
      pass
    grd.rplate1 = [ 0 for i in range(0, grd.nplate1) ]
    grd.zplate1 = [ 0 for i in range(0, grd.nplate1) ]
    grd.rplate2 = [ 0 for i in range(0, grd.nplate2) ]
    grd.zplate2 = [ 0 for i in range(0, grd.nplate2) ]

    bbb.allocate()
    if bbb.ismpion == 1:
      bbb.init_par_meshg()
    if not bbb.mype:
      bbb.writeToLog( "Uedge.buildData: allocation complete." )

  def buildUpdaters(self):
    # this constructs the solver itself.  Assumptions about data can be
    # made, since buildData is guaranteed to be called before this.
    # However, no assumptions can be made about the state of the data, as
    # it is not guaranteed that initialize will have been called before
    # this
    pass

  def fluxSurfAv_CE(self,InArray):
    # utility function to take flux surface averages on core-edge interface.
    # InArray is either a 1-D array defined on UEDGE's core poloidal array, of
    # length com.ixpt2-com.ixpt1, or a full 2-D array on UEDGE's full grid

    ixpt11 = com.ixpt1[0]+1; ixpt21 = com.ixpt2[0]+1
    if (len(shape(InArray)) == 1):
        return sum(com.volv[ixpt11:ixpt21,0]*InArray)/   \
               sum(com.vol[ixpt11:ixpt21,0])
    elif (len(shape(InArray)) == 2):
        return sum(com.volv[ixpt11:ixpt21,0]*InArray[ixpt11:ixpt21,0])/   \
          sum(com.vol[ixpt11:ixpt21,0])
    else:
        raise "Improper dimension for InArray"


  # set boundary parameters in the UEDGE object.  Any data processing and
  # massaging that needs to be done should be performed in these
  # functions.

  def gradAv_CE(self,InArray):
    # utility function to take gradient-weighted averages over a flux surface
    # in other words, returns Int(inArray*(d/dr))/Int(d/dr).
    # InArray must be 2D array defined over UEDGE's full grid.
    ixpt11 = com.ixpt1[0]+1; ixpt21 = com.ixpt2[0]+1
    return sum(com.gy[ixpt11:ixpt21,0]*com.sy[ixpt11:ixpt21,0] * \
	InArray[ixpt11:ixpt21,0]) / sum(com.gy[ixpt11:ixpt21,0] * \
	com.sy[ixpt11:ixpt21,0])

  def initialize(self):
    # Initialize the solver and prepare for runs.
    # This assumes that the data has been completely constructed and that the
    # solver is constructed and ready for use.  This is called last in
    # the initialization callpath.
    # Per facets rules, initialize should contain only things that are
    # used for a new run but not a restart.
    #
    self.commonsetup(restart=0)

  def commonsetup(self,restart):
    # Set up mass and charge arrays and flux surface average arrays needed
    # for the variable setting methods.  Must have arrays allocated and so must
    # be called after buildData.  Must know whether this is a restart,
    # to set bbb.restart in order to be able to call this from
    # initialize() or from restore()
    bbb.restart=restart

    if bbb.ismpion and restart == 1:
    # if bbb.ismpion:	# Hangs
      bbb.writeToLog( "WARNING: For parallel run, averages not computed by commonsetup." )
      return

    bbb.ueinit()

    # ion atomic mass and charge state lists
    self.imasslist = bbb.minu.round(0).tolist()
    self.chargelist = bbb.ziin.round(0).tolist()
    # neutral atomic mass list, construct from unique members of imasslist.
    # this will fail (produce run-time error) if com.ngsp is bigger than
    # the length of this list; but it is also a failure of the physics
    # specification.
    self.nmasslist = uniquelist(self.imasslist)

    self.ixpt11 = com.ixpt1 + 1
    self.ixpt21 = com.ixpt2 + 1

    ixpt11 = self.ixpt11
    ixpt21 = self.ixpt21
    self.Raverage = bbb.fluxsurfav1(com.rm[ixpt11:ixpt21,0,0])
    # In above, making use of thin-ness of guard cell so  needn't worry
    #  abut using centered cell value of r.

  def finalize(self):
    # this finalizes the UEDGE object.  Note that this call should NOT
    # destroy data.  Instead, all postprocessing of data and solver state
    # should be performed here.  Data destruction should be left up to
    # the class destructor.
    pass

  def getDouble(self, name):
    # get data used by UEDGE.  This function should stay as is, all
    # implementation should be done in the functions registered in
    # the map.  Note this can be used by UEDGE experts to get any variable
    # in UEDGE's variable descriptor files, if name is of the form
    # packagename.variablename
    self.nisp=com.nisp
    self.ngsp=com.ngsp
    self.nisp1=self.nisp+1
    self.nusp=com.nusp

    try:
      rhs = self.getMap[name]
      if type(rhs)==types.ListType:
        # ions or neutrals
        if len(rhs) == 4:
          # ions
          return rhs[0](self.findindexi(rhs[1:]))
        if len(rhs) == 3:
          # neutrals
          return rhs[0](self.findindexn(rhs[1:]))
      else:
          # electrons
          return rhs()
    except:
      # The expert option.
      print eval(name)
      return eval(name)
      # An exception will be raised if name is invalid

  def setDouble(self, name, value):
    # set data used by UEDGE.  This function should stay as is, all
    # implementation should be done in the functions registered in
    # the map.  Note this can be used by UEDGE experts to set any variable
    # in UEDGE's variable descriptor files, if name is of the form
    # packagename.variablename

    #bbb.writeToLog( "Inside uefacets::setDouble with name ", name, " value ", value )

    try:
      rhs = self.setMap[name]
      if type(rhs)==types.ListType:
        # ions or neutrals
        if len(rhs) == 4:
          rhs[0](value,self.findindexi(rhs[1:]))
        if len(rhs) == 3:
          # neutrals
          rhs[0](value,self.findindexn(rhs[1:]))
      else:
        # electrons
        rhs(value)
    except:
      # The expert option.
      bbb.writeToLog( "Using expert option to set a parameter" )
      bbb.writeToLog( "EXPERT OPTION WITH NAME '%s'" % name )
      exec(name+"="+`value`)
      # An exception will be raised if name is invalid

  def setString(self, name, value):
    # set data used by UEDGE.  This function should stay as is, all
    # implementation should be done in the functions registered in
    # the map.  Note this can be used by UEDGE experts to set any variable
    # in UEDGE's variable descriptor files, if name is of the form
    # packagename.variablename
    print name, value
    # self.setMap[name](value)

  def getString(self, name):
    # get data used by UEDGE.  This function should stay as is, all
    # implementation should be done in the functions registered in
    # the map.  Note this can be used by UEDGE experts to get any variable
    # in UEDGE's variable descriptor files, if name is of the form
    # packagename.variablename
    self.getMap[name](value)

  def getDoubleAtIndex(self, name,ndim,loc_index):
    # get data used by UEDGE.  This function should stay as is, all
    # implementation should be done in the functions registered in
    # the map.  It is assumed that any invocation of getDoubleAtIndex
    # is targeted at quantities for communication with the wall,
    # so name has a "_EW" attached.
    if name.find("_EW") < 0:
      raise name+" is invalid edge-wall interface name"
    self.nisp=com.nisp
    self.ngsp=com.ngsp
    self.nisp1=self.nisp+1
    self.nusp=com.nusp

    rhs = self.uewall.getMap[name]
    if type(rhs)==types.ListType:
      # ions or neutrals
      if len(rhs) == 4:
        # ions
        return rhs[0](loc_index,self.findindexi(rhs[1:]))
      if len(rhs) == 3:
        # neutrals
        return rhs[0](loc_index,self.findindexn(rhs[1:]))
    else:
      # electrons
      return rhs(loc_index)

  def setDoubleAtIndex(self, name, ndim,loc_index,value):
    # set wall-relevant data used by UEDGE.  This function should stay as is,
    # all implementation should be done in the functions registered in
    # the map, which are in the file uewall.py.
    # It is assumed that any invocation of setDoubleAtIndex is targeted at
    # quantities for communication with the wall, so name has a "_EW" attached.
    #
    if name.find("_EW") < 0:
      raise name+" is invalid edge-wall interface name"
   #bbb.writeToLog( "Inside uefacets::setDoubleAtIndex with name ", name, " value ", value )
    rhs =self.uewall.setMap[name]
    if type(rhs)==types.ListType:
      # ions or neutrals
      if len(rhs) == 4:
        rhs[0](value,loc_index,self.findindexi(rhs[1:]))
      if len(rhs) == 3:
        # neutrals
        rhs[0](value,loc_index,self.findindexn(rhs[1:]))
    else:
      # electrons
      rhs(value,loc_index)


  def advance(self, tnew):
    # sys.stderr.write("Uedge.advance: entered.\n")
    # tentatively advances state from current time to time + dt = tnew, and
    # does NOT return the updated values on the boundary of the quantities
    # that were not input.
    # To proceed to a next step, call ueAdvance again.
    # To retake this step, call ueReset and then this method again.
    # variables returned are those from the end of the step, in sequence:
    #  t,nivar,ngvar,upvar,tivar,tevar
    if not bbb.mype:
      # bbb.writeToLog("Uedge.advance: entered." )
      pass
    self.t_old = self.t
    self.t = tnew
    dt = self.t - self.t_old
    dt_min = 1.e-12

    bbb.dtreal = max(dt, dt_min)
    bbb.dt_tot = 0.
    bbb.t_stop = bbb.dtreal
    # bbb.iprtrundt = 1
    if (self.initcall == 1):
      if not bbb.mype:
        bbb.writeToLog("Calling bbb.exmain.")
      # sys.stderr.write("Uedge.advance: calling bbb.exmain.\n")
      bbb.exmain()
      # sys.stderr.write("Uedge.advance: bbb.exmain returned.\n")
    else:
      if not bbb.mype:
        bbb.writeToLog( "Calling bbb.rundt." )
      bbb.rundt()
    self.initcall = 0

    if not bbb.iterm==1:
      return 1
    else:
      return 0

  def interprett(self):
    # reads plasma profile data from DIII-D in the form of tanh fits for
    # ne (=ni) and te, and spline fits for ti. Data files must be called
    # netanh1.dat, tetanh1.dat, tislpline1.dat. Then calculates transport
    # coefficients dif_use, kye_use, and kyi_use (for particle D, Chi_e,
    # and Chi_i) using the uedge interpretive mode; valid only in core
    # with values outside taken as continuation of separatrix values. Leaves
    # calc transport coeff and plasma profiles updated for further simulation.
    bbb.interptrans()

  def reset(self):
    # resets Uedge variables to their state prior to last call to ueAdvance.
    # returns time at beginning of last call to ueAdvance
    bbb.nis=copy(bbb.ni0)
    bbb.ups=copy(bbb.up0)
    bbb.ngs=copy(bbb.ng0)
    bbb.tes=copy(bbb.te0)
    bbb.tis=copy(bbb.ti0)
    self.t = self.t_old
    # return self.t_old
    return

  def setShareDir(self, val):
    aph.sharedir = val

  def setNi_CE(self, val, index):
    bbb.isnicore = 1
    bbb.ncore[index] = val

  def setNg_CE(self, val, index):
    bbb.isngcore = 1
    bbb.ngcore[index] = val

  def setUp_CE(self, val, index):
    bbb.isupcore = 5
    # Facets agrees to use toroidal angular velocity; UEDGE isupcore=5
    # assumes ang. velocity is const but inut quantity is flux surface
    # av of vtor.   Connection is <vtor> = Omega*<R>, and
    # <R> = (Integral dArea/Bpol)/(Integral dArea/(RBol))
    # uedge's upcore does not.
    #bbb.upcore[index] = val*self.Raverage
    bbb.utorave = val*self.Raverage

  def setTi_CE(self, val, index):
    bbb.iflcore = 0
    # eventually may want to set an array of ion temperatures
    bbb.tcorei = val

  def setTe_CE(self, val):
    bbb.iflcore = 0
    bbb.tcoree = val

  def setTg_CE(self, val, index):
    raise "Setting of gas temperature not currently implemented in UEDGE"

  def setNiFlux_CE(self, val, index):
    bbb.isnicore = 3
    bbb.curcore[index] = val*bbb.qe*com.area_core

  def setNgFlux_CE(self, val, index):
    raise "Gas flux b.c. from core is not implemented in UEDGE"

  def setUpFlux_CE(self, val, index):
    bbb.isupcore = 4
    # The following assumes that FACETS adopts UEDGE's convention:
    # The flux used on the interface is the flux of toroidal angular momentum.
    # Unlike uedge's other flux inputs, this one is average per unit area, not
    # integrated over the flux surface
    bbb.lzflux[index] = val

  def setTiFlux_CE(self, val, index):
    bbb.iflcore = 1
    # eventually may want to set an array of ion heat fluxes
    bbb.pcorei = val * com.area_core
    if not bbb.mype:
      bbb.writeToLog( "setTiFlux_CE: val = " + `val` + ", global_area = " + `com.area_core` + ", global_flux = " + `bbb.pcorei` )

  def setTeFlux_CE(self, val):
    bbb.iflcore = 1
    bbb.pcoree = val * com.area_core
    if not bbb.mype:
      bbb.writeToLog( "setTeFlux_CE: val = " + `val` + ", global_area = " + `com.area_core` + ", global_flux = " + `bbb.pcoree` )


  def setTgFlux_CE(self, val, index):
    raise "Gas heat flux b.c. not currenetly implemented in UEDGE"

  # get boundary parameters in the UEDGE object.  Any data processing and
  # massaging that needs to be done should be performed in these
  # functions.

  def getNi_CE(self, index):
    ns_CEBdry = bbb.fluxsurfav1(bbb.nis[self.ixpt11:self.ixpt21,0,index])
    return ns_CEBdry

  def getNg_CE(self, index):
    nsgas_CEBdry = bbb.fluxsurfav1(bbb.ngs[self.ixpt11:self.ixpt21,0,index])
    return nsgas_CEBdry

  def getUp_CE(self, index):
    # convert from uedge's parallel velocity array to toroidal rotation
    # frequency.   And then do a flux surface average.
    ixpt11=self.ixpt11;ixpt21=self.ixpt21
    OmegaTor = bbb.ups[ixpt11:ixpt21,0,index]*bbb.rbfbt[ixpt11:ixpt21,0]/  \
               com.rm[ixpt11:ixpt21,0,0]
    return bbb.fluxsurfav1(OmegaTor)

  def getTi_CE(self, index):
    ti_CEBdry = bbb.fluxsurfav2(bbb.tis)/bbb.qe
    return ti_CEBdry

  def getTe_CE(self):
    te_CEBdry = bbb.fluxsurfav2(bbb.tes)/bbb.qe
    return te_CEBdry

  def getTg_CE(self, index):
    tsgas_CEBdry = bbb.fluxsurfav1(bbb.tg[self.ixpt11:self.ixpt21,0,index])
    return tsgas_CEBdry

  def getNiFlux_CE(self, index):

    nsFlux_CEBdry = sum(bbb.fniy[self.ixpt11:self.ixpt21,0,index])/   \
                    com.area_core
    return nsFlux_CEBdry

  def getNgFlux_CE(self, index):
    ngFlux_CEBdry = sum(bbb.fngy[self.ixpt11:self.ixpt21,0,index])/   \
                    com.area_core
    return ngFlux_CEBdry

  def getUpFlux_CE(self, index):
      # z angular momentum density flux
      ixpt11=self.ixpt11;ixpt21=self.ixpt21
      # fmity appears to be broken, so use fmiy for now, ignoring distinction
      # between toroidal momentum and parallel momentum
      v_parFlux_CEBdry = sum(bbb.fmiy[ixpt11:ixpt21,0,index]*com.rm[ixpt11:ixpt21,0,0])/ \
                         com.area_core
      return v_parFlux_CEBdry

  def getTiFlux_CE(self, index):
    # eventually may want to get an array of ion heat fluxes
      tsFlux_CEBdry = sum(bbb.feiy[self.ixpt11:self.ixpt21,0])/(com.area_core)
      return tsFlux_CEBdry


  def getTeFlux_CE(self):
      tsFlux_CEBdry = sum(bbb.feey[self.ixpt11:self.ixpt21,0])/(com.area_core)
      return tsFlux_CEBdry

  def getTgFlux_CE(self, index):
    # For now the neutral energy flux is just taken as the gas temperature
    # times the neutral particle flux
    return self.getNgFlux_CE(index)*self.getTg_CE(index)*joule_per_ev

  def getDensDif_CE(self, index):
    return self.gradAv_CE(bbb.fcdif*bbb.difni[index]+bbb.dif_use[:,:,index])
  def getDensConv_CE(self, index):
    return self.gradAv_CE(bbb.vy_use[:,:,index])
  def getParmomDif_CE(self, index):
    return self.gradAv_CE(bbb.fcdif*bbb.travis[index]+bbb.tra_use[:,:,index])
  def getParmomConv_CE(self, index):
    return self.gradAv_CE(bbb.vyup_use)
  def getTiDif_CE(self, index):
    return self.gradAv_CE(bbb.fcdif*bbb.kyi+bbb.kyi_use)
  def getTiConv_CE(self, index):
    return self.gradAv_CE(bbb.vyti_use)
  def getTeDif_CE(self):
    return self.gradAv_CE(bbb.fcdif*bbb.kye+bbb.kye_use)
  def getTeConv_CE(self):
    return self.gradAv_CE(bbb.vyte_use)

  # Utility routines to fetch the 2D+species arrays of the variables
  def getNi2D(self, name):
    return bbb.nis
  def getNg2D(self, name):
    return bbb.ngs
  def getUp2D(self, name):
    return bbb.ups
  def getTe2D(self, name):
    return bbb.tes
  def getTi2D(self, name):
    return bbb.tis

  """
  Dump the mesh and variables to a file.
  The type of file is determined by the suffix:
    pdb for PDB, h5 for HDF5, hdf for HDF4.
  If hdf5, the output is vizschema compliant.
  """
  def dump(self, fname):
    # bbb.writeToLog( "Entered dump with fname =", fname )

# Dump to hdf5
    if fname[len(fname) - 3:] == ".h5":
      doTransposeFlag=1
      flx.dump_vshdf5(fname, self.t, doTransposeFlag)
      bbb.writeToLog( "Dumped to hdf5" )
      return

# Dump to hdf4
    if fname[len(fname) - 4:] == ".hdf":
      dumpedobjs = Forthon.pydump(fname, hdf=1, returnfobjlist=1)
      bbb.writeToLog( "Dumped " + `dumpedobjs` + "to hdf" )
      return

# Default: dump to pdb
    dumpedobjs = Forthon.pydump(fname, hdf=0, returnfobjlist=1)
    bbb.writeToLog( "Dumped " + `dumpedobjs` + "to pdb" )

  """
  Restore the variables (not the mesh) from a file.
  The type of file is determined by the suffix:
    pdb for PDB, h5 for HDF5, hdf for HDF4.
  """
  def restore(self, fname):
    if not bbb.mype:
      # bbb.writeToLog( "Uedge.restore: entered." )
      pass
    # Before we do anything, restore from python file if specified
    pyrestart=bbb.pyrestart_file[0].strip()
    if pyrestart != "":
      if os.path.exists(pyrestart):
        bbb.writeToLog( "Including py file in restart:  " + pyrestart )
        execfile(pyrestart)
      else:
        bbb.writeToLog( "ERROR: "+ pyrestart + " does not exist." )

    try:
# Restore from hdf5
      if fname[len(fname) - 3:] == ".h5":
        # raise UedgeError("In restore.  have_tables =" + `have_tables`)
        # msg = "In restore in uefacets.py to open"
        if not bbb.mype:
          bbb.writeToLog( "Uedge.restore: Opening %s on rank %d."%(fname,bbb.mype) )
          bbb.writeToLog( "time = " + `self.t` )
        # For python, I want to do the transpose
        doTransposeFlag=1
        errval = 0
        flx.restore_vshdf5(fname, doTransposeFlag, errval)
        if ( errval != 0 ):
          return 1
        if not bbb.mype:
          pass
          # bbb.writeToLog( "Uedge.restore: flx.restore_vshdf5 returned." )

# Restore from hdf4
      elif fname[len(fname) - 4:] == ".hdf":
        bbb.writeToLog( "Restoring following variables from hdf" )
        dumpedobjs = Forthon.restore(fname, ls=1)

# Default: restore from pdb
      else:
        msg = "Restoring from pdb"
        bbb.writeToLog( msg )
        # raise UedgeError(msg)
        Forthon.restore(fname, ls=0)

# Iterate once to a solution
#       bbb.exmain()

    except string, ex:
      ex += ".  File name = fname."
      raise ex

# set bbb.restart=1 and set up the flux surface average,
# mass and charge arrays for name setting
# For hdf5 files, we handle the commonsetup.
    if not fname[len(fname) - 3:] == ".h5":
         self.commonsetup(restart=1)
    else:
         # ion atomic mass and charge state lists
         self.imasslist = bbb.minu.round(0).tolist()
         self.chargelist = bbb.ziin.round(0).tolist()
         # neutral atomic mass list, construct from unique members of imasslist.
         # this will fail (produce run-time error) if com.ngsp is bigger than
         # the length of this list; but it is also a failure of the physics
         # specification.
         self.nmasslist = uniquelist(self.imasslist)

    self.ixpt11 = com.ixpt1 + 1
    self.ixpt21 = com.ixpt2 + 1

    ixpt11 = self.ixpt11;ixpt21=self.ixpt21
    # bbb.writeToLog( "sy =", com.sy[ixpt11:ixpt21,0] )
    bbb.writeToLog( "area = "+ `com.area_core` )
    self.Raverage = bbb.fluxsurfav1(com.rm[ixpt11:ixpt21,0,0])
    # In above, making use of thin-ness of guard cell so needn't worry
    # about using centered cell value of r.

  def DoCommand(self,command):
    # Executes an arbitrary python command expressed as a string.
    # This allows skilled UEDGE users to have full access to uedge's
    # python interface via the FACETS interface
    exec(command) in __main__.__dict__

  def ListPackages(self):
    # returns a list of uedge's packages.   Variable names within
    # a package are requested as packagename.variablename
    return package()

  def ListVariables(self,packagename):
    # returns a list of the settable variables in the package whose name
    # ispackagename.
    # if packagename is an empty string then a list of settable variables in
    # all packages is returned
    if packagename != "":
      exec("vars = "+packagename+".varlist()")
      return vars
    else:
      packages = package()
      varlist = []
      for packagename in packages:
        exec("vars = "+packagename+".varlist()")
        varlist = varlist+vars
      return varlist

  def ListVarsAndVals(self,packagename):
    # returns a dictionary (pairs of variable names and values) for package
    # packagename.  In contrast to ListVariables, a packagename must be
    # specified or an exception is thrown.
    exec("varsvals = "+packagename+".getdict()")
    return varsvals

  def Dump(self,fnameroot):
    #dump variables to both HDF and PDB files, for now, until we find a good HDF reader for Python.
    # file named fname.  Files are named by appending .hdf or .pdb to fnameroot
    # For now use pydump until we figure out where uedge's dump method
    # comes from.
    pdbname = fnameroot+".pdb"
    hdfname = fnameroot+".hdf"
    pydump(hdfname,hdf=0)
    pydump(pdbname,hdf=1)

  def findindexi(self,params):
    # Returns uedge plasma species index corresponding to atomicmass and charge.
    # params is a list containing atomicnumber,atomicmass,charge
    # Requires masslist and chargelist to be set up as lists of masses and
    # charge states for each species index to be set up in advance of calling
    # For now, assuming species uniquely identified by mass and charge state.
    # Really we should allow for isotopes of different elements having same
    # atomic mass; this will come later.
    #
    # find all indices with specified mass
    # number of entries with proper mass
    atomicnumber=params[0]
    mass=params[1]
    charge=params[2]
    masslist=self.imasslist
    chargelist=self.chargelist
    numberwithmass = masslist.count(mass)
    massindices = []
    istart=0
    for i in range(numberwithmass):
        imass=masslist[istart:].index(mass)+istart
        massindices += [imass]
        istart=imass+1
    for ind in massindices:
        if chargelist[ind] == charge:
            return ind
    # if we've gotten this far, our mass, charge combo is not in list
    raise "Index error, specified mass and charge not in species list"

  def findindexn(self,params):
    # Find index for neutrals.  for now, this is just a lookup by atomic mass
    # params is a list containing atomicnumber,atomicmass
    return self.nmasslist.index(params[1])

def restore(fname):
  ueobj = Uedge()
  ueobj.restore(fname)

def dump(fname):
  ueobj = Uedge()
  ueobj.dump(fname)

def uniquelist(list1):
  # returns a list of the unique elements in list
  templist=sorted(list1)
  lenlist=len(templist)
  for i in range(lenlist-1):
    index=lenlist-i-1
    if templist[index-1]==templist[index]:
      del templist[index]
  return templist

def set_uedgeComm(communicator):
  ireturn = bbb.set_uedgeComm(communicator)

