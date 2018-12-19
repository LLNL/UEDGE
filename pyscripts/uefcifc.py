#
# @file        uefcifc.py
#
# @brief       Interface for facets interface coupling
#
# @version     $Id: uefcifc.py,v 7.0 2018/02/28 18:43:48 meyer8 Exp $
#

# required imports
import uefacets

# This class defines interfaces to Uedge in a fashion that closely
# resembles the FACETS interface in order to facilitate easy coupling
class UedgeFacetsIfc:

    # Constructor.  Functions should get registered in here so as to
    # ensure that the set and get maps are available for the complete
    # lifetime of the UedgeFacetsIfc object.  Internal variables to
    # the ifc class should also be declared here.
    def __init__(self):
        self.ueObj = uefacets.Uedge()

        # create an empty setter function dictionary, and register
        # all the set functions into this dict with their proper
        # string names
        
        self.setMap = {}

        self.setMap["init_file"]    = self.setInFileName
        self.setMap["restore_file"] = self.setRestoreName

        self.setMap["density_CE_H2+1"]                 = self.setNi
        self.setMap["density_CE_neutral"]              = self.setNg
        self.setMap["parallelMomentumDensity_CE_H2+1"] = self.setUp
        self.setMap["temperature_CE_H2+1"]             = self.setTi
        self.setMap["temperature_CE_electron"]         = self.setTe

        # create an empty getter function dictionary and register
        # all the get functions into this dict with their proper
        # string names.
        
        self.getMap = {}

        self.getMap["temperature_CE_electron"]         = self.getTe
        self.getMap["density_CE_neutral"]              = self.getNg
        self.getMap["density_CE_H2+1"]                 = self.getNi
        self.getMap["temperature_CE_H2+1"]             = self.getTi
        self.getMap["parallelMomentumDensity_CE_H2+1"] = self.getUp

    # this constructs the data required for the solver to function
    # no assumptions to the actual solver can be made here, because it
    # is not guaranteed that buildSolver will have been called before
    # this.
    def buildData(self):
        pass

    # this constructs the solver itself.  Assumptions about data can be
    # made, since buildData is guaranteed to be called before this.
    # However, no assumptions can be made about the state of the data, as
    # it is not guaranteed that initialize will have been called before
    # this
    def buildSolver(self):
        pass

    # this initializes the solver and prepares for runs.  This can safely
    # assume that the data has been completely constructed and that the
    # solver is constructed and ready for use.  This is called last in
    # the initialization callpath.
    def initialize(self):
        self.ueObj.readParams(self.infile)
        self.ueObj.restore(self.rsfile)

    # this advances UEDGE to time t
    def advance(self, t):
        self.ueObj.advance(t)

    # this finalizes the UEDGE object.  Note that this call should NOT
    # destroy data.  Instead, all postprocessing of data and solver state
    # should be performed here.  Data destruction should be left up to
    # the class destructor.
    def finalize(self):
        pass

    # set data used by UEDGE.  This function should stay as is, all
    # implementation should be done in the functions registered in
    # the map
    def setData(self, name, val):
        self.setMap[name](name, val)

    # get data used by UEDGE.  This function should stay as is, all
    # implementation should be done in the functions registered in
    # the map
    def getData(self, name):
        return self.getMap[name](name)

    # set the file name from which to read parameters
    def setInFileName(self, name, fname):
        self.infile = fname

    # set the file name from which to restore
    def setRestoreName(self, name, rname):
        self.rsfile = rname

    # set parameters in the UEDGE object.  Any data processing and
    # massaging that needs to be done should be performed in these
    # functions.
    def setNi(self, name, val):
        n = len(val)
        bbb.isnicore = 1
        bbb.ncore[0:n] = val

    def setNg(self, name, val):
        n = len(val)
        bbb.isngcore = 1
        bbb.ngcore[0:n] = val

    def setUp(self, name, val):
        n = len(val)
        bbb.isupcore = 0
        bbb.upcore[0:n] = val
  
    def setTi(self, name, val):
        n = len(val)
        bbb.iflcore = 0
        bbb.tcorei[0:n] = val * bbb.qe

    def setTe(self, name, val):
        n = len(val)
        bbb.iflcore = 0
        bbb.tcoree[0:n] = val * bbb.qe

    # get parameters from the UEDGE object.  Any data processing and
    # massaging that needs to be done should be performed in these
    # functions
    def getNi(self, name):
        return self.ueObj.getNi()

    def getNg(self, name):
        return self.ueObj.getNg()

    def getUp(self, name):
        return self.ueObj.getUp()

    def getTe(self, name):
        return self.ueObj.getTe()

    def getTi(self, name):
        return self.ueObj.getTi()
