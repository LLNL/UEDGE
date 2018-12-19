"""
This script contains all the mechanisms required for the Uedge HDF5
conversion/annotation script to function.  The class defined within
(UeConvertVsH5) will take an HDF5 file and annotate it or create a new
file that is vsh5 compliant.
"""

# Import all required libraries
import numpy
import tables

import ueplotdata

class UeConvertVsH5:
    """
    The purpose of this class is to convert a specified uedge file
    into the new vsh5 schema system.  Note that the old file is
    modified instead of a new file being created.  This is because the
    conversion script actually simply adds new annotation to the old
    HDF5 file.  The reasoning for this is for backwards compatibility
    with other scripts that have not yet been converted to read the
    new vizschema.
    """

    def __init__(self, fname=''):
        """
        Constructor, takes no action, unless an input file is
        specified.  If an input file is specified then that file will
        be annotated automatically.  Otherwise files can be annotated
        individually

        *fname: optional filename on which to perform processing
        """

        self.fname = fname
        if ( fname != '' ):
            self.readData(fname)
            self.createUnstructGrid()
            self.annotateFile(fname)

    def readData(self, fname=''):
        """
        This function reads data in from the specified filename.  It
        stores the raw data in preparation for converstion and output
        into the new vsh5 visualization format

        *fname: Optional name of the file to be annotated.  If
         unspecified it is assumed it was set in the constructor
        """
        
        if ( fname == '' ):
            fname = self.fname
        if ( fname == '' ):
            raise ValueError( 'UeConvertVsH5::readData: ' +
                              'file name was unspecified' )

        self.infile = fname
        fh = tables.openFile(fname, 'r')

        # Read the Uedge grid information
        uen = fh.getNode('/bbb/uegrid')
        self.ueg = uen.uegrid.read()
        self.psi = uen.psiVals.read()

        # Read the Uedge variable data
        uen = fh.getNode('/bbb')
        self.nis = uen.nis.read()
        self.ngs = uen.ngs.read()
        self.ups = uen.ups.read()
        self.tes = uen.tes.read()
        self.tis = uen.tis.read()
        self.phis = uen.phis.read()

        fh.close()

    def annotateFile(self, fname=''):
        """
        This function will write data out to the specified filename.
        If no filename is specified then it assumes that data should
        be annotated to the already existing file attached to this
        object.

        *fname: Optional name of the file to be annotated.  If
         unspecified it is assumed it was set previously
        """

        if ( fname == '' ):
            fname = self.fname
        if ( fname == '' ):
            raise ValueError( 'UeConvertVsH5::annotateFile: ' +
                              'file name was unspecified' )
        fh = tables.openFile(fname, 'a')

        # Write the unstructured data to the file
        gn = fh.createGroup( '/', 'unstructGrid', 
                             'Uedge grid annotated for vsh5' )
        fh.setNodeAttr(gn, 'vsType', 'mesh')
        fh.setNodeAttr(gn, 'vsKind', 'unstructured')
        fh.createArray(gn, 'points', self.nodes)
        fh.createArray(gn, 'cells', self.cells)

        nds = []
        # Write all the actual data to the file now
        nds.append( fh.createArray('/', 'nis', self.newNis) )
        nds.append( fh.createArray('/', 'ngs', self.newNgs) )
        nds.append( fh.createArray('/', 'ups', self.newUps) )
        nds.append( fh.createArray('/', 'tes', self.newTes) )
        nds.append( fh.createArray('/', 'tis', self.newTis) )
        nds.append( fh.createArray('/', 'phis', self.newPhis) )
        nds.append( fh.createArray('/', 'psi', self.newPsi) )

        for i in nds:
            fh.setNodeAttr( i, 'vsType',   'variable'     )
            fh.setNodeAttr( i, 'vsMesh',   'unstructGrid' )
            fh.setNodeAttr( i, 'vsCellOffset', 'center'   )

        fh.close()

    def annotateStructured(self, fname=''):
        """
        This function is almost identical to the normal annotateFile
        routine.  However, this one write the data in a 3 part block
        structured format instead of fully unstructured.

        *fname: optional name of the file to be annotated.  If
         unspecified it is assumed that it was set previously
        """

        if ( fname == '' ):
            fname = self.fname
        if ( fname == '' ):
            raise ValueError( 'UeConvertVsH5::annotateStructured: ' +
                              'file name was unspecified' )

        pdata = ueplotdata.UePlotData(self.infile)

        fh = tables.openFile(fname, 'a')

        # Write the structured grid data to the file
        gnr = fh.createGroup( '/', 'structuredData',
                             'Uedge grids annotated for structured vsh5' )

        xx, yy = pdata.getNodalCoreGrid()
        arr = numpy.zeros( (xx.shape[0], xx.shape[1], 2), numpy.float )
        arr[:,:,0] = xx
        arr[:,:,1] = yy
        gn = fh.createGroup( gnr, 'coreGrid',
                             'Uedge core section structured grid' )
        fh.setNodeAttr(gn, 'vsType', 'mesh')
        fh.setNodeAttr(gn, 'vsKind', 'structured')
        fh.createArray(gn, 'points', arr)

        xx, yy = pdata.getNodalSOLGrid()
        arr = numpy.zeros( (xx.shape[0], xx.shape[1], 2), numpy.float )
        arr[:,:,0] = xx
        arr[:,:,1] = yy
        gn = fh.createGroup( gnr, 'solGrid',
                             'Uedge scrape-off layer section grid' )
        fh.setNodeAttr(gn, 'vsType', 'mesh')
        fh.setNodeAttr(gn, 'vsKind', 'structured')
        fh.createArray(gn, 'points', arr)

        xx, yy = pdata.getNodalPFGrid()
        arr = numpy.zeros( (xx.shape[0], xx.shape[1], 2), numpy.float )
        arr[:,:,0] = xx
        arr[:,:,1] = yy
        gn = fh.createGroup( gnr, 'pfGrid',
                             'Uedge private flux layer section grid' )
        fh.setNodeAttr(gn, 'vsType', 'mesh')
        fh.setNodeAttr(gn, 'vsKind', 'structured')
        fh.createArray(gn, 'points', arr)

        nds = []
        ars = ['nis', 'ngs', 'ups', 'tes', 'tis', 'phis']
        for i in ars:
            gn = fh.createArray( gnr, i + 'Core', pdata.getCenteredCoreData(i) )
            fh.setNodeAttr( gn, 'vsType',' variable' )
            fh.setNodeAttr( gn, 'vsMesh', 'coreGrid' )
            fh.setNodeAttr( gn, 'vsCellOffset', 'center' )

            gn = fh.createArray( gnr, i + 'Sol', pdata.getCenteredSOLData(i) )
            fh.setNodeAttr( gn, 'vsType', 'variable' )
            fh.setNodeAttr( gn, 'vsMesh', 'solGrid' )
            fh.setNodeAttr( gn, 'vsCellOffset', 'center' )

            gn = fh.createArray( gnr, i + 'Pf', pdata.getCenteredPFData(i) )
            fh.setNodeAttr( gn, 'vsType', 'variable' )
            fh.setNodeAttr( gn, 'vsMesh', 'pfGrid' )
            fh.setNodeAttr( gn, 'vsCellOffset', 'center' )

        gn = fh.createArray( gnr, 'psiCore', pdata.getCenteredCorePsi() )
        fh.setNodeAttr( gn, 'vsType', 'variable' )
        fh.setNodeAttr( gn, 'vsMesh', 'coreGrid' )
        fh.setNodeAttr( gn, 'vsCellOffset', 'center' )

        gn = fh.createArray( gnr, 'psiSol', pdata.getCenteredSOLPsi() )
        fh.setNodeAttr( gn, 'vsType', 'variable' )
        fh.setNodeAttr( gn, 'vsMesh', 'solGrid' )
        fh.setNodeAttr( gn, 'vsCellOffset', 'center' )

        gn = fh.createArray( gnr, 'psiPf', pdata.getCenteredPFPsi() )
        fh.setNodeAttr( gn, 'vsType', 'variable' )
        fh.setNodeAttr( gn, 'vsMesh', 'pfGrid' )
        fh.setNodeAttr( gn, 'vsCellOffset', 'center' )

        fh.close()

    def createUnstructGrid(self):
        """
        The mapped rectangular Uedge grid gets converted into a
        full-on unstructured grid here.  This is necessary because the
        branch cut mapped rectangular grid of Uedge is incompatible
        with vsh5, however it can be made compliant as an unstruct
        quadrilateral mesh.

        returns: A pair of numpy arrays the first of which is the
        nodelist and the second of which is the unstructured cell list
        """

        nodes = []
        cells = []
        
        for i in range(0, self.psi.shape[0]):
            for j in range(0, self.psi.shape[1]):

                narr = range(0, 5)
                q = 0
                for k in [1, 2, 4, 3]:
                    q = q+1
                    inst = [ self.ueg[i, j, k, 0], self.ueg[i, j, k, 1] ]
                    if ( nodes.count(inst) == 0 ):
                        narr[q] = nodes.__len__()
                        nodes.append(inst)
                    else:
                        narr[q] = nodes.index(inst)

                narr[0] = 4
                cells.append(narr)

        self.nodes = numpy.array(nodes, numpy.float)
        self.cells = numpy.array(cells, numpy.int)

        self.newNis  = numpy.zeros( (self.cells.shape[0], 2), numpy.float )
        self.newNgs  = numpy.zeros( (self.cells.shape[0], 1), numpy.float )
        self.newUps  = numpy.zeros( (self.cells.shape[0], 2), numpy.float )
        self.newTes  = numpy.zeros( (self.cells.shape[0], 1), numpy.float )
        self.newTis  = numpy.zeros( (self.cells.shape[0], 1), numpy.float )
        self.newPhis = numpy.zeros( (self.cells.shape[0], 1), numpy.float )
        self.newPsi  = numpy.zeros( (self.cells.shape[0], 1), numpy.float )
        q = 0
        for i in range(0, self.psi.shape[0]):
            for j in range(0, self.psi.shape[1]):
                self.newNis[q]  = self.nis[i, j]
                self.newNgs[q]  = self.ngs[i, j]
                self.newUps[q]  = self.ups[i, j]
                self.newTes[q]  = self.tes[i, j]
                self.newTis[q]  = self.tis[i, j]
                self.newPhis[q] = self.phis[i, j]
                self.newPsi[q]  = self.psi[i, j, 0]
                q = q+1

        self.nisWMesh  = numpy.zeros( (self.cells.shape[0], 4), numpy.float )
        self.ngsWMesh  = numpy.zeros( (self.cells.shape[0], 3), numpy.float )
        self.upsWMesh  = numpy.zeros( (self.cells.shape[0], 4), numpy.float )
        self.tesWMesh  = numpy.zeros( (self.cells.shape[0], 3), numpy.float )
        self.tisWMesh  = numpy.zeros( (self.cells.shape[0], 3), numpy.float )
        self.phisWMesh = numpy.zeros( (self.cells.shape[0], 3), numpy.float )
        self.psiWMesh  = numpy.zeros( (self.cells.shape[0], 3), numpy.float )
        q = 0;
        for i in range(0, self.psi.shape[0]):
            for j in range(0, self.psi.shape[1]):
                self.nisWMesh[q, :] = numpy.array( [ self.ueg[i, j, 0, 0],
                                                     self.ueg[i, j, 0, 1],
                                                     self.nis[i, j, 0],
                                                     self.nis[i, j, 1 ] ],
                                                   numpy.float )
                self.ngsWMesh[q, :] = numpy.array( [ self.ueg[i, j, 0, 0],
                                                     self.ueg[i, j, 0, 1],
                                                     self.ngs[i, j] ],
                                                   numpy.float )
                self.upsWMesh[q, :] = numpy.array( [ self.ueg[i, j, 0, 0],
                                                     self.ueg[i, j, 0, 1],
                                                     self.ups[i, j, 0],
                                                     self.ups[i, j, 1 ] ],
                                                   numpy.float )
                self.tesWMesh[q, :] = numpy.array( [ self.ueg[i, j, 0, 0],
                                                     self.ueg[i, j, 0, 1],
                                                     self.tes[i, j] ],
                                                   numpy.float )
                self.tisWMesh[q, :] = numpy.array( [ self.ueg[i, j, 0, 0],
                                                     self.ueg[i, j, 0, 1],
                                                     self.tis[i, j] ],
                                                   numpy.float )
                self.phisWMesh[q, :] = numpy.array( [ self.ueg[i, j, 0, 0],
                                                      self.ueg[i, j, 0, 1],
                                                      self.phis[i, j] ],
                                                    numpy.float )
                self.psiWMesh[q, :] = numpy.array( [ self.ueg[i, j, 0, 0],
                                                     self.ueg[i, j, 0, 1],
                                                     self.psi[i, j, 0] ],
                                                   numpy.float )
                q = q+1

        return self.nodes, self.cells
