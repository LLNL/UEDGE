ppp
{
}

***** ParallelSettings:
OMPParallelPandf1 integer /0/      # [0]: serial pandf1 rhs calc [1] omp parallel pandf1 rhs calc
OMPParallelJac    integer /0/     # [0]: serial jacobian calc [1] omp parallel jacobian calc
ParallelWarning   integer /1/     # Warning for users who wish to use it
CheckJac          integer /0/      # [0/1]: Turn on on-the-fly comparison of parallel vs serial evaluation of Jacobian.
                                  # If differences between para and serial Jacobians, dump both Jacs in serialjac.dat and paralleljac.dat with routine jac_write in current working folder. See UEDGEToolBox docs for analysis tools.
Nthreads          integer /64/ # Number of threads to be used to calculate the Jacobian
CheckPandf1       integer /1/      # [0/1]: Turn on on-the-fly comparison of parallel vs serial evaluation of pandf1.

***** ParallelDebug:
OMPJacDebug       integer /0/ #Print debug info for omp constructs
iidebugprint      integer /-1/ # index ii of jacobian dyldot(ii)/yl(iv) at which threadprivate variables are printed after calculation of the jacobian element. iv is determined by ivdebugprint
ivdebugprint      integer /-1/ # index iv of jacobian dyldot(ii)/yl(iv) at which threadprivate variables are printed after calculation of the jacobian element. ii is determined by iidebugprint
ForceSerialCheck  integer /0/      # [0/1]: Force two sequential serial evaluations of the Jacobian to verify that Jacobian evaluation is reproducible (fail when e.g. variables are not properly initialized in pandf).
DebugJac          integer /0/
DumpJac           integer /0/      # [0/1]: Turn on dumping of data for the diverging element of serial and parallel jacobian (only available when CheckJac is on). See UEDGEToolBox docs for analysis tools.
DumpFullJac       integer /0/      # [0/1]: Turn on dumping of full serial jacobian for analysis of bandwidth (dumping in file ). See UEDGEToolBox docs docs for analysis tools.
WriteJacobian     integer /0/ # Write jacobian in an ascii text file
OMPCopyArray      integer /1/ # For Debug purpose: turn on/off(0/1) copy of threadprivate arrays before jacobian calculation (WARNING:could cause numerical inacurarry if turned on)
OMPCopyScalar     integer /1/ # For Debug purpose: turn on/off copy(0/1) of threadprivate scalar before jacobian calculation (WARNING:could cause numerical inacurarry if turned on)

***** OMPJacSettings:
OMPJacNchunks     integer /0/# Number of chunks to be used to calculate the Jacobian. If Nchunks.lt.0, Nchunks=Nthreads elif Nchunks=0, Nchunks=neq
OMPJacVerbose integer /0/ #Print info for omp jacobian calculation
OMPlenpfac       integer /5/ # Factor to increase nnzmxperchunk
OMPCheckNaN       integer /0/ #Check whether jacobian terms are NaN after jacobian calculation
OMPLoadBalance integer /0/ # Enable user defined weights for each OMP tasks (overrided by MPIAutoBalance)
OMPAutoBalance integer /1/ # Automatic load balancing for OMP thread tasks (if OMPLoadWeight=)
OMPBalanceStrength real /1.0/ # Strenght s of the load balance (Loadweight=Loadweight*(t_thread/<t_thread>)**s)
OMPTimingJacRow    integer /0/ # Profile execution time of calculatation of each row of the jacobian
OMPLoopJacNchunk     integer /1/
OMPJacStamp       character*20 /"*OMPJac* "/ # Stamp for hybrid output (not an user input)

***** OMPJac:
NchunksJac          integer /1/
nnzmxperchunk   integer # Maximum number of jacobian elements which 
                        # can be stored per thread. Can be increased 
                        # with omplenpfac
OMPivmin(NchunksJac)   _integer # jacobian rows with 
                                # ivmin(ithread)<=iv<=ivmax(ithread) are calculated
                                # on thread ithread
OMPivmax(NchunksJac)   _integer # jacobian rows with ivmin(ithread)<=iv<=ivmax(ithread) 
                                # are calculated on thread ithread
OMPLoadWeight(1:NchunksJac)  _real  # weight for load distribution of jacobian 
                                    # calculation among threads
OMPTimeLocalJac(1:NchunksJac)  _real    # runtime for jac calculation on each threads. 
                                        # Used to optimize load distribution of 
                                        # jacobian calculation among threads when 
                                        # AutoBalance=1
iJacRow(neq) _integer  #
OMPTimeJacRow(neq) _real  #
iJacCol(nnzmxperchunk,NchunksJac) _integer #
rJacElem(nnzmxperchunk,NchunksJac) _real #
nnz(NchunksJac) _integer #
nnzcum(NchunksJac) _integer #

**** OMPPandf1Settings:
OMPPandf1Stamp character*20 /"*OMPPandf1* "/ #
OMPPandf1Debug integer /0/
OMPPandf1Verbose integer /0/
OMPTimeParallelPandf1 real /0.0/
OMPTimeSerialPandf1 real /0.0/
OMPPandf1LoopNchunk     integer /1/
OMPPandf1Nychunks integer /0/
OMPPandf1Nxchunks integer /1/
xpadding integer /2/
ypadding integer /2/

**** OMPPandf1:
Nychunks integer /0/
Nxchunks integer /1/
NchunksPandf1 integer /1/
Nychunks_old integer /-1/
Nxchunks_old integer /-1/
yincchunk(NchunksPandf1) _integer
xincchunk(NchunksPandf1) _integer
ixchunk(NchunksPandf1) _integer
iychunk(NchunksPandf1) _integer
Nivchunk(NchunksPandf1) _integer
ivchunk(NchunksPandf1,neq) _integer
iymaxchunk(NchunksPandf1)_integer
ixmaxchunk(NchunksPandf1)_integer
iyminchunk(NchunksPandf1)_integer
ixminchunk(NchunksPandf1)_integer

neq_old integer /0/

**** OMPTiming:
DebugTime  integer         /0/    # Display execution times of various subroutines
ShowTime       integer /1/     # Show execution time of routines
SerialDebug       integer /0/     # Show execution time of routines

#OMPTotJacCalc real  /0./ # time to calculate jacobian in jac_calc_omp
#OMPTotTimeCollect real  /0./ # time to collect jacobian elements in jac_calc_omp
#OMPTotTimeBuild real  /0./ # time to calculate elements of jacobian in jac_calc_omp

***** JacDebug:
EvalDumpJac(FileName:string) subroutine
