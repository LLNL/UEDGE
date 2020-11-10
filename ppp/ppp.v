ppp
{
}

***** ParallelSettings:
OMPParallelPandf1 integer /0/      # [0]: serial pandf1 rhs calc [1] omp parallel pandf1 rhs calc
OMPParallelJac     integer /0/     # [0]: serial jacobian calc [1] omp parallel jacobian calc
MPIParallelJac     integer /0/     # [0]: serial jacobian calc [1] omp parallel jacobian calc
ParallelWarning         integer /1/     # Warning for users who wish to use it

***** OMPSettings:
OMPAllocDebug integer /0/ #Print info on allocation and association of variables
OMPneq  integer # number of equation (=neq)
OMPCopyArray integer /1/ # For Debug purpose: turn on/off(0/1) copy of threadprivate arrays before jacobian calculation (WARNING:could cause numerical inacurarry if turned on)
OMPCopyScalar integer /1/ # For Debug purpose: turn on/off copy(0/1) of threadprivate scalar before jacobian calculation (WARNING:could cause numerical inacurarry if turned on)
Nthreads          integer /64/ # Number of threads to be used to calculate the Jacobian
OMPCopyDebug integer /0/ #Print debug info for omp constructs

***** HybridSettings:
HybridOMPMPI    integer /0/        # Flag for HybridOMPMPI. Automatically turn on with OMPParallelJac and PMIParallelJac
HybridDebug integer /0/ #Print debug info for omp constructs
HybridVerbose integer /1/ #Print info for omp jacobian calculation
HybridCheckNaN       integer /0/ #Check whether jacobian terms are NaN after jacobian calculation
HybridStamp       character*20 # Stamp for hybrid output (not an user input)

***** MPIJacSettings:
MPIRank integer /0/ # Rank of the processor
ComSize integer /1/ # Size of the common world
MPIDebug integer /0/ #Print debug info for omp constructs
MPIVerbose integer /1/ #Print info for omp jacobian calculation
MPIWriteJacobian     integer /0/ # Write jacobian in an ascii text file
MPIlenpfac       integer /1/ # Factor to increase nnzmxperthread
nnzmxperproc   integer # Maximum number of jacobian elements which can be stored per thread. Can be increased with omplenpfac
MPIneq  integer # number of equation (=neq)
Nprocs          integer /64/ # Number of threads to be used to calculate the Jacobian
MPICheckNaN       integer /0/ # Check whether jacobian terms are NaN after jacobian calculation
ioutmpi           integer /6/ # Unit for stdout for common mpi write statements
MPILoadBalance integer /0/ #  User defined weights for each MPI tasks (overrided by MPIAutoBalance)
MPIAutoBalance integer /1/ # Automatic load balancing for MPI tasks
MPIBalanceStrength real /1.0/ # Strenght s of the load balance (Loadweight=Loadweight*(tproc/<tproc>)**s)
MPIStamp          character*20 # Stamp for MPI output (not an user input)

***** MPIJacobian:
MPIivmin(0:Nprocs-1)   _integer # jacobian rows with ivmin(ithread)<=iv<=ivmax(ithread) are calculated on thread ithread (not an user input)
MPIivmax(0:Nprocs-1)   _integer # jacobian rows with ivmin(ithread)<=iv<=ivmax(ithread) are calculated on thread ithread (not an user input)
MPIiJacRow(MPIneq) _integer #
MPIiJacCol(nnzmxperproc) _integer #
MPIrJacElem(nnzmxperproc) _real #
MPILoadWeight(0:Nprocs-1)  _real  # weight for load distribution of jacobian calculation among threads
MPITimeLocalJac(0:Nprocs-1)  _real  # runtime for jac calculation on each threads. Used to optimize load distribution of jacobian calculation among threads when AutoBalance=1


***** OMPJacSettings:
OMPJacDebug integer /0/ #Print debug info for omp constructs
OMPJacVerbose integer /1/ #Print info for omp jacobian calculation
iidebugprint      integer /-1/ # index ii of jacobian dyldot(ii)/yl(iv) at which threadprivate variables are printed after calculation of the jacobian element. iv is determined by ivdebugprint
ivdebugprint      integer /-1/ # index iv of jacobian dyldot(ii)/yl(iv) at which threadprivate variables are printed after calculation of the jacobian element. ii is determined by iidebugprint
WriteJacobian     integer /0/ # Write jacobian in an ascii text file
OMPlenpfac       integer /1/ # Factor to increase nnzmxperthread
nnzmxperthread   integer # Maximum number of jacobian elements which can be stored per thread. Can be increased with omplenpfac
OMPCheckNaN       integer /0/ #Check whether jacobian terms are NaN after jacobian calculation
OMPLoadBalance integer /0/ # Enable user defined weights for each OMP tasks (overrided by MPIAutoBalance)
OMPAutoBalance integer /1/ # Automatic load balancing for OMP thread tasks (if OMPLoadWeight=)
OMPBalanceStrength real /1.0/ # Strenght s of the load balance (Loadweight=Loadweight*(t_thread/<t_thread>)**s)
OMPJacStamp       character*20 /"*OMPJac* "/ # Stamp for hybrid output (not an user input)
OMPTimingJacRow    integer /0/ # Profile execution time of calculatation of each row of the jacobian
DebugJac          integer /0/
CheckJac          integer /0/      # [0/1]: Turn on on-the-fly comparison of parallel vs serial evaluation of Jacobian.
                                  # If differences between para and serial Jacobians, dump both Jacs in serialjac.dat and paralleljac.dat with routine jac_write in current working folder. See UEDGEToolBox docs for analysis tools.
DumpJac           integer /0/      # [0/1]: Turn on dumping of data for the diverging element of serial and parallel jacobian (only available when CheckJac is on). See UEDGEToolBox docs for analysis tools.
DumpFullJac  integer /0/      # [0/1]: Turn on dumping of full serial jacobian for analysis of bandwidth (dumping in file ). See UEDGEToolBox docs docs for analysis tools.
ForceSerialCheck  integer /0/      # [0/1]: Force two sequential serial evaluations of the Jacobian to verify that Jacobian evaluation is reproducible (fail when e.g. variables are not properly initialized in pandf).

***** OMPJacobian:
OMPivmin(Nthreads)   _integer # jacobian rows with ivmin(ithread)<=iv<=ivmax(ithread) are calculated on thread ithread
OMPivmax(Nthreads)   _integer # jacobian rows with ivmin(ithread)<=iv<=ivmax(ithread) are calculated on thread ithread
OMPLoadWeight(1:Nthreads)  _real  # weight for load distribution of jacobian calculation among threads
OMPTimeLocalJac(1:Nthreads)  _real  # runtime for jac calculation on each threads. Used to optimize load distribution of jacobian calculation among threads when AutoBalance=1
iJacRow(OMPneq,Nthreads) _integer  #
OMPTimeJacRow(OMPneq) _real  #
iJacCol(nnzmxperthread,Nthreads) _integer #
rJacElem(nnzmxperthread,Nthreads) _real #
nnz(Nthreads) _integer
nnzcum(Nthreads) _integer

**** OMPPandf1Settings:
OMPPandf1Stamp character*20 /"*OMPPandf1* "/ # Stamp for hybrid output (not an user input)
OMPPandf1Debug integer /0/
OMPPandf1Verbose integer /0/
OMPPandf1FlagVerbose integer /0/
OMPPandf1Check integer /1/
NthreadsPandf1 integer /8/ # Nthreads for threaded pandf1.
OMPPandf1yinc integer /3/
OMPPandf1padyinc integer /2/
OMPTimeParallelPandf1 real /0.0/
OMPTimeSerialPandf1 real /0.0/
OMPPandf1RunPara integer /1/
OMPPandf1FirstRun integer /1/

**** OMPPandf1:
OMPic(1:Nthreads) _integer
OMPyinc(1:Nthreads) _integer
OMPivthread(1:OMPneq) _integer
OMPTimeLocalPandf1(1:Nthreads) _real #
OMPTimeCollectPandf1(1:Nthreads) _real #

**** OMPTiming:
DebugTime  integer         /0/    # Display execution times of various subroutines
ShowTime       integer /1/     # Show execution time of routines
SerialDebug       integer /0/     # Show execution time of routines

OMPTotJacCalc real  /0./ # time to calculate jacobian in jac_calc_omp
OMPTotTimeCollect real  /0./ # time to collect jacobian elements in jac_calc_omp
OMPTotTimeBuild real  /0./ # time to calculate elements of jacobian in jac_calc_omp

**** MPITiming:
MPITotTimeCollect real  /0./ # time to collect jacobian elements in jac_calc_mpi/hybrid
MPITotTimeBuild real  /0./ # time to calculate elements of jacobian in jac_calc_mpi/hybrid
MPITotJacCalc real  /0./ # time to calculate jacobian in jac_calc_mpi/jac_calc_hybrid



#TimeSerialPandf real /0.0/
#TimeParallelPandf real /0.0/
#TimingParaPandf integer /0/
#TimingConvert integer /0/
#TimingParaConvert integer /0/
#OMPTimePandf real /0.0/
#OMPPandfVerbose integer /0/
#RhsEval integer /0/
#TimeConvert0 real /0/
#TimeConvert1 real /0/
#TimeNeudif real /0/
#TimeBlock1 real /0/
#TimeBlock2 real /0/
#TimeBlock3 real /0/
#TimeBlock4 real /0/
#TimeBlock5 real /0/
#TimeBlock6 real /0/
#TimeBlock7 real /0/
#TimeBlock8 real /0/
#TimeBlock9 real /0/
#TimeBlock10 real /0/
#TimeBlock11 real /0/
#TimeBlock12 real /0/
#TimeBlock13 real /0/
#TimeBlock14 real /0/
#TimeBlock15 real /0/
#TimeBlock16 real /0/
#TimeConv0 real /0/
#TimeConv1 real /0/
#TimeParaConv1 real /0/
#TimeConv2 real /0/
#TimeConv3 real /0/
#TimeConv4 real /0/
#TimeConv5 real /0/
#TimeConv6 real /0/
#TimeParaBlock1 real /0/
#TimeParaBlock2 real /0/
#TimeParaBlock3 real /0/
#TimeParaBlock4 real /0/
#TimeParaBlock5 real /0/
#TimeParaBlock6 real /0/
#TimeParaBlock7 real /0/
#TimeParaBlock8 real /0/
#TimeParaBlock9 real /0/
#TimeParaBlock10 real /0/
#TimeParaBlock11 real /0/
#TimeParaBlock12 real /0/
#TimeParaBlock13 real /0/
#TimeParalock14 real /0/
#TimeParaBlock15 real /0/
#TimeParaBlock16 real /0/
#TimeCopy0 real /0/
#TimeCopy1 real /0/
#TimeCopy2 real /0/

***** JacDebug:
EvalDumpJac(FileName:string) subroutine
