#include "../bbb/bbb.h"
#include "../com/com.h"
!#include "../mppl.h"
!#include "../sptodp.h"
!-------------------------------------------------------------------------------------------------
subroutine InitParallel

    Use ParallelSettings,only:OMPParallelJac,MPIParallelJac,ParallelJac,OMPParallelPandf1
    Use HybridSettings,only:HybridOMPMPI
    implicit none
    if ((OMPParallelJac==1 .or. OMPParallelPandf1.gt.0) .and. MPIParallelJac==0) then
        call InitOMP
        ParallelJac=1
    elseif (OMPParallelJac==0 .and. MPIParallelJac==1) then
        call InitMPI
        ParallelJac=1
    elseif (OMPParallelJac==1 .and. MPIParallelJac==1) then
        HybridOMPMPI=1
        call InitMPI() ! MPI first to setup nnzmxperproc used by InitOMP
        call InitOMP()
        ParallelJac=1
    else
        write(*,'(a)') '*OMPJac* Jacobian calculation: OMP not enabled'
        write(*,'(a)') '*MPIJac* Jacobian calculation: MPI not enabled'
        ParallelJac=0
    endif

end subroutine InitParallel
!-------------------------------------------------------------------------------------------------
subroutine jac_calc_parallel(neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)


    Use ParallelSettings,only:OMPParallelJac,MPIParallelJac
    Use OMPJacSettings,only:iidebugprint,ivdebugprint,DebugJac,ForceSerialCheck,CheckJac,DumpFullJac, DumpJac
    Use Cdv,only:exmain_aborted
    implicit none
    ! ... Input arguments:
    integer,intent(in):: neq      !      total number of equations (all grid points)
    real,intent(in)   :: t              ! physical time
    real,intent(inout)   :: yl(*)          ! dependent variables
    real,intent(in)   :: yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu         ! lower and upper bandwidths
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian

    ! ... Output arguments:
    real,intent(out)   :: jac(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: ja(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: ia(neq+1)   ! pointers to beginning of each row in jac,ja
    real   :: jaccopy(nnzmx)     ! nonzero Jacobian elements
    integer:: jacopy(nnzmx)   ! col indices of nonzero Jacobian elements
    integer:: iacopy(neq+1)   ! pointers to beginning of each row in jac,ja
    integer i,iv,t_start
    real,external::tock
    ! ... Work-array argument:
    real wk(neq)     ! work space available to this subroutine
    if (exmain_aborted) call xerrab('exmain aborted...')
    if (OMPParallelJac==1 .and. MPIParallelJac==0) then
        call jac_calc_omp (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
    elseif (OMPParallelJac==0 .and. MPIParallelJac==1) then
        call jac_calc_mpi (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
    elseif (OMPParallelJac==1 .and. MPIParallelJac==1) then
        call jac_calc_hybrid (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
    else
        call xerrab('Cannot call calc_jac_parallel OMP and MPI jac calc are not enabled')
    endif

    if (DumpFullJac.gt.0) then
    write(*,*) 'Dumping full serial jacobian for analysis of bandwidth'
    call jac_calc (neq, t, yl, yldot00, neq, neq, wk,neq*neq, jaccopy, jacopy, iacopy)
    call jac_write('serialfulljac.dat',neq, jaccopy, jacopy, iacopy)
    endif

    if (ForceSerialCheck.gt.0) then
    write(*,*) 'Force check of serial evaluation of jacobian'
    write(*,*) '----- Performing first serial Evaluation of jacobian...'
    call jac_calc (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
    write(*,*) '----- Performing second serial Evaluation of jacobian...'
      call jac_calc (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jaccopy, jacopy, iacopy)
     write(*,*) '----- Comparing jacobians...: nzmax=',ia(neq)-1
      do i=1,ia(neq)-1
       if (abs(jaccopy(i)-jac(i)).gt.1e-14) then
           write(*,*) ' ****** diff between jacs :',i,jac(i),jaccopy(i),'|',ja(i),jacopy(i)
           do iv=1,neq
            if (i<iacopy(iv)) then
            write(*,*) 'idx row serial:',iv-1
            exit
            endif
           enddo
           do iv=1,neq
            if (i<ia(iv)) then
            write(*,*) 'idx row parallel:',iv-1
            exit
            endif
           enddo
        call xerrab('Parallel evaluation of Jacobian differs from serial evaluation... ')
        exit
      endif
      enddo

    endif

    if (CheckJac.gt.0) then


      write(*,*) '--- Checking if parallel and serial evaluations of jacobian are the same...'
       write(*,*) '----- Performing serial Evaluation of jacobian...'
       call tick(t_start)
      call jac_calc (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jaccopy, jacopy, iacopy)
      write(*,*) '----- Serial Evaluation of jacobian done in: ',tock(t_start) ,' s'
     write(*,*) '----- Writing jacobian into serialjac.dat and paralleljac.dat'
      call jac_write('paralleljac.dat',neq, jac, ja, ia)
      call jac_write('serialjac.dat',neq, jaccopy, jacopy, iacopy)
      write(*,*) '----- Comparing jacobians...: nzmax=',ia(neq)-1
      do i=1,ia(neq)-1
       if (abs(jaccopy(i)-jac(i)).gt.1e-14) then
           write(*,*) ' ****** diff between jacs :',i,jac(i),jaccopy(i),'|',ja(i),jacopy(i)
           do iv=1,neq
            if (i<iacopy(iv)) then
            write(*,*) 'idx row serial:',iv-1
            exit
            endif
           enddo
           do iv=1,neq
            if (i<ia(iv)) then
            write(*,*) 'idx row parallel:',iv-1
            exit
            endif
           enddo
           if (DumpJac.gt.0) then
           iidebugprint=ja(i)
           ivdebugprint=iv-1
           write(*,*) 'recalculating jacobian to dump debug data for iv=',ivdebugprint,' ii=',iidebugprint
           if (OMPParallelJac==1 .and. MPIParallelJac==0) then
            call jac_calc_omp (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
           elseif (OMPParallelJac==0 .and. MPIParallelJac==1) then
            call jac_calc_mpi (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
           elseif (OMPParallelJac==1 .and. MPIParallelJac==1) then
            call jac_calc_hybrid (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
           else
            call xerrab('Cannot call calc_jac_parallel OMP and MPI jac calc are not enabled')
           endif
           call jac_calc (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jaccopy, jacopy, iacopy)
           endif
           call xerrab('Parallel evaluation of Jacobian differs from serial evaluation... ')
           exit
       endif
      enddo

      endif

end subroutine jac_calc_parallel
!-------------------------------------------------------------------------------------------------
#ifdef _OPENMP
!subroutine InitPandfOMP()
!use OMPJacSettings,only: Nthreads
!use OMPPandf
!use Dim,only:ny
!OMPPandf_Nthreads=min(Nthreads,ny)
!OMPPandf_Nthreads=min(OMPPandf_Nthreads,OMPPandfNthreads)
!if (OMPPandf_Nthreads<1) call xerrab('OMPPandf_Nthreads<1')
!if (OMPParallelPandf.gt.0) then
!   write(*,*) 'Parallel RHS eval Pandf activated: Pandf_Nthreads=',OMPPandf_Nthreads
!endif
!call gchange('OMPPandf',0)
!if (OMPThreadedPandfngxflux.lt.1 .and. OMPThreadedPandfngyflux.gt.0) then
!call xerrab('OMPThreadedPandfngyflux cannot be on with OMPThreadedPandfngxflux')
!endif
!
!if ((OMPThreadedPandfngxflux.lt.1 .or. OMPThreadedPandfngyflux.lt.1) .and. OMPThreadedPandfngxyflux.gt.0) then
!call xerrab('OMPThreadedPandfngxyflux cannot be on with OMPThreadedPandfngxflux and OMPThreadedPandfngyflux')
!endif
!
!end subroutine InitPandfOMP

subroutine InitOMP()

    Use OMPJacSettings,only:nnzmxperthread,omplenpfac,OMPJacVerbose,OMPJacStamp
    Use OMPSettings,only:Nthreads,ompneq
    Use ParallelSettings,only: OMPParallelPandf1
    Use OMPPandf1Settings,only: OMPPandf1FlagVerbose,OMPPandf1FirstRun
    Use MPIJacSettings,only:nnzmxperproc,MPIrank
    Use HybridSettings,only: HybridOMPMPI,Hybridstamp
    Use Jacobian,only:nnzmx
    Use Lsode, only:neq

    implicit none
    integer:: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,OMP_GET_MAX_THREADS
    character*8 :: MPIRankTag

    ! prepare MPI/Hybrid stamps for output
    if (HybridOMPMPI>0) then
        write(MPIRankTag,'(I4)') MPIrank
        write(OMPJacStamp,'(a,a,a)') '[',trim(adjustl(trim(MPIRankTag))),'] OMPJac* '
    else
        !write(OMPJacStamp,'(a)') '*OMPJac* '
    endif
    if (OMPJacVerbose.gt.1) write(*,*) OMPJacStamp,' Max number of threads available:',OMP_GET_MAX_THREADS()
    call OMP_SET_NUM_THREADS(OMP_GET_MAX_THREADS())

    !$omp parallel
    if (OMP_GET_THREAD_NUM().eq.0) then
        if (Nthreads.gt.OMP_GET_MAX_THREADS()) then
            if (OMPJacVerbose.gt.1) write(*,*) OMPJacStamp,' Warning: # threads requested > # threads available'

            Nthreads=OMP_GET_NUM_THREADS()

            write(*,*) OMPJacStamp,' Nthreads:', Nthreads
        endif
        if (Nthreads.le.0) then
            call xerrab('Nthread must be >0')
        endif
        if (OMPJacVerbose.gt.0) write(*,'(a,a,i3)') OMPJacStamp,' Number of threads for omp calculations:',Nthreads
    endif
    !$omp END parallel

    if (Nthreads.gt.1) then
        if (HybridOMPMPI>0) then
            nnzmxperthread=ceiling(real(nnzmxperproc)/real(Nthreads-1)*omplenpfac)
        else
            nnzmxperthread=ceiling(real(nnzmx)/real(Nthreads-1)*omplenpfac)
        endif
    else
        if (HybridOMPMPI>0) then
            nnzmxperthread=nnzmxperproc
        else
            nnzmxperthread=nnzmx
        endif
    endif
    ompneq=neq
    call gchange('OMPJacobian',0)
    if (OMPParallelPandf1.gt.0) then
      OMPPandf1FlagVerbose=1
      OMPPandf1FirstRun=1
    endif

end subroutine InitOMP
!-------------------------------------------------------------------------------------------------
subroutine OMPCollectJacobian(neq,nnzmx,rcsc,icsc,jcsc,nnzcumout)

    use OMPJacobian,only:iJacCol,rJacElem,iJacRow,OMPivmin,OMPivmax,nnz,nnzcum,OMPTimeJacRow
    use OMPJacSettings,only:OMPJacDebug,OMPJacVerbose,OMPJacStamp,OMPTimingJacRow
    use OMPSettings,only:Nthreads
    integer,intent(in):: neq
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian
    real,intent(out)   :: rcsc(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: icsc(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: jcsc(neq+1)   ! pointers to beginning of each row in jac,ja
    integer,intent(out):: nnzcumout
    integer ith
    integer:: iunit,iv
    nnzcum(1:Nthreads)=-1
    nnzcum(1)=nnz(1)-1

    do ith=2,Nthreads
        nnzcum(ith)=nnzcum(ith-1)+nnz(ith)-1
    enddo
    if (OMPJacDebug.gt.0) then
        write(*,*) OMPJacStamp,' nnz:',nnz(1:Nthreads)
        write(*,*) OMPJacStamp,' nnzcum:',nnzcum(1:Nthreads)
    endif
    if (OMPJacVerbose.gt.0) write(*,'(a,i9)') '**** Number of non-zero Jacobian elems:',nnzcum(Nthreads)

    if (nnzcum(Nthreads).gt.nnzmx) then
        write(*,*) 'nnzcum=',nnzcum
        write(*,*) nnzmx
        call xerrab(' Problem: nnzcum > nnzmx...')
    endif



    jcsc(OMPivmin(1):OMPivmax(1))= iJacRow(OMPivmin(1):OMPivmax(1),1)
    do ith=2,Nthreads
        jcsc(OMPivmin(ith):OMPivmax(ith))= iJacRow(OMPivmin(ith):OMPivmax(ith),ith)+nnzcum(ith-1)
    enddo

    rcsc(1:nnz(1)-1)= rJacElem(1:nnz(1)-1,1)
    icsc(1:nnz(1)-1)= iJacCol(1:nnz(1)-1,1)
    do ith=2,Nthreads
        rcsc(nnzcum(ith-1)+1:nnzcum(ith))=rJacElem(1:nnz(ith)-1,ith)
        icsc(nnzcum(ith-1)+1:nnzcum(ith))=iJacCol(1:nnz(ith)-1,ith)
    enddo
    nnzcumout=nnzcum(Nthreads)
    if (OMPTimingJacRow.gt.0) then
    open (newunit = iunit, file = 'omptiming.dat')
    do iv=1,neq
    write(iunit,*) iv,OMPTimeJacRow(iv)
    enddo
    close(iunit)
    endif
end subroutine OMPCOllectJacobian
!-------------------------------------------------------------------------------------------------
subroutine jac_calc_omp (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)

    !   Calculate Jacobian matrix (derivatives with respect to each
    !   dependent variable of the right-hand side of each rate equation).
    !   Lower and upper bandwidths are used to select for computation
    !   only those Jacobian elements that may be nonzero.
    !   Estimates of Jacobian elements are computed by finite differences.
    !   The Jacobian is stored in compressed sparse row format.

    use Timing,only:istimingon,ttjstor,ttotjf,ttimpjf
    use OMPTiming,only:OMPTotTimeCollect,OMPTotTimeBuild,OMPTotJacCalc
    use Grid,only:ngrid,ig,ijac,ijactot
    use Jacobian_csc,only:rcsc,jcsc,icsc,yldot_pert,yldot_unpt
    use OMPSettings,only:Nthreads
    use OMPJacSettings,only:OMPJacDebug,OMPJacVerbose,OMPJacDebug,nnzmxperthread,OMPCheckNaN,WriteJacobian,&
        OMPLoadBalance,OMPAutoBalance,OMPJacStamp,OMPBalanceStrength
    use OMPJacobian,only:iJacCol,rJacElem,iJacRow,OMPivmin,OMPivmax,nnz,nnzcum,OMPLoadWeight,OMPTimeLocalJac
    use UEpar, only: svrpkg
    use Math_problem_size,only:neqmx
    implicit none
    ! ... Input arguments:
    integer,intent(in):: neq      !      total number of equations (all grid points)
    real,intent(in)   :: t              ! physical time
    real,intent(in)   :: yl(*)          ! dependent variables
    real,intent(in)   :: yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu         ! lower and upper bandwidths
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian

    ! ... Output arguments:
    real,intent(out)   :: jac(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: ja(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: ia(neq+1)   ! pointers to beginning of each row in jac,ja

    ! ... Work-array argument:
    real wk(neq)     ! work space available to this subroutine
    integer,allocatable :: iJacConstructor(:,:)
    real,allocatable:: rJacConstructor(:,:)
    integer:: nnzcumout
    ! ... Functions
    logical tstguardc
    real(kind=4) gettime

    ! ... Local variables:
    real(kind=4) sec4, tsjstor, tsimpjf, dtimpjf,time0,time1,OMPTimeBuild,OMPTimeCollect,OMPTimeJacCalc
    integer:: i,thread,ith,iv,TID, OMP_GET_THREAD_NUM
    character(len = 80) ::  filename
    ! Calculate load distribution for threads
    if (OMPLoadBalance.ne.1 .and. OMPAutoBalance.ne.1) then
        OMPLoadWeight(1:Nthreads)=1.0
    endif
    if (OMPAutoBalance.eq.1) then
        !Check that time are not zero
        if (OMPBalanceStrength<=0) call xerrab('OMPBalanceStrength must be >0')
        if (minval(OMPTimeLocalJac).gt.0.0) then
            do i=1,Nthreads
                OMPLoadWeight(i)=OMPLoadWeight(i)*1/(OMPTimeLocalJac(i)/sum(OMPTimeLocalJac)*real(Nthreads))**OMPBalanceStrength
            enddo
        else
            OMPLoadWeight(1:Nthreads)=1.0
        endif
    endif
    !   Get the range of the iv index for each thread
    call OMPSplitIndex(1,neq,Nthreads,OMPivmin,OMPivmax,OMPLoadWeight)

    if (OMPJacVerbose.gt.0) then
        write(*,*)' *OMPJac* neq=',neq,neqmx
        write(*,*)' *OMPJac* Ivmin(ith),Ivmax(ith), OMPLoadWeight(ith) : OMPTimeLocalJac(ith) ***'
        do ith=1,Nthreads
            write(*,'(a,I3,a,I7,I7,f6.2,a,f6.2)') '  *    ithread ', ith,':',OMPivmin(ith),OMPivmax(ith),OMPLoadWeight(ith)&
            ,' : ',OMPTimeLocalJac(ith)
        enddo
    endif


    !    iJacCol(1:nnzmxperthread,1:Nthreads)=0
    !    rJacElem(1:nnzmxperthread,1:Nthreads)=0.0
    !    iJacRow(1:neq,1:Nthreads)=0
    !    if (OMPJacDebug.gt.0) then
    !        write(*,*) OMPJacStamp,' Jacobian arrays set to zero'
    !    endif
    OMPTimeJacCalc= gettime(sec4)

    !   Get initial value of system cpu timer.
    if (istimingon .eq. 1) tsjstor = gettime(sec4)

    !   Count Jacobian evaluations, both for total and for this case
    ijactot = ijactot + 1
    ijac(ig) = ijac(ig) + 1
    if (svrpkg.eq.'nksol') write(*,*) ' Updating Jacobian, npe =  ',ijac(ig)

    !   Set up diagnostic arrays for debugging
    do iv = 1, neq
        yldot_unpt(iv) = yldot00(iv)  ! for diagnostic only
        yldot_pert(iv) = 0.
    enddo

    !   build jacobian ##############################################################
    OMPTimeBuild=gettime(sec4)
    nnz(1:Nthreads)=-1
    call OMPJacBuilder(neq, t, yl,yldot00, ml, mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    OMPTimeBuild=gettime(sec4)-OMPTimeBuild
    if (istimingon .eq. 1) OMPTotTimebuild = OMPTimeBuild+OMPTotTimebuild
    if (OMPJacVerbose.gt.0) write(*,*)OMPJacStamp,' Time to build jac:',OMPTimeBuild
    !   end build jacobian ##############################################################

    !   collect jacobian ##############################################################
    OMPTimeCollect=gettime(sec4)
    call OMPCollectJacobian(neq,nnzmx,rcsc,icsc,jcsc,nnzcumout)
    OMPTimeCollect=gettime(sec4)-OMPTimeCollect
    if (istimingon .eq. 1) OMPTotTimeCollect = OMPTimeCollect+OMPTotTimeCollect
    if (OMPJacVerbose.gt.0) write(*,*)OMPJacStamp,' Time to collect jac:',OMPTimeCollect
    !   end collect jacobian ##############################################################

    jcsc(neq+1) = nnzcumout+1 ! This is set here out of OMPJAcCollect for compatibility with hybrid jac_calc

    !   for Debug purpose
    if (WriteJacobian.eq.1) then
        write(filename,'(a,3i3,a)') "jac_omp_",ijac(ig),".txt"
        call jac_write(filename,neq, rcsc, icsc, jcsc)
    endif

    if (OMPCheckNaN.gt.0) then
        do i=1,nnzmx
            if (isnan(rcsc(i))) then
                write(*,*) 'rcsc is NaN at i=',i,rcsc(i)
            endif
            if (icsc(i).ne.icsc(i)) then
                write(*,*) 'icsc is NaN at i=',i,icsc(i)
            endif
        enddo
        do i=1,neq+1
            if (jcsc(i).ne.jcsc(i)) then
                write(*,*) 'jcsc is NaN at i=',i
            endif
        enddo
    endif

    !   Convert Jacobian from compressed sparse column to compressedsparse row format.
    time1=gettime(sec4)
    call csrcsc (neq, 1, 1, rcsc, icsc, jcsc, jac, ja, ia)

    !! ... If desired, calculate Jacobian elements for ion continuity
    !!     equation by using INEL routines, and combine elements with
    !!     those calculated above.
    !      if (iondenseqn .eq. "inel") then
    !         call iondens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !      endif

    !! ... If necessary, calculate Jacobian elements for rows corresponding
    !!     to impurity-density equations (for cells other than guard cells),
    !!     and combine elements with those calculated above.
    !      if (isimpon .eq. 3 .or. isimpon .eq. 4) then
    !         if (istimingon .eq. 1) tsimpjf = gettime(sec4)
    !         call impdens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !         if (istimingon .eq. 1) then
    !            dtimpjf = gettime(sec4) - tsimpjf
    !            ttimpjf = ttimpjf + dtimpjf
    !            ttotjf = ttotjf + dtimpjf
    !         endif
    !      endif

    if (istimingon .eq. 1) ttjstor = ttjstor + gettime(sec4) - tsjstor

    OMPTimeJacCalc=gettime(sec4)-OMPTimeJacCalc
    if (istimingon .eq. 1) OMPTotJacCalc = OMPTimeJacCalc+OMPTotJacCalc
    write(*,'(a,f8.2,a,I3,a)') 'Time in omp jac_calc:',OMPTimeJacCalc,'[',Nthreads,']'
    return
end subroutine jac_calc_omp
!-------------------------------------------------------------------------------------------------
subroutine OMPJacBuilder(neq, t, yl,yldot00, ml,mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    use OMPJacSettings,only:OMPJacDebug,nnzmxperthread,OMPJacStamp,OMPJacVerbose
    use OMPSettings,only: OMPCopyArray,OMPCopyScalar,Nthreads
    use OMPJacobian, only:OMPivmin,OMPivmax,OMPTimeLocalJac,OMPTimeJacRow
    use OmpCopybbb
    use OmpCopycom
    use OmpCopyapi
    use omp_lib

    implicit none
    integer,intent(inout)::nnz(Nthreads)
    integer,intent(in):: neq      ! total number of equations (all grid points)
    integer,intent(in):: ml, mu   ! lower and upper bandwidths
    real,intent(in):: t           ! physical time
    real,intent(in) ::yl(*)       ! dependent variables
    real,intent(in) ::yldot00(neq+2) ! right-hand sides evaluated at yl
    real,intent(inout) :: wk(neq)
    integer,intent(out)::iJacCol(nnzmxperthread,nthreads)
    integer,intent(out):: iJacRow(neq,nthreads)
    real,intent(out):: rJacElem(nnzmxperthread,nthreads)
    real ::wkcopy(neq)
    real::ylcopy(neq+2)

    integer ::iJacColCopy(nnzmxperthread),iJacRowCopy(neq)
    integer ::ivmincopy(Nthreads),ivmaxcopy(Nthreads),nnzcopy(Nthreads)
    integer ::Nthreadscopy,nnzmxperthreadcopy
    real :: rJacElemCopy(nnzmxperthread),TimeJacRowcopy(neq)
    integer:: ith,tid,nnzlocal,ithcopy
    DOUBLE PRECISION :: TimeThread

    if (OMPJacDebug.gt.0)write(*,*) OMPJacStamp,' Copying data....'
    call pandf1 (-1, -1, 0, neq, 0, yl, ylcopy)
    if (OMPCopyArray.gt.0) then
        if (OMPJacDebug.gt.0)write(*,*) OMPJacStamp,' Copying array....'
        call OmpCopyPointerbbb
        call OmpCopyPointercom
        call OmpCopyPointerapi
    endif

    if (OMPCopyScalar.gt.0) then
        if (OMPJacDebug.gt.0)write(*,*) OMPJacStamp,' Copying scalar....'
        call OmpCopyScalarbbb
        call OmpCopyScalarcom
        call OmpCopyScalarapi
    endif


    !   We cannot use variables in the parallel construct declarations below when these variables are not in the scope of the subroutine
    Nthreadscopy=Nthreads
    nnzmxperthreadcopy=nnzmxperthread
    ivmincopy(1:Nthreads)=OMPivmin(1:Nthreads)
    ivmaxcopy(1:Nthreads)=OMPivmax(1:Nthreads)
    iJacColCopy(1:nnzmxperthread)=0
    rJacElemCopy(1:nnzmxperthread)=0.0
    TimeJacRowcopy(1:neq)=0
    iJacRowCopy(1:neq)=0
    ylcopy(1:neq+2)=yl(1:neq+2) ! a very barbarian use of yl(neq+1) is implemented as a switch in pandf... Error-prone!
    wkcopy(1:neq)=wk(1:neq) ! Could be set equal to zero as well. The worker wk is not an output...

    if (OMPJacDebug.gt.0) then
        write(*,*) OMPJacStamp,' Starting parallel loop'
    endif
    tid=-1
    nnzlocal=-10000

    ! ivmincopy,ivmaxcopy,yldot00, neq an t  could be shared as well as well as
    !$omp parallel do default(shared)&
    !$omp& firstprivate(ithcopy,ivmincopy,ivmaxcopy,tid,nnzlocal,ylcopy,wkcopy,ml,mu,yldot00,t,neq)&
    !$omp& firstprivate(nnzmxperthreadcopy,nthreadscopy,iJacRowCopy,iJacColCopy,rJacElemCopy,TimeJacRowcopy)&
    !$omp& private(TimeThread)

    loopthread: do ith=1,Nthreads !ith from 1 to Nthread, tid from 0 to Nthread-1
        Timethread = omp_get_wtime()
        tid=omp_get_thread_num()
        ithcopy=ith
        if (OMPJacDebug.gt.0) write(*,*) OMPJacStamp,' Thread id:',tid,' <-> ith:',ithcopy
        ! we keep all these parameters as it is easier to debug LocalJacBuilder and deal with private/shared attributes

        call LocalJacBuilder(ivmincopy(ithcopy),ivmaxcopy(ithcopy),neq, t, ylcopy,yldot00,ml,mu,wkcopy,&
            iJacColcopy,rJacElemcopy,iJacRowcopy,ithcopy,nnzlocal,nnzmxperthreadcopy,nthreadscopy,TimeJacRowcopy)
        if (OMPJacDebug.gt.0) write(*,*) OMPJacStamp,',',tid,' nzlocal:',nnzlocal

        !$omp  critical
        iJacCol(1:nnzlocal,ithcopy)=iJacColCopy(1:nnzlocal)
        rJacElem(1:nnzlocal,ithcopy)=rJacElemCopy(1:nnzlocal)
        iJacRow(1:neq,ithcopy)=iJacRowCopy(1:neq)
        OMPTimeJacRow(ivmincopy(ithcopy):ivmaxcopy(ithcopy))=TimeJacRowcopy(ivmincopy(ithcopy):ivmaxcopy(ithcopy))
        nnzcopy(ithcopy)=nnzlocal
        !$omp  end critical

        OMPTimeLocalJac(ithcopy)=omp_get_wtime() - Timethread
        if (OMPJacVerbose.gt.1) write(*,*) OMPJacStamp,' Time in thread #', tid,':',OMPTimeLocalJac(ithcopy)
        if (OMPJacVerbose.gt.1) write(*,'(a,I3,a)') 'OMP thread ',tid,' exiting...'
    enddo loopthread
    !$omp  END PARALLEL DO


    nnz(1:Nthreads)=nnzcopy(1:Nthreads) !nnzcopy is not necssary as nnz would be shared anyway in the parallel construct

    if (OMPJacDebug.gt.0) then
        write(*,*) OMPJacStamp,' End of parallel loop....'
    endif


end subroutine OMPJacBuilder
#endif
!-------------------------------------------------------------------------------------------------
    subroutine LocalJacBuilder(ivmin,ivmax,neq, t, yl,yldot00, ml, mu, wk,iJacCol,rJacElem,iJacRow,ith,nnz,nnzmxperthread,nthreads&
,TimeJacRow)

    ! ... Calculate Jacobian matrix (derivatives with respect to each
    ! ... Calculate Jacobian matrix (derivatives with respect to each
    !     dependent variable of the right-hand side of each rate equation).
    !     Lower and upper bandwidths are used to select for computation
    !     only those Jacobian elements that may be nonzero.
    !     Estimates of Jacobian elements are computed by finite differences.
    !     The Jacobian is stored in compressed sparse row format.

    use Dim, only:nx
    use Timing,only:istimingon,ttjstor,ttotjf,ttimpjf
    use Math_problem_size,only:neqmx,numvar
    use Indexes,only:igyl,iseqalg,idxphi
    use Variable_perturbation,only:delperturb,dylconst,isjacreset
    use Jacobian_clipping,only:jaccliplim,istopjac,irstop,icstop
    use Ynorm,only:suscal,sfscal
    use UEpar,only:isphion,isnewpot,svrpkg,isbcwdt
    use Model_choice,only:iondenseqn
    use Bcond,only:isextrnpf,isextrtpf,isextrngc,isextrnw,isextrtw
    use Time_dep_nwt,only:nufak,dtreal,ylodt,dtuse,dtphi
    use OMPJacSettings,only:iidebugprint,ivdebugprint,OMPTimingJacRow
    use Jacobian_csc,only:yldot_pert
    use omp_lib

    implicit none
    integer,intent(in):: ith,nnzmxperthread,nthreads,ivmin,ivmax,neq
    integer,intent(inout)::nnz
    real,intent(in):: t           ! physical time
    real,intent(inout) ::yl(*)       ! dependent variables
    integer,intent(inout)::iJacCol(nnzmxperthread),iJacRow(neq)
    real,intent(inout):: rJacElem(nnzmxperthread)
    real,intent(in) ::yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu   ! lower and upper bandwidths
    real,intent(out)::TimeJacRow(neq)
    ! ... Work-array argument:
    real,intent(inout) :: wk(neq)     ! work space available to this subroutine

    ! ... Functions:
    logical ::tstguardc
    real(kind=4) gettime
    !     real(kind=4) ranf

    ! ... Local variables:
    real ::yold, dyl, jacelem,Time
    real(kind=4) sec4, tsjstor, tsimpjf, dtimpjf
    integer:: ii, iv, jv,ii1, ii2, xc, yc, ix, iy,tid
    !       time0=gettime(sec4) #here

    ! ... Only perturb variables that are being solved for (for Daspk option)
    !      if (iseqon(iv) .eq. 0) goto 18

    ! ... Set beginning and ending indices of right-hand sides that might be
    !     perturbed.

    nnz=1
    loopiv: do iv=ivmin,ivmax

        if (OMPTimingJacRow.gt.0) Time=omp_get_wtime()
        ii1 = max(iv-mu, 1)
        ii2 = min(iv+ml, neq)
        ! ... Reset range if this is a potential perturbation with isnewpot=1
        !         if (isphion*isnewpot.eq.1 .and. mod(iv,numvar).eq.0) then
        if (isphion*isnewpot.eq.1 .and. mod(iv,numvar).eq.0) then
            ii1 = max(iv-4*numvar*nx, 1)      ! 3*nx may be excessive
            ii2 = min(iv+4*numvar*nx, neq)    ! 3*nx may be excessive
        endif
        ! ... Reset range if extrapolation boundary conditions are used
        if (isextrnpf+isextrtpf+isextrngc+isextrnw+isextrtw.gt.0) then
            ii1 = max(iv-2*numvar*(nx+3), 1)      ! guess to include extrap. bc
            ii2 = min(iv+2*numvar*(nx+3), neq)    ! guess to include extrap. bc
        endif

        ! ... Initialize all of those right-hand sides to their unperturbed
        !     values.

        do ii = ii1, ii2   ! a) below wk is reset, but only over limited range
            wk(ii) = yldot00(ii)
        enddo

        ! ... Set spatial-location indices for this dependent variable.
        xc = igyl(iv,1)
        yc = igyl(iv,2)

        ! ... Save value of dependent variable, then perturb it.
        !     The perturbation to the variable is proportional to parameter
        !     del and to a measure of the size of the variable.  That measure
        !     increases with the absolute value of the variable if it exceeds
        !     the typical size given by dylconst/suscal but can never be less
        !     than that typical size.

        yold = yl(iv)
        dyl = delperturb * (abs(yold) + dylconst / suscal(iv))
        yl(iv) = yold + dyl

        !Calculate right-hand sides near location of perturbed variable.
        !call convsr_vo (xc, yc, yl)  ! test new convsr placement
        !call convsr_aux (xc, yc, yl) ! test new convsr placement
        call pandf1 (xc, yc, iv, neq, t, yl, wk)

        !yl(iv) = yold - dyl    ! for 2nd order Jac
        !call pandf1 (xc, yc, iv, neq, t, yl, yldot0) ! for 2nd order Jac

        !Calculate possibly nonzero Jacobian elements for this variable,
        !and store nonzero elements in compressed sparse column format.

        iJacRow(iv) = nnz      ! sets index for first Jac. elem. of var. iv
                                       ! the index is set to start at 0

        loopii: do ii = ii1, ii2
            jacelem = (wk(ii) - yldot00(ii)) / dyl
            !jacelem = (wk/(ii) - yldot0(ii)) / (2*dyl)  ! for 2nd order Jac

            !Add diagonal 1/dt for nksol
            ifdiagonal:if (((svrpkg.eq."nksol") .or. (svrpkg.eq."petsc")) .and. iv.eq.ii) then
                if (iseqalg(iv)*(1-isbcwdt).eq.0) then
                    jacelem = jacelem - 1/dtuse(iv)
                endif
                ix = igyl(iv,1)
                iy = igyl(iv,2)
                if (idxphi(ix,iy)==iv .and. dtphi<1e10) then ! selects phi eqn
                    jacelem = jacelem - 1/dtphi
                endif
            endif ifdiagonal

            ! ...  Add a pseudo timestep to the diagonal #! if eqn is not algebraic
            if (svrpkg .ne. "cvode" .and. nufak .gt. 0) then
                if (iv.eq.ii .and. yl(neq+1).eq.1) jacelem = jacelem - nufak  !omit .and. iseqalg(iv).eq.0)
            !     .                   jacelem = jacelem - nufak*suscal(iv)/sfscal(iv)
            endif

            ! Debug
            debug: if (ii==iidebugprint.and.iv==ivdebugprint) then
                tid=omp_get_thread_num()
                write(*,'(a,i3,a,i7,i7,i7, E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)') '#', tid,&
                    ' : ',iv,nnz,ii,jacelem,wk(ii),yldot00(ii),sfscal(iv),jaccliplim,nufak
                call DebugHelper('dumpomp.txt') ! See at the end of this file
            endif debug

            storage: if (abs(jacelem*sfscal(iv)) .gt. jaccliplim) then
                if (nnz .gt. nnzmxperthread) then
                    write(*,*) 'nnz=',nnz,'nnzmxperthread=',nnzmxperthread
                    write(*,*) 'ith',ith,'*** jac_calc -- More storage needed for Jacobian. Increase lenpfac or omplenpfac.'
                    call xerrab("")
                endif
                !              if (rdoff.ne.0.e0) jacelem=jacelem*(1.0e0+ranf()*rdoff)
                iJacCol(nnz)=ii
                rJacElem(nnz)=jacelem
                nnz = nnz + 1
            endif storage

            check: if (istopjac.gt.0 .and. ii.eq.irstop .and. iv.eq.icstop) then
                !               yldot_pert(ii) = wk(ii)      ! for diagnostic only
                if (istopjac == 2) then
                    yl(iv) = yold
                    call pandf1 (xc, yc, iv, neq, t, yl, wk)
                endif
                !               call convsr_vo (xc, yc, yl)  ! was one call to convsr
                !               call convsr_aux (xc, yc, yl)
                call remark("***** non-zero jac_elem at irstop,icstop")
                write(*,*) 'irstop = ', irstop, ', icstop = ', icstop
                call xerrab("")
            endif check

        enddo loopii
        ! ... Restore dependent variable and plasma variables near its location.
        yl(iv) = yold
        !         call convsr_vo (xc, yc, yl)  ! was one call to convsr
        !         call convsr_aux (xc, yc, yl)

        call pandf1 (xc, yc, iv, neq, t, yl, wk)

        ! 18   continue
        !...  If this is the last variable before jumping to new cell, reset pandf

        reset: if( isjacreset.ge.1) then
            yldot_pert(1:neq)=wk(1:neq)
            ! 18   continue
            !...  If this is the last variable before jumping to new cell, reset pandf
            !JG this call to pandf1 can be safely ignored with ijacreset=0 (and save some time...)
            if (mod(iv,numvar).eq.0) then
                call pandf1 (xc, yc, iv, neq, t, yl, wk)
            endif

            do ii=1,neq
                if (yldot_pert(ii).ne.wk(ii)) then
                    write(*,'(a,i5,e20.12,e20.12)') ' *** wk modified on second call to pandf1 at ii=', ii,yldot_pert(ii),wk(ii)
                    call xerrab('*** Stop ***')
                endif
                if (isnan(yldot_pert(ii))) then
                    write(*,*) 'NaN at ii=',ii
                    call xerrab('*** Nan in wk array in jac_calc ***')
                endif

            enddo
        endif reset
        if (OMPTimingJacRow.gt.0) then
        Time=omp_get_wtime()-Time
        TimeJacRow(iv)=Time
        endif
    enddo loopiv
! ... End loop over dependent variables and finish Jacobian storage.
end subroutine LocalJacBuilder
!-------------------------------------------------------------------------------------------------
#ifdef MPIJAC
subroutine InitMPI
    Use MPIJacSettings,only:Nprocs,ComSize,MPIRank,mpilenpfac,mpineq,MPIVerbose,nnzmxperproc,MPIstamp
    Use HybridSettings,only:Hybridstamp,HybridOMPMPI
    use mpi

    use LSode,only:neq
    Use Jacobian,only:nnzmx
    implicit none
    integer(kind=4):: ierr
    integer(kind=4)::irank,icomsize
    character*8 :: MPIRankTag
    ! define rank and size

    ! check the size of the common world
    call MPI_COMM_SIZE(MPI_COMM_WORLD, iComSize, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, iRank, ierr)
    MPIrank=irank       ! irank is integer(kind=4)
    ComSize=iComSize !

    ! prepare MPI/Hybrid stamps for output
    write(MPIRankTag,'(I4)') MPIrank
    write(MPIstamp,'(a,a,a)') '[',trim(adjustl(trim(MPIRankTag))),'] MPIJac*'
    if (HybridOMPMPI>0) then
        write(Hybridstamp,'(a,a,a)') '[',trim(adjustl(trim(MPIRankTag))),'] HybridJac*'
    endif

    if (MPIVerbose>1) write(*,*) MPIStamp,' MPI enabled with ComSize=',ComSize
    !
    if (Nprocs.gt.ComSize) then
        if (MPIVerbose.gt.1) write(*,*) MPIStamp,' # processors requested > # processors available.',&
            'Nprocmax:',ComSize
        Nprocs=ComSize
        if (MPIVerbose.gt.0) write(*,*) MPIStamp,' Resetting Nproc to Nprocmax: Nproc=',Nprocs
    endif
    if (Nprocs.le.0) then
        call xerrab('Nprocs must be >0')
    endif

    if (Nprocs.gt.1) then
        nnzmxperproc=ceiling(real(nnzmx)/real(Nprocs-1)*mpilenpfac)
    else
        nnzmxperproc=nnzmx
    endif
    if (MPIVerbose.gt.0) write(*,'(a,a,i10)') MPIStamp,' nnzmxperproc:',nnzmxperproc
    mpineq=neq
    call gchange('MPIJacobian',0)
    if (MPIVerbose.gt.0) write(*,'(a,a)') MPIStamp,'  MPI variables allocated'

end subroutine InitMPI
!-------------------------------------------------------------------------------------------------
subroutine jac_calc_mpi (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)

    !   Calculate Jacobian matrix (derivatives with respect to each
    !   dependent variable of the right-hand side of each rate equation).
    !   Lower and upper bandwidths are used to select for computation
    !   only those Jacobian elements that may be nonzero.
    !   Estimates of Jacobian elements are computed by finite differences.
    !   The Jacobian is stored in compressed sparse row format.

    use Timing,only:istimingon,ttjstor,ttotjf,ttimpjf,&
    MPITotJacCalc,MPITotTimeCollect,MPITotTimeBuild
    use Grid,only:ngrid,ig,ijac,ijactot
    use Jacobian_csc,only:rcsc,jcsc,icsc,yldot_pert,yldot_unpt
    use mpi
    use MPIJacSettings,only:Nprocs,MPIVerbose,MPIDebug,nnzmxperproc,MPICheckNaN,MPIWriteJacobian,MPIRank,&
        MPILoadBalance,MPIAutoBalance,MPIStamp,MPIBalanceStrength
    use MPIJacobian,only:MPIivmin,MPIivmax,MPIiJacCol,MPIrJacElem,MPIiJacRow,MPILoadWeight,MPITimeLocalJac
    use UEpar, only: svrpkg
    implicit none
    ! ... Input arguments:
    integer,intent(in):: neq      !      total number of equations (all grid points)
    real,intent(in)   :: t              ! physical time
    real,intent(inout)   :: yl(*)          ! dependent variables
    real,intent(in)   :: yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu         ! lower and upper bandwidths
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian

    ! ... Output arguments:
    real,intent(out)   :: jac(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: ja(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: ia(neq+1)   ! pointers to beginning of each row in jac,ja

    ! ... Work-array argument:
    real wk(neq)     ! work space available to this subroutine

    ! ... Functions
    logical tstguardc
    real(kind=4) gettime

    ! ... Local variables:
    real factor
    real(kind=4) sec4, tsjstor, tsimpjf, dtimpjf,time0,time1,MPITimeBuild,MPITimeCollect,MPITimeJacCalc
    integer:: i,iv,iproc,nnz
    integer (kind=4) :: ierr
    character(len = 80) ::  filename
    !integer::iJacCol(1:nnzmxperproc)
    !real ::rJacElem(1:nnzmxperproc)
    if (MPILoadBalance.ne.1 .and. MPIAutoBalance.ne.1) then
        MPILoadWeight(0:Nprocs-1)=1.0
    endif
    if (MPIAutoBalance.eq.1) then
    if (MPIBalanceStrength<=0) call xerrab('MPIBalanceStrength must be >0')
        !Check that time are not zero
        if (minval(MPITimeLocalJac(0:Nprocs-1)).gt.0.0) then
            do i=0,Nprocs-1
                factor=1/(MPITimeLocalJac(i)/sum(MPITimeLocalJac(0:Nprocs-1))*real(Nprocs))**MPIBalanceStrength
                MPILoadWeight(i)=MPILoadWeight(i)*factor
            enddo
        else
            MPILoadWeight(0:Nprocs-1)=1.0
        endif
    endif
    if (MPIDebug.gt.0) write(*,*) MPIStamp,' Starting jac_calc_mpi'
    !   Get the range of the iv index for each thread
    call MPISplitIndex(neq,Nprocs,MPIivmin,MPIivmax,MPILoadWeight)

    if (MPIVerbose.gt.1) then
        write(*,*)MPIStamp,'neq=',neq
        write(*,*)MPIStamp,' MPIivmin | MPIivmax | MPILoadWeight | MPITimeLocalJac'
        do iproc=0,Nprocs-1
    write(*,'(a7,I3,a,I7,I7,f5.1,f8.3)') 'rank:', iproc,':',MPIivmin(iproc),MPIivmax(iproc),&
    MPILoadWeight(iproc),MPITimeLocalJac(iproc)
        enddo
    endif
    MPITimeJacCalc= gettime(sec4)

    !   Get initial value of system cpu timer.
    if (istimingon .eq. 1) tsjstor = gettime(sec4)

    !   Count Jacobian evaluations, both for total and for this case
    ijactot = ijactot + 1
    ijac(ig) = ijac(ig) + 1
    if (svrpkg.eq.'nksol') write(*,*) ' Updating Jacobian, npe =  ',ijac(ig)

    !   Set up diagnostic arrays for debugging
    do iv = 1, neq
        yldot_unpt(iv) = yldot00(iv)  ! for diagnostic only
        yldot_pert(iv) = 0.
    enddo

    !   build jacobian ##############################################################
    MPITimeBuild=gettime(sec4)

    call MPIJacBuilder(neq, t, yl,yldot00, ml, mu,wk,MPIiJacCol,MPIrJacElem,MPIiJacRow,nnz)

    MPITimeBuild=gettime(sec4)-MPITimeBuild
    if (istimingon .eq. 1) MPITotTimebuild = MPITimeBuild+MPITotTimebuild
    if (MPIVerbose.gt.0) write(*,*)MPIStamp,' Rank:',MPIRank ,'Time to build jac:',MPITimeBuild

    call MPICollectBroadCastTime(real(MPITimeBuild,kind=8))
    !   end build jacobian ##############################################################
    !   collect jacobian ##############################################################

    MPITimeCollect=gettime(sec4)
    call MPICollectBroadCastJacobian(MPIiJacRow,MPIiJacCol,MPIrJacElem,nnz)
    MPITimeCollect=gettime(sec4)-MPITimeCollect
    if (istimingon .eq. 1) MPITotTimeCollect = MPITimeCollect+MPITotTimeCollect
    if (MPIVerbose.gt.0) write(*,*)MPIStamp,' Time to collect/broadcast jac:',MPITimeCollect
    !   end collect jacobian ##############################################################


    !   for Debug purpose
    if (MPIWriteJacobian.eq.1) then
        write(filename,'(a,3i3,a)') "jac_MPI_",ijac(ig),".txt"
        call jac_write(filename,neq, rcsc, icsc, jcsc)
    endif

    if (MPICheckNaN.gt.0) then
        do i=1,nnzmx
            if (isnan(rcsc(i))) then
                write(*,*) 'rcsc is NaN at i=',i,rcsc(i)
            endif
            if (icsc(i).ne.icsc(i)) then
                write(*,*) 'icsc is NaN at i=',i,icsc(i)
            endif
        enddo
        do i=1,neq+1
            if (jcsc(i).ne.jcsc(i)) then
                write(*,*) 'jcsc is NaN at i=',i
            endif
        enddo
    endif

    !   Convert Jacobian from cMPIressed sparse column to compressedsparse row format.
    time1=gettime(sec4)
    call csrcsc (neq, 1, 1, rcsc, icsc, jcsc, jac, ja, ia)

    !! ... If desired, calculate Jacobian elements for ion continuity
    !!     equation by using INEL routines, and combine elements with
    !!     those calculated above.
    !      if (iondenseqn .eq. "inel") then
    !         call iondens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !      endif

    !! ... If necessary, calculate Jacobian elements for rows corresponding
    !!     to impurity-density equations (for cells other than guard cells),
    !!     and combine elements with those calculated above.
    !      if (isimpon .eq. 3 .or. isimpon .eq. 4) then
    !         if (istimingon .eq. 1) tsimpjf = gettime(sec4)
    !         call impdens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !         if (istimingon .eq. 1) then
    !            dtimpjf = gettime(sec4) - tsimpjf
    !            ttimpjf = ttimpjf + dtimpjf
    !            ttotjf = ttotjf + dtimpjf
    !         endif
    !      endif

    if (istimingon .eq. 1) ttjstor = ttjstor + gettime(sec4) - tsjstor
    MPITimeJacCalc=gettime(sec4)-MPITimeJacCalc
    if (istimingon .eq. 1) MPITotJacCalc = MPITimeJacCalc+MPITotJacCalc
     write(*,'(a,f8.2)') 'Time in jac_calc:',MPITimeJacCalc
    call MPI_barrier(MPI_COMM_WORLD,ierr)
    return
end subroutine jac_calc_mpi
!-------------------------------------------------------------------------------------------------
subroutine MPICollectBroadCastTime(TimeLocal)

    use mpi
    use MPIJacSettings,only:Nprocs,MpiRank,MPIVerbose,MPIDebug,ioutmpi,MPIStamp
    use MPIJacobian,only: MPITimeLocalJac

    implicit none
    integer(kind=4)::iproc
    integer(kind=4) :: ierr,req0
    real,intent(in) :: TimeLocal
    real(kind=4) gettime
    real(kind=4) sec4, TimeCollect
    TimeCollect=gettime(sec4)
    if (MPIRank==0) MPITimeLocalJac(0)=Timelocal
    loopproc:do iproc=1,Nprocs-1
        if (MPIRank.eq.iproc) then
            ! send to the master proc
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank',MPIRAnk,' sending data to 0'
            endif
            CALL MPI_SEND(TimeLocal,1,MPI_REAL8,0,7,MPI_COMM_WORLD,ierr)
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank',iproc,'has sent data to 0'
            endif
            ! collect on the master proc
        elseif (MPIRank.eq.0) then
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank 0 getting data from:',iproc
            endif
            CALL MPI_RECV(Timelocal,1,MPI_REAL8,iproc,7,MPI_COMM_WORLD,req0,ierr)
            MPITimeLocalJac(iproc)=Timelocal
        endif
    enddo loopproc
    call MPI_barrier(MPI_COMM_WORLD,ierr)
    TimeCollect=gettime(sec4)-TimeCollect
    if (MPIVerbose.gt.0) write(*,*)MPIStamp,' Time to collect time:',TimeCollect
    TimeCollect=gettime(sec4)
    ! And now broadcast it to all the procs
    call MPI_BCAST(MPITimeLocalJac(0:Nprocs-1), int(Nprocs,kind=4), MPI_Real8, 0, MPI_COMM_WORLD, IERR)
    call MPI_barrier(MPI_COMM_WORLD,ierr)

    TimeCollect=gettime(sec4)-TimeCollect
    if (MPIVerbose.gt.0) write(*,*)MPIStamp,' Time to broadcast time:',TimeCollect
end subroutine MPICollectBroadCastTime
!-----------------------------------------------------------------------
subroutine MPICollectBroadCastJacobian(iJacRow,iJacCol,rJacElem,nnz)

    use mpi
    use MPIJacSettings,only:Nprocs,MpiRank,mpineq,nnzmxperproc,MPIVerbose,MPIDebug,ioutmpi,MPIStamp
    use MPIJacobian,only: MPIivmin,MPIivmax
    use Jacobian_csc,only:rcsc,jcsc,icsc
    use LSode,only:neq
    Use Jacobian,only:nnzmx

    implicit none
    integer:: nnzcumlocal(Nprocs),ith
    integer(kind=4)::iproc,nnzmxperproclocal,mpineqlocal
    integer ,intent(in):: nnz
    integer,intent(in)::iJacCol(nnzmxperproc)
    integer,intent(in):: iJacRow(mpineq+1)
    real,intent(in):: rJacElem(nnzmxperproc)
    integer(kind=4) :: ierr,req0,req1,req2,req3
    real :: TimeCollect
    real(kind=4) gettime
    real(kind=4) sec4, TimeJacCalc
    TimeCollect=gettime(sec4)
    nnzmxperproclocal=int(nnzmxperproc,kind=4)
    mpineqlocal=int(mpineq,kind=4)
    ! First we collect data from the master proc (Mpirank=0)
    if (MPIDebug.gt.0) write(ioutmpi,*) 'Rank',MPIrank,'nnz-1:',nnz-1
    if (MPIRank.eq.0) then
        nnzcumlocal(1)=nnz-1
        jcsc(MPIivmin(0):MPIivmax(0))= iJacRow(MPIivmin(0):MPIivmax(0))
        rcsc(1:nnz-1)= rJacElem(1:nnz-1)
        icsc(1:nnz-1)= iJacCol(1:nnz-1)
    endif
    ! Then we collect from each non-master process (MPIrank>0)

    loopproc:do iproc=1,Nprocs-1
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        ith=iproc+1
        if (MPIRank.eq.iproc) then
            ! send to the master proc
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank',MPIRAnk,' sending data to 0'
            endif
            CALL MPI_SEND(nnz,1,MPI_INTEGER8,0,7,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(iJacRow,mpineqlocal,MPI_INTEGER8,0,9,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(iJacCol,nnzmxperproclocal,MPI_INTEGER8,0,10,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(rJacElem,nnzmxperproclocal,MPI_REAL8,0,11,MPI_COMM_WORLD,ierr)
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank',iproc,'has sent data to 0'
            endif
            ! collect on the master proc
        elseif (MPIRank.eq.0) then
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank 0 getting data from:',iproc
            endif
            CALL MPI_RECV(nnz,1,MPI_INTEGER8,iproc,7,MPI_COMM_WORLD,req0,ierr)
            CALL MPI_RECV(iJacRow,mpineqlocal,MPI_INTEGER8,iproc,9,MPI_COMM_WORLD,req1,ierr)
            CALL MPI_RECV(iJacCol,nnzmxperproclocal,MPI_INTEGER8,iproc,10,MPI_COMM_WORLD,req2,ierr)
            CALL MPI_RECV(rJacElem,nnzmxperproclocal,MPI_REAL8,iproc,11,MPI_COMM_WORLD,req3,ierr)
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank 0 got data from:',iproc, '; ivmin:ivmax=',MPIivmin(iproc),MPIivmax(iproc)
            endif
            CALL MPI_WAIT(req0, MPI_STATUS_IGNORE, ierr)
            CALL MPI_WAIT(req1, MPI_STATUS_IGNORE, ierr)
            CALL MPI_WAIT(req2, MPI_STATUS_IGNORE, ierr)
            CALL MPI_WAIT(req3, MPI_STATUS_IGNORE, ierr)
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank 0 going next'
            endif

            ! Now build the jacobian vector iteratively
            nnzcumlocal(ith)=nnzcumlocal(ith-1)+nnz-1
            if (MPIDebug.gt.0) then
                write(ioutmpi,*) MPIStamp,' Rank:',iproc,' sent: nnz:',nnz, ';nnzcum:',nnzcumlocal(ith)
            endif
            if (nnzcumlocal(ith).gt.nnzmx) then
                write(*,*) 'nnzcum=',nnzcumlocal(ith),'nnzmx=',nnzmx
                call xerrab(' Problem: nnzcum > nnzmx... Increase lenpfac')
            endif
            jcsc(MPIivmin(iproc):MPIivmax(iproc))= iJacRow(MPIivmin(iproc):MPIivmax(iproc))+nnzcumlocal(ith-1)
            rcsc(nnzcumlocal(ith-1)+1:nnzcumlocal(ith))=rJacElem(1:nnz-1)
            icsc(nnzcumlocal(ith-1)+1:nnzcumlocal(ith))=iJacCol(1:nnz-1)

        endif
        if (MPIDebug.gt.0) then
            write(ioutmpi,*) MPIStamp,' Rank:',MPIRank,' at end of loop'
        endif
        call MPI_barrier(MPI_COMM_WORLD,ierr)
    enddo loopproc

    if (MPIRank.eq.0) then
        jcsc(neq+1) = nnzcumlocal(Nprocs)+1
        if (MPIVerbose.gt.0) write(*,'(a,i9)') '**** Number of non-zero Jacobian elems:',nnzcumlocal(Nprocs)
    endif

    call MPI_barrier(MPI_COMM_WORLD,ierr)
    TimeCollect=gettime(sec4)-TimeCollect
    if (MPIVerbose.gt.0) write(*,*)MPIStamp,' Time to collect jac:',TimeCollect
    TimeCollect=gettime(sec4)
    ! And now broadcast it to all the procs
    call MPI_BCAST(rcsc, nnzmx, MPI_REAL8, 0, MPI_COMM_WORLD, IERR)
    call MPI_BCAST(jcsc, mpineq+1, MPI_INTEGER8, 0, MPI_COMM_WORLD, IERR)
    call MPI_BCAST(icsc, nnzmx, MPI_INTEGER8, 0, MPI_COMM_WORLD, IERR)

    call MPI_barrier(MPI_COMM_WORLD,ierr)
    TimeCollect=gettime(sec4)-TimeCollect
    if (MPIVerbose.gt.0) write(*,*)MPIStamp,' Time to broadcast jac:',TimeCollect
end subroutine MPICollectBroadCastJacobian
!-----------------------------------------------------------------------
subroutine MPIJacBuilder(neq, t, yl,yldot00, ml,mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    use MPIJacSettings,only:MPIDebug,Nprocs,nnzmxperproc,MPIRank,ioutmpi,MPIStamp
    use MPIJacobian, only:MPIivmin,MPIivmax

    Use Jacobian,only:nnzmx


    implicit none
    integer,intent(inout)::nnz
    integer,intent(in):: neq      ! total number of equations (all grid points)
    integer,intent(in):: ml, mu   ! lower and upper bandwidths
    real,intent(in):: t           ! physical time
    real,intent(inout) ::yl(*)       ! dependent variables
    real,intent(in) ::yldot00(neq+2) ! right-hand sides evaluated at yl
    real,intent(inout) :: wk(neq)
    integer,intent(out)::iJacCol(nnzmxperproc)
    integer,intent(out):: iJacRow(neq)
    real,intent(out):: rJacElem(nnzmxperproc)
    real:: TimeJacRow(neq)
    integer :: iproc

    if (MPIDebug.gt.0) write(ioutmpi,*) MPIStamp,' Building jacobian on proc:',MPIrank

    call LocalJacBuilder(MPIivmin(MPIrank),MPIivmax(MPIrank),neq, t, yl,yldot00,ml,mu,wk,&
        iJacCol,rJacElem,iJacRow,MPIrank,nnz,nnzmxperproc,Nprocs,TimeJacRow)


end subroutine MPIJacBuilder
!-----------------------------------------------------------------------
subroutine jac_calc_hybrid (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)

    !   Calculate Jacobian matrix (derivatives with respect to each
    !   dependent variable of the right-hand side of each rate equation).
    !   Lower and upper bandwidths are used to select for computation
    !   only those Jacobian elements that may be nonzero.
    !   Estimates of Jacobian elements are computed by finite differences.
    !   The Jacobian is stored in compressed sparse row format.

    use Timing,only:istimingon,ttjstor,ttotjf,ttimpjf
    use OMPTiming,only:OMPTotTimeCollect,OMPTotTimeBuild,OMPTotJacCalc
    use MPITiming,only:MPITotTimeBuild,MPITotJacCalc,MPITotTimeCollect
    use Grid,only:ngrid,ig,ijac,ijactot
    use Jacobian_csc,only:rcsc,jcsc,icsc,yldot_pert,yldot_unpt
    use mpi
    use MPIJacSettings,only:Nprocs,nnzmxperproc,MPICheckNaN,MPIWriteJacobian,MPIRank,&
        MPILoadBalance,MPIAutoBalance,MPIBalanceStrength
    use OMPSettings,only:Nthreads
    use OMPJacSettings,only:OMPLoadBalance,OMPAutoBalance,OMPBalanceStrength
    use MPIJacobian,only:MPIivmin,MPIivmax,MPIiJacCol,MPIrJacElem,MPIiJacRow,MPILoadWeight,MPITimeLocalJac
    use OMPJacobian,only:iJacCol,rJacElem,iJacRow,OMPivmin,OMPivmax,nnz,OMPLoadWeight,OMPTimeLocalJac
    use HybridSettings,only:HybridDebug,HybridVerbose,HybridCheckNaN,Hybridstamp
    use UEpar, only: svrpkg
    implicit none
    ! ... Input arguments:
    integer,intent(in):: neq      !      total number of equations (all grid points)
    real,intent(in)   :: t              ! physical time
    real,intent(in)   :: yl(*)          ! dependent variables
    real,intent(in)   :: yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu         ! lower and upper bandwidths
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian

    ! ... Output arguments:
    real,intent(out)   :: jac(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: ja(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: ia(neq+1)   ! pointers to beginning of each row in jac,ja

    ! ... Work-array argument:
    real wk(neq)     ! work space available to this subroutine

    ! ... Functions
    logical tstguardc
    real(kind=4) gettime

    ! ... Local variables:
    real factor
    real(kind=4)        :: sec4, tsjstor, tsimpjf, dtimpjf,time0,time1
    real(kind=4)        :: OMPTimeBuild,MPITimeJacCalc,OMPTimeCollect,MPITimeBuild,MPITimeCollect
    integer             :: i,iv,iproc,ithread,nnzcumout
    integer(kind=4)     :: ierr
    character(len = 80) ::  filename

    if (MPILoadBalance.ne.1 .and. MPIAutoBalance.ne.1) then
        MPILoadWeight(0:Nprocs-1)=1.0
    endif
    if (MPIAutoBalance.eq.1) then
        if (MPIBalanceStrength<=0) call xerrab('MPIBalanceStrength must be >0')
        !Check that time are not zero
        if (minval(MPITimeLocalJac(0:Nprocs-1)).gt.0.0) then
            do i=0,Nprocs-1
                factor=1/(MPITimeLocalJac(i)/sum(MPITimeLocalJac(0:Nprocs-1))*real(Nprocs))**MPIBalanceStrength
                MPILoadWeight(i)=MPILoadWeight(i)*factor
            enddo
        else
            MPILoadWeight(0:Nprocs-1)=1.0
        endif
    endif

    if (OMPLoadBalance.ne.1 .and. OMPAutoBalance.ne.1) then
        OMPLoadWeight(1:Nthreads)=1.0
    endif
    if (OMPAutoBalance.eq.1) then
    if (OMPBalanceStrength<=0) call xerrab('OMPBalanceStrength must be >0')
        !Check that time are not zero
        if (minval(OMPTimeLocalJac).gt.0.0) then
            do i=1,Nthreads
                factor=1/(OMPTimeLocalJac(i)/sum(OMPTimeLocalJac)*real(Nthreads))**OMPBalanceStrength
                OMPLoadWeight(i)=OMPLoadWeight(i)*factor
            enddo
        else
            OMPLoadWeight(1:Nthreads)=1.0
        endif
    endif


    if (HybridDebug.gt.0) write(*,*) Hybridstamp,' Starting jac_calc_hybrid'
    !   Get the range of the iv index for each thread
    call MPISplitIndex(neq,Nprocs,MPIivmin,MPIivmax,MPILoadWeight)

    call OMPSplitIndex(MPIivmin(MPIRank),MPIivmax(MPIRank),Nthreads,OMPivmin,OMPivmax,OMPLoadWeight)

    if (HybridVerbose.gt.1) then
        write(*,*)HybridStamp, 'neq=',neq
        write(*,*)HybridStamp,'| MPIivmin MPIivmax MPILoadWeight MPITimeJac | OMPivmin OMPivmax OMPLoadWeight '
        do iproc=0,Nprocs-1
            do ithread=1,Nthreads
                if (ithread==1) then
                    if (iproc==MPIRank) then
                    write(*,'(a6,I3,a7,I3,a3,I10,I10,f10.1,f10.3,a3,I8,I8,f8.1)') 'rank', iproc,'thread', ithread,'|',&
                        MPIivmin(iproc),MPIivmax(iproc),MPILoadWeight(iproc),MPITimeLocalJac(iproc),&
                        '| ',OMPivmin(ithread),OMPivmax(ithread),OMPLoadWeight(ithread)
                    else
                    write(*,'(a6,I3,a7,I3,a3,I10,I10,f10.1,f10.3,a3,a8,a8,a8)') 'rank', iproc,'thread', ithread,'|',&
                        MPIivmin(iproc),MPIivmax(iproc),MPILoadWeight(iproc),MPITimeLocalJac(iproc),&
                        '| ','-','-','-'
                        endif
                else
                if (HybridVerbose.gt.3) then
                    write(*,'(a6,a3,a7,I3,a3,a10,a10,a10,a10,a3,I8,I8,f8.1)') ' ', ' ',  'thread', ithread,'|',&
                        ' ',' ',' ',' ',' | ',OMPivmin(ithread),OMPivmax(ithread),OMPLoadWeight(ithread)
                endif
                endif
            enddo
        enddo
    endif


    MPITimeJacCalc= gettime(sec4)

    !   Get initial value of system cpu timer.
    if (istimingon .eq. 1) tsjstor = gettime(sec4)

    !   Count Jacobian evaluations, both for total and for this case
    ijactot = ijactot + 1
    ijac(ig) = ijac(ig) + 1
    if (svrpkg.eq.'nksol') write(*,*) ' Updating Jacobian, npe =  ',ijac(ig)

    !   Set up diagnostic arrays for debugging
    do iv = 1, neq
        yldot_unpt(iv) = yldot00(iv)  ! for diagnostic only
        yldot_pert(iv) = 0.
    enddo

    !   build jacobian #####################################################################
    OMPTimeBuild=gettime(sec4)
    call OMPJacBuilder(neq, t, yl,yldot00, ml,mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    OMPTimeBuild=gettime(sec4)-OMPTimeBuild
    if (istimingon .eq. 1) OMPTotTimebuild = OMPTimeBuild+OMPTotTimebuild
    if (HybridVerbose.gt.2) write(*,*)Hybridstamp,' Time to build OMP jac:',OMPTimeBuild

    ! Collect OMP part of the Jacobian ######################################################
    OMPTimeCollect=gettime(sec4)
    call OMPCollectJacobian(neq,nnzmxperproc,MPIrJacElem,MPIiJacCol,MPIiJacRow,nnzcumout)
    OMPTimeCollect=gettime(sec4)-OMPTimeCollect
    if (istimingon .eq. 1) OMPTotTimeCollect = OMPTimeCollect+OMPTotTimeCollect
    if (HybridVerbose.gt.2) write(*,*)Hybridstamp,' Time to collect OMP jac:',OMPTimeCollect

    ! broadcast  MPItimebuild for load balancing #############################################
    MPITimeBuild=OMPTimeCollect+OMPTimeBuild
    if (istimingon .eq. 1) MPITotTimebuild = MPITimeBuild+MPITotTimebuild
    if (HybridVerbose.gt.2) write(*,*)Hybridstamp,' Time to build MPI jac:',MPITimeBuild
    call MPICollectBroadCastTime(real(MPITimeBuild,kind=8))


    !   collect MPI jacobian    ##############################################################
    MPITimeCollect=gettime(sec4)
    !nnzcumout is the total amount of non-zero Jac elements stored in iJacCol and rJacElem
    ! but nnz=nnzcumout+1 is passed in the collector
    call MPICollectBroadCastJacobian(MPIiJacRow,MPIiJacCol,MPIrJacElem,nnzcumout+1)
    MPITimeCollect=gettime(sec4)-MPITimeCollect
    if (istimingon .eq. 1) MPITotTimeCollect = MPITimeBuild+MPITotTimeCollect
    if (HybridVerbose.gt.2) write(*,*)Hybridstamp,' MPI: Time to collect/broadcast jac:',MPITimeCollect


    !   for Debug purpose
    if (MPIWriteJacobian.eq.1) then
        write(filename,'(a,3i3,a)') "jac_hybrid_",ijac(ig),".txt"
        call jac_write(filename,neq, rcsc, icsc, jcsc)
    endif

    if (HybridCheckNaN.gt.0) then
        do i=1,nnzmx
            if (isnan(rcsc(i))) then
                write(*,*) 'rcsc is NaN at i=',i,rcsc(i)
            endif
            if (icsc(i).ne.icsc(i)) then
                write(*,*) 'icsc is NaN at i=',i,icsc(i)
            endif
        enddo
        do i=1,neq+1
            if (jcsc(i).ne.jcsc(i)) then
                write(*,*) 'jcsc is NaN at i=',i
            endif
        enddo
    endif

    !   Convert Jacobian from cMPIressed sparse column to compressedsparse row format.
    time1=gettime(sec4)
    call csrcsc (neq, 1, 1, rcsc, icsc, jcsc, jac, ja, ia)

    !! ... If desired, calculate Jacobian elements for ion continuity
    !!     equation by using INEL routines, and combine elements with
    !!     those calculated above.
    !      if (iondenseqn .eq. "inel") then
    !         call iondens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !      endif

    !! ... If necessary, calculate Jacobian elements for rows corresponding
    !!     to impurity-density equations (for cells other than guard cells),
    !!     and combine elements with those calculated above.
    !      if (isimpon .eq. 3 .or. isimpon .eq. 4) then
    !         if (istimingon .eq. 1) tsimpjf = gettime(sec4)
    !         call impdens2 (neq, nnzmx, jac, ja, ia, rcsc, icsc, jcsc)
    !         if (istimingon .eq. 1) then
    !            dtimpjf = gettime(sec4) - tsimpjf
    !            ttimpjf = ttimpjf + dtimpjf
    !            ttotjf = ttotjf + dtimpjf
    !         endif
    !      endif

    if (istimingon .eq. 1) ttjstor = ttjstor + gettime(sec4) - tsjstor
    MPITimeJacCalc=gettime(sec4)-MPITimeJacCalc
    if (istimingon .eq. 1) MPITotJacCalc = MPITimeJacCalc+MPITotJacCalc
    write(*,*)HybridStamp,'Time in hybrid jac_calc:',MPITimeJacCalc

    return
end subroutine jac_calc_hybrid
#endif
!-------------------------------------------------------------------------------------------------
subroutine MPISplitIndex(neq,Nprocs,ivmin,ivmax,weight)
    implicit none
    integer,intent(in) ::neq,Nprocs
    real,intent(inout) :: weight(0:Nprocs-1)
    integer,intent(out)::ivmin(0:Nprocs-1),ivmax(0:Nprocs-1)
    integer:: Nsize(0:Nprocs-1),imax,i

    if (Nprocs.gt.1) then
        do i=0,Nprocs-1
                if (weight(i)<=0) call xerrab('MPISplitIndex: weight <0')
        enddo
        ! Normalized weights
        weight(0:Nprocs-1)=weight(0:Nprocs-1)/sum(weight(0:Nprocs-1))*real(Nprocs)
        do i=0,Nprocs-1
            Nsize(i)=int(real(neq/Nprocs)*weight(i))
        enddo

        do i=0,Nprocs-1
            if (Nsize(i)<0) call xerrab('Nsize<0')
            if (Nsize(i)<2) Nsize(i)=Nsize(i)+1
        enddo

        if (neq.ne.sum(Nsize(0:Nprocs-1))) then
            imax=0
            do i=1,Nprocs-1
                if (Nsize(i)>Nsize(i-1)) then
                    imax=i
                endif
            enddo
            Nsize(imax) = Nsize(imax) + (neq-sum(Nsize(0:Nprocs-1)))
        endif
        if (neq.ne.sum(Nsize(0:Nprocs-1))) call xerrab('Nsize .ne. neq!!!')

        ivmin(0)=1
        ivmax(0)=1+Nsize(0)-1
        do i=1,Nprocs-1
            ivmin(i)=ivmax(i-1)+1
            ivmax(i)=ivmin(i)+Nsize(i)-1
        enddo
        if (ivmax(Nprocs-1)-ivmin(0)+1.ne.neq) call xerrab('ivmax(Nprocs-1)!=neq')
    else
        ivmin(0)=1
        ivmax(0)=neq
    endif

end subroutine MPISplitIndex
!-------------------------------------------------------------------------------------------------
subroutine OMPSplitIndex(ieqmin,ieqmax,Nthreads,ivmin,ivmax,weight)
    implicit none
    integer,intent(in) ::ieqmin,ieqmax,Nthreads
    real::weight(Nthreads)
    integer,intent(out)::ivmin(Nthreads),ivmax(Nthreads)
    integer:: Nsize(Nthreads),Msize,R,i,imax
    if (ieqmax-ieqmin+1<2) call xerrab('Number of equations to solve <2')
    if (Nthreads.gt.1) then
        !        if (OMPLoadWeight.eq.1) then

        do i=1,Nthreads
            if (weight(i)<=0) call xerrab('OMPSplitIndex: weight <0')
            !write(*,*) weight(i)
        enddo

        ! Normalized weights
        weight(1:Nthreads)=weight(1:Nthreads)/sum(weight(1:Nthreads))*real(Nthreads)
        do i=1,Nthreads
            Nsize(i)=int(real((ieqmax-ieqmin+1)/Nthreads)*weight(i))
            !write(*,*) Nsize(i),weight(i)
        enddo

        do i=1,Nthreads
            if (Nsize(i)<0) call xerrab('Nsize<0')
            if (Nsize(i)<2) Nsize(i)=Nsize(i)+1
        enddo
        if (ieqmax-ieqmin+1.ne.sum(Nsize)) then
            imax=1
            do i=2,Nthreads
                if (Nsize(i)>Nsize(i-1)) then
                    imax=i
                endif
            enddo
            Nsize(imax) = Nsize(imax) + ((ieqmax-ieqmin+1)-sum(Nsize))
        endif
        !write(*,*) Nsize,neq,sum(Nsize)
        if (ieqmax-ieqmin+1.ne.sum(Nsize)) call xerrab('Nsize .ne. neq!!!')
        ivmin(1)=ieqmin
        ivmax(1)=ieqmin+Nsize(1)-1
        do i=2,Nthreads
            ivmin(i)=ivmax(i-1)+1
            ivmax(i)=ivmin(i)+Nsize(i)-1
        enddo
        if (ivmax(Nthreads)-ivmin(1)+1.ne.(ieqmax-ieqmin+1)) call xerrab('ivmax(Nthreads)!=neq')
    else
        ivmin(Nthreads)=ieqmin
        ivmax(Nthreads)=ieqmax
    endif

end subroutine OMPSplitIndex
!-----------------------------------------------------------------------
! The routines below are for debugging of OMP implementation
subroutine WriteArrayReal(array,s,iu)
    implicit none
    real:: array(*)
    integer:: i,s,iu
    do i=1,s
        write(iu,*) array(i)
    enddo
end subroutine WriteArrayReal

subroutine WriteArrayInteger(array,s,iu)
    implicit none
    integer:: array(*)
    integer:: i,s,iu
    do i=1,s
        write(iu,*) array(i)
    enddo
end subroutine WriteArrayInteger
!JG 03/2020 Subroutine for debugging purpose: print out all the threadprivate variable.
!JG This routine was created with the OMPDebugger.py script, which also contains methods to analyze the outpout of this routine.
!JG This routine can be called within the jacobian constructor to determine which variables behave differenlty in serial vs openmp jacobian calculations
      subroutine DebugHelper(FileName)

      Use Bcond
      Use Cfric
      Use Coefeq
      Use Comflo
      Use Comgeo
      Use Compla
      Use Comtra
      Use Conduc
      Use Dim
      Use Gradients
      Use Imprad
      Use Imslwrk
      Use Indices_domain_dcl
      Use Jacobian_csc
      Use Locflux
      Use MCN_sources
      Use PNC_params
      Use Poten
      Use Rccoef
      Use Reduced_ion_interface
      Use Rhsides
      Use Save_terms
      Use Selec
      Use Time_dep_nwt
      Use Timespl
      Use Timing
      Use UEpar
      Use Wkspace
      implicit none
      integer:: iunit
      character(len = *) ::  filename
      open (newunit = iunit, file = trim(filename))
      write(iunit,*) "alfe"
      call WriteArrayReal(alfe,size(alfe),iunit)
      write(iunit,*) "alfneo"
      call WriteArrayReal(alfneo,size(alfneo),iunit)
      write(iunit,*) "bcel"
      call WriteArrayReal(bcel,size(bcel),iunit)
      write(iunit,*) "bcer"
      call WriteArrayReal(bcer,size(bcer),iunit)
      write(iunit,*) "bcil"
      call WriteArrayReal(bcil,size(bcil),iunit)
      write(iunit,*) "bcir"
      call WriteArrayReal(bcir,size(bcir),iunit)
      write(iunit,*) "betai"
      call WriteArrayReal(betai,size(betai),iunit)
      write(iunit,*) "betap"
      call WriteArrayReal(betap,size(betap),iunit)
      write(iunit,*) "cfneut"
      write(iunit,*) cfneut
      write(iunit,*) "cfneutdiv"
      write(iunit,*) cfneutdiv
      write(iunit,*) "cfvcsx"
      call WriteArrayReal(cfvcsx,size(cfvcsx),iunit)
      write(iunit,*) "cfvcsy"
      call WriteArrayReal(cfvcsy,size(cfvcsy),iunit)
      write(iunit,*) "cfvgpx"
      call WriteArrayReal(cfvgpx,size(cfvgpx),iunit)
      write(iunit,*) "cfvgpy"
      call WriteArrayReal(cfvgpy,size(cfvgpy),iunit)
      write(iunit,*) "cmneut"
      write(iunit,*) cmneut
      write(iunit,*) "cmneutdiv"
      write(iunit,*) cmneutdiv
      write(iunit,*) "coef1"
      write(iunit,*) coef1
      write(iunit,*) "coll_fe"
      call WriteArrayReal(coll_fe,size(coll_fe),iunit)
      write(iunit,*) "coll_fi"
      call WriteArrayReal(coll_fi,size(coll_fi),iunit)
      write(iunit,*) "conx"
      call WriteArrayReal(conx,size(conx),iunit)
      write(iunit,*) "conxe"
      call WriteArrayReal(conxe,size(conxe),iunit)
      write(iunit,*) "conxg"
      call WriteArrayReal(conxg,size(conxg),iunit)
      write(iunit,*) "conxge"
      call WriteArrayReal(conxge,size(conxge),iunit)
      write(iunit,*) "conxi"
      call WriteArrayReal(conxi,size(conxi),iunit)
      write(iunit,*) "cony"
      call WriteArrayReal(cony,size(cony),iunit)
      write(iunit,*) "conye"
      call WriteArrayReal(conye,size(conye),iunit)
      write(iunit,*) "conyg"
      call WriteArrayReal(conyg,size(conyg),iunit)
      write(iunit,*) "conyge"
      call WriteArrayReal(conyge,size(conyge),iunit)
      write(iunit,*) "conyi"
      call WriteArrayReal(conyi,size(conyi),iunit)
      write(iunit,*) "cs"
      write(iunit,*) cs
      write(iunit,*) "ctaue"
      call WriteArrayReal(ctaue,size(ctaue),iunit)
      write(iunit,*) "ctaui"
      call WriteArrayReal(ctaui,size(ctaui),iunit)
      write(iunit,*) "dclass_e"
      call WriteArrayReal(dclass_e,size(dclass_e),iunit)
      write(iunit,*) "dclass_i"
      call WriteArrayReal(dclass_i,size(dclass_i),iunit)
      write(iunit,*) "den"
      call WriteArrayReal(den,size(den),iunit)
      write(iunit,*) "dif2_use"
      call WriteArrayReal(dif2_use,size(dif2_use),iunit)
      write(iunit,*) "dif_use"
      call WriteArrayReal(dif_use,size(dif_use),iunit)
      write(iunit,*) "diffusivwrk"
      call WriteArrayReal(diffusivwrk,size(diffusivwrk),iunit)
      write(iunit,*) "difp_use"
      call WriteArrayReal(difp_use,size(difp_use),iunit)
      write(iunit,*) "dp1"
      write(iunit,*) dp1
      write(iunit,*) "dphi_iy1"
      call WriteArrayReal(dphi_iy1,size(dphi_iy1),iunit)
      write(iunit,*) "dtold"
      write(iunit,*) dtold
      write(iunit,*) "dtreal"
      write(iunit,*) dtreal
      write(iunit,*) "dutm_use"
      call WriteArrayReal(dutm_use,size(dutm_use),iunit)
      write(iunit,*) "dztot"
      call WriteArrayReal(dztot,size(dztot),iunit)
      write(iunit,*) "eeli"
      call WriteArrayReal(eeli,size(eeli),iunit)
      write(iunit,*) "eqp"
      call WriteArrayReal(eqp,size(eqp),iunit)
      write(iunit,*) "eqpg"
      call WriteArrayReal(eqpg,size(eqpg),iunit)
      write(iunit,*) "erliz"
      call WriteArrayReal(erliz,size(erliz),iunit)
      write(iunit,*) "erlrc"
      call WriteArrayReal(erlrc,size(erlrc),iunit)
      write(iunit,*) "eta1"
      call WriteArrayReal(eta1,size(eta1),iunit)
      write(iunit,*) "ex"
      call WriteArrayReal(ex,size(ex),iunit)
      write(iunit,*) "ey"
      call WriteArrayReal(ey,size(ey),iunit)
      write(iunit,*) "fcdif"
      write(iunit,*) fcdif
      write(iunit,*) "fdiaxlb"
      call WriteArrayReal(fdiaxlb,size(fdiaxlb),iunit)
      write(iunit,*) "fdiaxrb"
      call WriteArrayReal(fdiaxrb,size(fdiaxrb),iunit)
      write(iunit,*) "feex"
      call WriteArrayReal(feex,size(feex),iunit)
      write(iunit,*) "feexy"
      call WriteArrayReal(feexy,size(feexy),iunit)
      write(iunit,*) "feey"
      call WriteArrayReal(feey,size(feey),iunit)
      write(iunit,*) "feey4ord"
      call WriteArrayReal(feey4ord,size(feey4ord),iunit)
      write(iunit,*) "feeycbo"
      call WriteArrayReal(feeycbo,size(feeycbo),iunit)
      write(iunit,*) "fegx"
      call WriteArrayReal(fegx,size(fegx),iunit)
      write(iunit,*) "fegy"
      call WriteArrayReal(fegy,size(fegy),iunit)
      write(iunit,*) "feix"
      call WriteArrayReal(feix,size(feix),iunit)
      write(iunit,*) "feixy"
      call WriteArrayReal(feixy,size(feixy),iunit)
      write(iunit,*) "feiy"
      call WriteArrayReal(feiy,size(feiy),iunit)
      write(iunit,*) "feiy4ord"
      call WriteArrayReal(feiy4ord,size(feiy4ord),iunit)
      write(iunit,*) "feiycbo"
      call WriteArrayReal(feiycbo,size(feiycbo),iunit)
      write(iunit,*) "flngx"
      call WriteArrayReal(flngx,size(flngx),iunit)
      write(iunit,*) "flngxy"
      call WriteArrayReal(flngxy,size(flngxy),iunit)
      write(iunit,*) "flngy"
      call WriteArrayReal(flngy,size(flngy),iunit)
      write(iunit,*) "flox"
      call WriteArrayReal(flox,size(flox),iunit)
      write(iunit,*) "floxe"
      call WriteArrayReal(floxe,size(floxe),iunit)
      write(iunit,*) "floxebgt"
      call WriteArrayReal(floxebgt,size(floxebgt),iunit)
      write(iunit,*) "floxg"
      call WriteArrayReal(floxg,size(floxg),iunit)
      write(iunit,*) "floxge"
      call WriteArrayReal(floxge,size(floxge),iunit)
      write(iunit,*) "floxi"
      call WriteArrayReal(floxi,size(floxi),iunit)
      write(iunit,*) "floxibgt"
      call WriteArrayReal(floxibgt,size(floxibgt),iunit)
      write(iunit,*) "floy"
      call WriteArrayReal(floy,size(floy),iunit)
      write(iunit,*) "floye"
      call WriteArrayReal(floye,size(floye),iunit)
      write(iunit,*) "floyg"
      call WriteArrayReal(floyg,size(floyg),iunit)
      write(iunit,*) "floyge"
      call WriteArrayReal(floyge,size(floyge),iunit)
      write(iunit,*) "floyi"
      call WriteArrayReal(floyi,size(floyi),iunit)
      write(iunit,*) "fmihxpt"
      call WriteArrayReal(fmihxpt,size(fmihxpt),iunit)
      write(iunit,*) "fmity"
      call WriteArrayReal(fmity,size(fmity),iunit)
      write(iunit,*) "fmivxpt"
      call WriteArrayReal(fmivxpt,size(fmivxpt),iunit)
      write(iunit,*) "fmix"
      call WriteArrayReal(fmix,size(fmix),iunit)
      write(iunit,*) "fmixy"
      call WriteArrayReal(fmixy,size(fmixy),iunit)
      write(iunit,*) "fmiy"
      call WriteArrayReal(fmiy,size(fmiy),iunit)
      write(iunit,*) "fngx"
      call WriteArrayReal(fngx,size(fngx),iunit)
      write(iunit,*) "fngx4ord"
      call WriteArrayReal(fngx4ord,size(fngx4ord),iunit)
      write(iunit,*) "fngxy"
      call WriteArrayReal(fngxy,size(fngxy),iunit)
      write(iunit,*) "fngy"
      call WriteArrayReal(fngy,size(fngy),iunit)
      write(iunit,*) "fngy4ord"
      call WriteArrayReal(fngy4ord,size(fngy4ord),iunit)
      write(iunit,*) "fngysi"
      call WriteArrayReal(fngysi,size(fngysi),iunit)
      write(iunit,*) "fngyso"
      call WriteArrayReal(fngyso,size(fngyso),iunit)
      write(iunit,*) "fnix"
      call WriteArrayReal(fnix,size(fnix),iunit)
      write(iunit,*) "fnixcb"
      call WriteArrayReal(fnixcb,size(fnixcb),iunit)
      write(iunit,*) "fniy"
      call WriteArrayReal(fniy,size(fniy),iunit)
      write(iunit,*) "fniy4ord"
      call WriteArrayReal(fniy4ord,size(fniy4ord),iunit)
      write(iunit,*) "fniycb"
      call WriteArrayReal(fniycb,size(fniycb),iunit)
      write(iunit,*) "fniycbo"
      call WriteArrayReal(fniycbo,size(fniycbo),iunit)
      write(iunit,*) "fq2"
      call WriteArrayReal(fq2,size(fq2),iunit)
      write(iunit,*) "fq2d"
      call WriteArrayReal(fq2d,size(fq2d),iunit)
      write(iunit,*) "fqp"
      call WriteArrayReal(fqp,size(fqp),iunit)
      write(iunit,*) "fqpsatlb"
      call WriteArrayReal(fqpsatlb,size(fqpsatlb),iunit)
      write(iunit,*) "fqpsatrb"
      call WriteArrayReal(fqpsatrb,size(fqpsatrb),iunit)
      write(iunit,*) "fqx"
      call WriteArrayReal(fqx,size(fqx),iunit)
      write(iunit,*) "fqxb"
      call WriteArrayReal(fqxb,size(fqxb),iunit)
      write(iunit,*) "fqy"
      call WriteArrayReal(fqy,size(fqy),iunit)
      write(iunit,*) "fqya"
      call WriteArrayReal(fqya,size(fqya),iunit)
      write(iunit,*) "fqyae"
      call WriteArrayReal(fqyae,size(fqyae),iunit)
      write(iunit,*) "fqyai"
      call WriteArrayReal(fqyai,size(fqyai),iunit)
      write(iunit,*) "fqyao"
      call WriteArrayReal(fqyao,size(fqyao),iunit)
      write(iunit,*) "fqyb"
      call WriteArrayReal(fqyb,size(fqyb),iunit)
      write(iunit,*) "fqyd"
      call WriteArrayReal(fqyd,size(fqyd),iunit)
      write(iunit,*) "fqydt"
      call WriteArrayReal(fqydt,size(fqydt),iunit)
      write(iunit,*) "fqydti"
      call WriteArrayReal(fqydti,size(fqydti),iunit)
      write(iunit,*) "fqygp"
      call WriteArrayReal(fqygp,size(fqygp),iunit)
      write(iunit,*) "fqym"
      call WriteArrayReal(fqym,size(fqym),iunit)
      write(iunit,*) "fqymi"
      call WriteArrayReal(fqymi,size(fqymi),iunit)
      write(iunit,*) "fqyn"
      call WriteArrayReal(fqyn,size(fqyn),iunit)
      write(iunit,*) "frice"
      call WriteArrayReal(frice,size(frice),iunit)
      write(iunit,*) "frici"
      call WriteArrayReal(frici,size(frici),iunit)
      write(iunit,*) "fricnrl"
      call WriteArrayReal(fricnrl,size(fricnrl),iunit)
      write(iunit,*) "fxe"
      write(iunit,*) fxe
      write(iunit,*) "fxi"
      write(iunit,*) fxi
      write(iunit,*) "ghxpt"
      write(iunit,*) ghxpt
      write(iunit,*) "gpex"
      call WriteArrayReal(gpex,size(gpex),iunit)
      write(iunit,*) "gpey"
      call WriteArrayReal(gpey,size(gpey),iunit)
      write(iunit,*) "gpix"
      call WriteArrayReal(gpix,size(gpix),iunit)
      write(iunit,*) "gpiy"
      call WriteArrayReal(gpiy,size(gpiy),iunit)
      write(iunit,*) "gprx"
      call WriteArrayReal(gprx,size(gprx),iunit)
      write(iunit,*) "gpry"
      call WriteArrayReal(gpry,size(gpry),iunit)
      write(iunit,*) "gradp"
      call WriteArrayReal(gradp,size(gradp),iunit)
      write(iunit,*) "gradt"
      call WriteArrayReal(gradt,size(gradt),iunit)
      write(iunit,*) "gtex"
      call WriteArrayReal(gtex,size(gtex),iunit)
      write(iunit,*) "gtey"
      call WriteArrayReal(gtey,size(gtey),iunit)
      write(iunit,*) "gtix"
      call WriteArrayReal(gtix,size(gtix),iunit)
      write(iunit,*) "gtiy"
      call WriteArrayReal(gtiy,size(gtiy),iunit)
      write(iunit,*) "gvxpt"
      write(iunit,*) gvxpt
      write(iunit,*) "hcxe"
      call WriteArrayReal(hcxe,size(hcxe),iunit)
      write(iunit,*) "hcxg"
      call WriteArrayReal(hcxg,size(hcxg),iunit)
      write(iunit,*) "hcxi"
      call WriteArrayReal(hcxi,size(hcxi),iunit)
      write(iunit,*) "hcxij"
      call WriteArrayReal(hcxij,size(hcxij),iunit)
      write(iunit,*) "hcxineo"
      call WriteArrayReal(hcxineo,size(hcxineo),iunit)
      write(iunit,*) "hcxn"
      call WriteArrayReal(hcxn,size(hcxn),iunit)
      write(iunit,*) "hcye"
      call WriteArrayReal(hcye,size(hcye),iunit)
      write(iunit,*) "hcyg"
      call WriteArrayReal(hcyg,size(hcyg),iunit)
      write(iunit,*) "hcyi"
      call WriteArrayReal(hcyi,size(hcyi),iunit)
      write(iunit,*) "hcyij"
      call WriteArrayReal(hcyij,size(hcyij),iunit)
      write(iunit,*) "hcyn"
      call WriteArrayReal(hcyn,size(hcyn),iunit)
      write(iunit,*) "i1"
      write(iunit,*) i1
      write(iunit,*) "i2"
      write(iunit,*) i2
      write(iunit,*) "i2p"
      write(iunit,*) i2p
      write(iunit,*) "i3"
      write(iunit,*) i3
      write(iunit,*) "i4"
      write(iunit,*) i4
      write(iunit,*) "i5"
      write(iunit,*) i5
      write(iunit,*) "i5m"
      write(iunit,*) i5m
      write(iunit,*) "i6"
      write(iunit,*) i6
      write(iunit,*) "i7"
      write(iunit,*) i7
      write(iunit,*) "i8"
      write(iunit,*) i8
      write(iunit,*) "impradloc"
      call WriteArrayReal(impradloc,size(impradloc),iunit)
      write(iunit,*) "ispwrbcl"
      write(iunit,*) ispwrbcl
      write(iunit,*) "iwalli"
      call WriteArrayReal(iwalli,size(iwalli),iunit)
      write(iunit,*) "iwallo"
      call WriteArrayReal(iwallo,size(iwallo),iunit)
      write(iunit,*) "ixf6"
      write(iunit,*) ixf6
      write(iunit,*) "ixs1"
      write(iunit,*) ixs1
      write(iunit,*) "iyf6"
      write(iunit,*) iyf6
      write(iunit,*) "iys1"
      write(iunit,*) iys1
      write(iunit,*) "j1"
      write(iunit,*) j1
      write(iunit,*) "j1p"
      write(iunit,*) j1p
      write(iunit,*) "j2"
      write(iunit,*) j2
      write(iunit,*) "j2p"
      write(iunit,*) j2p
      write(iunit,*) "j3"
      write(iunit,*) j3
      write(iunit,*) "j4"
      write(iunit,*) j4
      write(iunit,*) "j5"
      write(iunit,*) j5
      write(iunit,*) "j5m"
      write(iunit,*) j5m
      write(iunit,*) "j5p"
      write(iunit,*) j5p
      write(iunit,*) "j6"
      write(iunit,*) j6
      write(iunit,*) "j6p"
      write(iunit,*) j6p
      write(iunit,*) "j7"
      write(iunit,*) j7
      write(iunit,*) "j8"
      write(iunit,*) j8
      write(iunit,*) "k2neo"
      call WriteArrayReal(k2neo,size(k2neo),iunit)
      write(iunit,*) "kappal"
      call WriteArrayReal(kappal,size(kappal),iunit)
      write(iunit,*) "kappar"
      call WriteArrayReal(kappar,size(kappar),iunit)
      write(iunit,*) "kincorlb"
      call WriteArrayReal(kincorlb,size(kincorlb),iunit)
      write(iunit,*) "kincorrb"
      call WriteArrayReal(kincorrb,size(kincorrb),iunit)
      write(iunit,*) "ktneo"
      call WriteArrayReal(ktneo,size(ktneo),iunit)
      write(iunit,*) "kxbohm"
      call WriteArrayReal(kxbohm,size(kxbohm),iunit)
      write(iunit,*) "kxe_use"
      call WriteArrayReal(kxe_use,size(kxe_use),iunit)
      write(iunit,*) "kxi_use"
      call WriteArrayReal(kxi_use,size(kxi_use),iunit)
      write(iunit,*) "kybohm"
      call WriteArrayReal(kybohm,size(kybohm),iunit)
      write(iunit,*) "kye_use"
      call WriteArrayReal(kye_use,size(kye_use),iunit)
      write(iunit,*) "kyi_use"
      call WriteArrayReal(kyi_use,size(kyi_use),iunit)
      write(iunit,*) "lng"
      call WriteArrayReal(lng,size(lng),iunit)
      write(iunit,*) "loglambda"
      call WriteArrayReal(loglambda,size(loglambda),iunit)
      write(iunit,*) "msor"
      call WriteArrayReal(msor,size(msor),iunit)
      write(iunit,*) "msorold"
      call WriteArrayReal(msorold,size(msorold),iunit)
      write(iunit,*) "msorxr"
      call WriteArrayReal(msorxr,size(msorxr),iunit)
      write(iunit,*) "msorxrold"
      call WriteArrayReal(msorxrold,size(msorxrold),iunit)
      write(iunit,*) "na"
      call WriteArrayReal(na,size(na),iunit)
      write(iunit,*) "ncrhs"
      write(iunit,*) ncrhs
      write(iunit,*) "ne"
      call WriteArrayReal(ne,size(ne),iunit)
      write(iunit,*) "netap"
      call WriteArrayReal(netap,size(netap),iunit)
      write(iunit,*) "ney0"
      call WriteArrayReal(ney0,size(ney0),iunit)
      write(iunit,*) "ney1"
      call WriteArrayReal(ney1,size(ney1),iunit)
      write(iunit,*) "nfsp"
      write(iunit,*) nfsp
      write(iunit,*) "ng"
      call WriteArrayReal(ng,size(ng),iunit)
      write(iunit,*) "ngy0"
      call WriteArrayReal(ngy0,size(ngy0),iunit)
      write(iunit,*) "ngy1"
      call WriteArrayReal(ngy1,size(ngy1),iunit)
      write(iunit,*) "ni"
      call WriteArrayReal(ni,size(ni),iunit)
      write(iunit,*) "nit"
      call WriteArrayReal(nit,size(nit),iunit)
      write(iunit,*) "nity0"
      call WriteArrayReal(nity0,size(nity0),iunit)
      write(iunit,*) "nity1"
      call WriteArrayReal(nity1,size(nity1),iunit)
      write(iunit,*) "nixpt"
      call WriteArrayReal(nixpt,size(nixpt),iunit)
      write(iunit,*) "niy0"
      call WriteArrayReal(niy0,size(niy0),iunit)
      write(iunit,*) "niy0s"
      call WriteArrayReal(niy0s,size(niy0s),iunit)
      write(iunit,*) "niy1"
      call WriteArrayReal(niy1,size(niy1),iunit)
      write(iunit,*) "niy1s"
      call WriteArrayReal(niy1s,size(niy1s),iunit)
      write(iunit,*) "nm"
      call WriteArrayReal(nm,size(nm),iunit)
      write(iunit,*) "nratio"
      call WriteArrayReal(nratio,size(nratio),iunit)
      write(iunit,*) "ntau"
      call WriteArrayReal(ntau,size(ntau),iunit)
      write(iunit,*) "nucx"
      call WriteArrayReal(nucx,size(nucx),iunit)
      write(iunit,*) "nucxi"
      call WriteArrayReal(nucxi,size(nucxi),iunit)
      write(iunit,*) "nuelg"
      call WriteArrayReal(nuelg,size(nuelg),iunit)
      write(iunit,*) "nueli"
      call WriteArrayReal(nueli,size(nueli),iunit)
      write(iunit,*) "nuii"
      call WriteArrayReal(nuii,size(nuii),iunit)
      write(iunit,*) "nuiistar"
      call WriteArrayReal(nuiistar,size(nuiistar),iunit)
      write(iunit,*) "nuix"
      call WriteArrayReal(nuix,size(nuix),iunit)
      write(iunit,*) "nuiz"
      call WriteArrayReal(nuiz,size(nuiz),iunit)
      write(iunit,*) "nurc"
      call WriteArrayReal(nurc,size(nurc),iunit)
      write(iunit,*) "nuvl"
      call WriteArrayReal(nuvl,size(nuvl),iunit)
      write(iunit,*) "nz2"
      call WriteArrayReal(nz2,size(nz2),iunit)
      write(iunit,*) "nzloc"
      call WriteArrayReal(nzloc,size(nzloc),iunit)
      write(iunit,*) "openbox"
      write(iunit,*) openbox
      write(iunit,*) "parvis"
      call WriteArrayReal(parvis,size(parvis),iunit)
      write(iunit,*) "pg"
      call WriteArrayReal(pg,size(pg),iunit)
      write(iunit,*) "pgy0"
      call WriteArrayReal(pgy0,size(pgy0),iunit)
      write(iunit,*) "pgy1"
      call WriteArrayReal(pgy1,size(pgy1),iunit)
      write(iunit,*) "phi"
      call WriteArrayReal(phi,size(phi),iunit)
      write(iunit,*) "phiv"
      call WriteArrayReal(phiv,size(phiv),iunit)
      write(iunit,*) "phiy0"
      call WriteArrayReal(phiy0,size(phiy0),iunit)
      write(iunit,*) "phiy0s"
      call WriteArrayReal(phiy0s,size(phiy0s),iunit)
      write(iunit,*) "phiy1"
      call WriteArrayReal(phiy1,size(phiy1),iunit)
      write(iunit,*) "phiy1s"
      call WriteArrayReal(phiy1s,size(phiy1s),iunit)
      write(iunit,*) "pr"
      call WriteArrayReal(pr,size(pr),iunit)
      write(iunit,*) "prad"
      call WriteArrayReal(prad,size(prad),iunit)
      write(iunit,*) "pradc"
      call WriteArrayReal(pradc,size(pradc),iunit)
      write(iunit,*) "pradcff"
      call WriteArrayReal(pradcff,size(pradcff),iunit)
      write(iunit,*) "pradhyd"
      call WriteArrayReal(pradhyd,size(pradhyd),iunit)
      write(iunit,*) "pradz"
      call WriteArrayReal(pradz,size(pradz),iunit)
      write(iunit,*) "pradzc"
      call WriteArrayReal(pradzc,size(pradzc),iunit)
      write(iunit,*) "pre"
      call WriteArrayReal(pre,size(pre),iunit)
      write(iunit,*) "prev"
      call WriteArrayReal(prev,size(prev),iunit)
      write(iunit,*) "pri"
      call WriteArrayReal(pri,size(pri),iunit)
      write(iunit,*) "priv"
      call WriteArrayReal(priv,size(priv),iunit)
      write(iunit,*) "priy0"
      call WriteArrayReal(priy0,size(priy0),iunit)
      write(iunit,*) "priy1"
      call WriteArrayReal(priy1,size(priy1),iunit)
      write(iunit,*) "prtv"
      call WriteArrayReal(prtv,size(prtv),iunit)
      write(iunit,*) "psor"
      call WriteArrayReal(psor,size(psor),iunit)
      write(iunit,*) "psor_tmpov"
      call WriteArrayReal(psor_tmpov,size(psor_tmpov),iunit)
      write(iunit,*) "psorbgg"
      call WriteArrayReal(psorbgg,size(psorbgg),iunit)
      write(iunit,*) "psorbgz"
      call WriteArrayReal(psorbgz,size(psorbgz),iunit)
      write(iunit,*) "psorc"
      call WriteArrayReal(psorc,size(psorc),iunit)
      write(iunit,*) "psorcxg"
      call WriteArrayReal(psorcxg,size(psorcxg),iunit)
      write(iunit,*) "psorcxgc"
      call WriteArrayReal(psorcxgc,size(psorcxgc),iunit)
      write(iunit,*) "psordis"
      call WriteArrayReal(psordis,size(psordis),iunit)
      write(iunit,*) "psorg"
      call WriteArrayReal(psorg,size(psorg),iunit)
      write(iunit,*) "psorgc"
      call WriteArrayReal(psorgc,size(psorgc),iunit)
      write(iunit,*) "psori"
      call WriteArrayReal(psori,size(psori),iunit)
      write(iunit,*) "psorold"
      call WriteArrayReal(psorold,size(psorold),iunit)
      write(iunit,*) "psorrg"
      call WriteArrayReal(psorrg,size(psorrg),iunit)
      write(iunit,*) "psorrgc"
      call WriteArrayReal(psorrgc,size(psorrgc),iunit)
      write(iunit,*) "psorxr"
      call WriteArrayReal(psorxr,size(psorxr),iunit)
      write(iunit,*) "psorxrc"
      call WriteArrayReal(psorxrc,size(psorxrc),iunit)
      write(iunit,*) "psorxrold"
      call WriteArrayReal(psorxrold,size(psorxrold),iunit)
      write(iunit,*) "pwrebkg"
      call WriteArrayReal(pwrebkg,size(pwrebkg),iunit)
      write(iunit,*) "pwribkg"
      call WriteArrayReal(pwribkg,size(pwribkg),iunit)
      write(iunit,*) "pwrze"
      call WriteArrayReal(pwrze,size(pwrze),iunit)
      write(iunit,*) "pwrzec"
      call WriteArrayReal(pwrzec,size(pwrzec),iunit)
      write(iunit,*) "q2cd"
      call WriteArrayReal(q2cd,size(q2cd),iunit)
      write(iunit,*) "qipar"
      call WriteArrayReal(qipar,size(qipar),iunit)
      write(iunit,*) "resco"
      call WriteArrayReal(resco,size(resco),iunit)
      write(iunit,*) "resee"
      call WriteArrayReal(resee,size(resee),iunit)
      write(iunit,*) "reseg"
      call WriteArrayReal(reseg,size(reseg),iunit)
      write(iunit,*) "resei"
      call WriteArrayReal(resei,size(resei),iunit)
      write(iunit,*) "resmo"
      call WriteArrayReal(resmo,size(resmo),iunit)
      write(iunit,*) "resng"
      call WriteArrayReal(resng,size(resng),iunit)
      write(iunit,*) "resphi"
      call WriteArrayReal(resphi,size(resphi),iunit)
      write(iunit,*) "rtau"
      call WriteArrayReal(rtau,size(rtau),iunit)
      write(iunit,*) "rtaue"
      call WriteArrayReal(rtaue,size(rtaue),iunit)
      write(iunit,*) "rtaux"
      call WriteArrayReal(rtaux,size(rtaux),iunit)
      write(iunit,*) "rtauy"
      call WriteArrayReal(rtauy,size(rtauy),iunit)
      write(iunit,*) "seec"
      call WriteArrayReal(seec,size(seec),iunit)
      write(iunit,*) "seev"
      call WriteArrayReal(seev,size(seev),iunit)
      write(iunit,*) "seg_ue"
      call WriteArrayReal(seg_ue,size(seg_ue),iunit)
      write(iunit,*) "seic"
      call WriteArrayReal(seic,size(seic),iunit)
      write(iunit,*) "seiv"
      call WriteArrayReal(seiv,size(seiv),iunit)
      write(iunit,*) "smoc"
      call WriteArrayReal(smoc,size(smoc),iunit)
      write(iunit,*) "smov"
      call WriteArrayReal(smov,size(smov),iunit)
      write(iunit,*) "sng_ue"
      call WriteArrayReal(sng_ue,size(sng_ue),iunit)
      write(iunit,*) "snic"
      call WriteArrayReal(snic,size(snic),iunit)
      write(iunit,*) "sniv"
      call WriteArrayReal(sniv,size(sniv),iunit)
      write(iunit,*) "sputflxlb"
      call WriteArrayReal(sputflxlb,size(sputflxlb),iunit)
      write(iunit,*) "sputflxpf"
      call WriteArrayReal(sputflxpf,size(sputflxpf),iunit)
      write(iunit,*) "sputflxrb"
      call WriteArrayReal(sputflxrb,size(sputflxrb),iunit)
      write(iunit,*) "sputflxw"
      call WriteArrayReal(sputflxw,size(sputflxw),iunit)
      write(iunit,*) "sxyxpt"
      write(iunit,*) sxyxpt
      write(iunit,*) "te"
      call WriteArrayReal(te,size(te),iunit)
      write(iunit,*) "tempa"
      call WriteArrayReal(tempa,size(tempa),iunit)
      write(iunit,*) "tev"
      call WriteArrayReal(tev,size(tev),iunit)
      write(iunit,*) "tey0"
      call WriteArrayReal(tey0,size(tey0),iunit)
      write(iunit,*) "tey1"
      call WriteArrayReal(tey1,size(tey1),iunit)
      write(iunit,*) "tg"
      call WriteArrayReal(tg,size(tg),iunit)
      write(iunit,*) "tgy0"
      call WriteArrayReal(tgy0,size(tgy0),iunit)
      write(iunit,*) "tgy1"
      call WriteArrayReal(tgy1,size(tgy1),iunit)
      write(iunit,*) "ti"
      call WriteArrayReal(ti,size(ti),iunit)
      write(iunit,*) "tiv"
      call WriteArrayReal(tiv,size(tiv),iunit)
      write(iunit,*) "tiy0"
      call WriteArrayReal(tiy0,size(tiy0),iunit)
      write(iunit,*) "tiy0s"
      call WriteArrayReal(tiy0s,size(tiy0s),iunit)
      write(iunit,*) "tiy1"
      call WriteArrayReal(tiy1,size(tiy1),iunit)
      write(iunit,*) "tiy1s"
      call WriteArrayReal(tiy1s,size(tiy1s),iunit)
      write(iunit,*) "totb2val"
      write(iunit,*) totb2val
      write(iunit,*) "totfeexl"
      call WriteArrayReal(totfeexl,size(totfeexl),iunit)
      write(iunit,*) "totfeexr"
      call WriteArrayReal(totfeexr,size(totfeexr),iunit)
      write(iunit,*) "totfeixl"
      call WriteArrayReal(totfeixl,size(totfeixl),iunit)
      write(iunit,*) "totfeixr"
      call WriteArrayReal(totfeixr,size(totfeixr),iunit)
      write(iunit,*) "travis"
      call WriteArrayReal(travis,size(travis),iunit)
      write(iunit,*) "trax_use"
      call WriteArrayReal(trax_use,size(trax_use),iunit)
      write(iunit,*) "tray_use"
      call WriteArrayReal(tray_use,size(tray_use),iunit)
      write(iunit,*) "ttimpfe"
      write(iunit,*) ttimpfe
      write(iunit,*) "ttimpjf"
      write(iunit,*) ttimpjf
      write(iunit,*) "ttngfd2"
      write(iunit,*) ttngfd2
      write(iunit,*) "ttngfxy"
      write(iunit,*) ttngfxy
      write(iunit,*) "ttngxlog"
      write(iunit,*) ttngxlog
      write(iunit,*) "ttngylog"
      write(iunit,*) ttngylog
      write(iunit,*) "ttnpg"
      write(iunit,*) ttnpg
      write(iunit,*) "ttotfe"
      write(iunit,*) ttotfe
      write(iunit,*) "ttotjf"
      write(iunit,*) ttotjf
      write(iunit,*) "up"
      call WriteArrayReal(up,size(up),iunit)
      write(iunit,*) "upe"
      call WriteArrayReal(upe,size(upe),iunit)
      write(iunit,*) "upi"
      call WriteArrayReal(upi,size(upi),iunit)
      write(iunit,*) "upxpt"
      call WriteArrayReal(upxpt,size(upxpt),iunit)
      write(iunit,*) "uu"
      call WriteArrayReal(uu,size(uu),iunit)
      write(iunit,*) "uug"
      call WriteArrayReal(uug,size(uug),iunit)
      write(iunit,*) "uup"
      call WriteArrayReal(uup,size(uup),iunit)
      write(iunit,*) "uz"
      call WriteArrayReal(uz,size(uz),iunit)
      write(iunit,*) "v2"
      call WriteArrayReal(v2,size(v2),iunit)
      write(iunit,*) "v2cb"
      call WriteArrayReal(v2cb,size(v2cb),iunit)
      write(iunit,*) "v2cd"
      call WriteArrayReal(v2cd,size(v2cd),iunit)
      write(iunit,*) "v2ce"
      call WriteArrayReal(v2ce,size(v2ce),iunit)
      write(iunit,*) "v2dd"
      call WriteArrayReal(v2dd,size(v2dd),iunit)
      write(iunit,*) "v2rd"
      call WriteArrayReal(v2rd,size(v2rd),iunit)
      write(iunit,*) "v2xgp"
      call WriteArrayReal(v2xgp,size(v2xgp),iunit)
      write(iunit,*) "ve2cb"
      call WriteArrayReal(ve2cb,size(ve2cb),iunit)
      write(iunit,*) "ve2cd"
      call WriteArrayReal(ve2cd,size(ve2cd),iunit)
      write(iunit,*) "vex"
      call WriteArrayReal(vex,size(vex),iunit)
      write(iunit,*) "vey"
      call WriteArrayReal(vey,size(vey),iunit)
      write(iunit,*) "veycb"
      call WriteArrayReal(veycb,size(veycb),iunit)
      write(iunit,*) "veycp"
      call WriteArrayReal(veycp,size(veycp),iunit)
      write(iunit,*) "visvol_q"
      call WriteArrayReal(visvol_q,size(visvol_q),iunit)
      write(iunit,*) "visvol_v"
      call WriteArrayReal(visvol_v,size(visvol_v),iunit)
      write(iunit,*) "visx"
      call WriteArrayReal(visx,size(visx),iunit)
      write(iunit,*) "visxneo"
      call WriteArrayReal(visxneo,size(visxneo),iunit)
      write(iunit,*) "visy"
      call WriteArrayReal(visy,size(visy),iunit)
      write(iunit,*) "visyxpt"
      call WriteArrayReal(visyxpt,size(visyxpt),iunit)
      write(iunit,*) "vsoree"
      call WriteArrayReal(vsoree,size(vsoree),iunit)
      write(iunit,*) "vsoreec"
      call WriteArrayReal(vsoreec,size(vsoreec),iunit)
      write(iunit,*) "vy"
      call WriteArrayReal(vy,size(vy),iunit)
      write(iunit,*) "vy_cft"
      call WriteArrayReal(vy_cft,size(vy_cft),iunit)
      write(iunit,*) "vy_use"
      call WriteArrayReal(vy_use,size(vy_use),iunit)
      write(iunit,*) "vyavis"
      call WriteArrayReal(vyavis,size(vyavis),iunit)
      write(iunit,*) "vycb"
      call WriteArrayReal(vycb,size(vycb),iunit)
      write(iunit,*) "vyce"
      call WriteArrayReal(vyce,size(vyce),iunit)
      write(iunit,*) "vycf"
      call WriteArrayReal(vycf,size(vycf),iunit)
      write(iunit,*) "vycp"
      call WriteArrayReal(vycp,size(vycp),iunit)
      write(iunit,*) "vycr"
      call WriteArrayReal(vycr,size(vycr),iunit)
      write(iunit,*) "vydd"
      call WriteArrayReal(vydd,size(vydd),iunit)
      write(iunit,*) "vyg"
      call WriteArrayReal(vyg,size(vyg),iunit)
      write(iunit,*) "vygp"
      call WriteArrayReal(vygp,size(vygp),iunit)
      write(iunit,*) "vygtan"
      call WriteArrayReal(vygtan,size(vygtan),iunit)
      write(iunit,*) "vyhxpt"
      call WriteArrayReal(vyhxpt,size(vyhxpt),iunit)
      write(iunit,*) "vyrd"
      call WriteArrayReal(vyrd,size(vyrd),iunit)
      write(iunit,*) "vytan"
      call WriteArrayReal(vytan,size(vytan),iunit)
      write(iunit,*) "vyte_cft"
      call WriteArrayReal(vyte_cft,size(vyte_cft),iunit)
      write(iunit,*) "vyti_cft"
      call WriteArrayReal(vyti_cft,size(vyti_cft),iunit)
      write(iunit,*) "vyvxpt"
      call WriteArrayReal(vyvxpt,size(vyvxpt),iunit)
      write(iunit,*) "w"
      call WriteArrayReal(w,size(w),iunit)
      write(iunit,*) "w0"
      call WriteArrayReal(w0,size(w0),iunit)
      write(iunit,*) "w1"
      call WriteArrayReal(w1,size(w1),iunit)
      write(iunit,*) "w2"
      call WriteArrayReal(w2,size(w2),iunit)
      write(iunit,*) "w3"
      call WriteArrayReal(w3,size(w3),iunit)
      write(iunit,*) "wjdote"
      call WriteArrayReal(wjdote,size(wjdote),iunit)
      write(iunit,*) "wvh"
      call WriteArrayReal(wvh,size(wvh),iunit)
      write(iunit,*) "xcnearlb"
      write(iunit,*) xcnearlb
      write(iunit,*) "xcnearrb"
      write(iunit,*) xcnearrb
      write(iunit,*) "yld_carbi"
      call WriteArrayReal(yld_carbi,size(yld_carbi),iunit)
      write(iunit,*) "yld_carbo"
      call WriteArrayReal(yld_carbo,size(yld_carbo),iunit)
      write(iunit,*) "yldot_pert"
      call WriteArrayReal(yldot_pert,size(yldot_pert),iunit)
      write(iunit,*) "zcoef"
      write(iunit,*) zcoef
      write(iunit,*) "zeff"
      call WriteArrayReal(zeff,size(zeff),iunit)
      write(iunit,*) "znot"
      call WriteArrayReal(znot,size(znot),iunit)
      close(iunit)
      end subroutine DebugHelper
!
!    subroutine TestParallelPandf1(neq,time,yl,yldot)
!    use omp_lib
!    use OMPPandf1
!    use OMPPandf1Settings
!    use OMPJacSettings
!    use Dim,only:nx
!    use Selec, only: xlinc,xrinc
!    integer xlinc_bkp,xrinc_bkp,iv
!    integer,intent(in)::neq
!    real,intent(in)::yl(*)
!    real,intent(out)::yldot(*)
!    real,intent(in)::time
!    real::yldotcopy(1:neq)
!    real::TimeStart,TimeEnd
!    real::walltime
!    real yldotsave(neq)
!    !call omp_set_num_threads(int(OMPPandf_Nthreads,kind=4))
!
!integer::ith,OMPic(Nthreads),OMPlinc(Nthreads),OMPrinc(Nthreads),OMPivthread(1:neq)
!
!    call OMPSplitIndexPandf(0,nx+1,Nthreads,OMPic,OMPlinc,OMPrinc,OMPivthread,neq)
!
!    !if (OMPJacVerbose.gt.0) then
!        write(*,*)' *OMPPandf1* nx=',nx
!        write(*,*)' *OMPPandf1* Ic(ith),linc(ith),rinc(ith) ***'
!        do ith=1,Nthreads
!            write(*,'(a,I3,a,I7,I7,I7)') '  *    ithread ', ith,':',OMPic(ith),OMPlinc(ith),OMPrinc(ith)
!        enddo
!      walltime=omp_get_wtime()
!!      if (OMPCheckThreadedPandf.gt.0) then
!!      TimingPandf=0
!!      TimingParaPandf=1
!!      TimingParaConvert=1
!!      endif
!      call pandf1 (-1, -1, 0, neq, time, yl, yldot)
!
!      TimeParallelPandf=omp_get_wtime()-walltime+TimeParallelPandf
!      write(*,*) "Timing Pandf1 serial",TimeParallelPandf
!      xlinc_bkp=xlinc
!      xrinc_bkp=xrinc
!      do ith=1,Nthreads
!        xlinc=OMPlinc(ith)
!        xrinc=OMPrinc(ith)
!        call pandf1 (ith, -1, 0, neq, time, yl, yldotcopy)
!    do iv=1,neq
!      if (OMPivthread(iv)==ith) then
!       yldotsave(iv)=yldotcopy(iv)
!      endif
!      enddo
!
!      enddo
!        xlinc=xlinc_bkp
!        xrinc=xrinc_bkp
!
!
!
!      call Compare(yldot,yldotsave,neq)
!      write(*,*) "yldot and yldotsave are identical"
!!      call omp_set_num_threads(int(Nthreads,kind=4))
!!
!!      if (OMPCheckThreadedPandf.gt.0.and.OMPParallelPandf.gt.0) then
!!      if (OMPThreadedPandfVerbose.gt.0) write(*,*) 'pandf checked'
!!      OMPParallelPandf=0
!!      walltime=omp_get_wtime()
!!      TimingPandf=1
!!      TimingParaPandf=0
!!      TimingParaConvert=0
!!      call pandf1 (-1, -1, 0, neq, time, yl, yldotsave)
!!      TimingPandf=0
!!      TimeSerialPandf=omp_get_wtime()-walltime+TimeSerialPandf
!!      OMPParallelPandf=1
!!
!
!
!!      endif
!
!
!end subroutine TestParallelPandf1
!

subroutine Compare(yldot,yldotsave,neq)
integer:: iv,neq
real yldot(neq+2),yldotsave(neq+2)
            do iv=1,neq
                if (abs(yldotsave(iv)-yldot(iv)).gt.1e-14) then
                    if (max(abs(yldot(iv)),abs(yldotsave(iv)))>0) then
                    if (abs(yldotsave(iv)-yldot(iv))/max(abs(yldot(iv)),abs(yldotsave(iv)))>1e-14) then
                    write(*,*) '>>>>',iv,yldotsave(iv),yldot(iv),abs(yldotsave(iv)-yldot(iv))/max(abs(yldot(iv)),abs(yldotsave(iv)))
                    call xerrab('stop')
                    endif
                    else
                    write(*,*) '>>>>',iv,0
                    endif
                    !call xerrab('diff in rhsnk')
                endif
            enddo



end subroutine Compare



subroutine AddTimeDerivative(neq,yl,yldot)
use Compla,only: zi
use UEpar,only:isbcwdt,isnionxy,isuponxy,isteonxy,istionxy,isngonxy,isphionxy,ineudif,fdtnixy,&
               fdtupxy,fdttexy,fdttixy,fdtngxy,fdttgxy,fdtphixy,istgonxy
use Dim, only:  nusp,nisp,ngsp,nx,ny
use Time_dep_nwt,only:dtreal,ylodt,dtuse,dtphi
use Indexes,only: idxn,idxg,idxu,idxti,idxte,idxphi,idxtg
use Ynorm,only:n0,n0g
!use Share    # geometry,nxc,isnonog,cutlo
use Indices_domain_dcl,only: ixmnbcl,ixmxbcl,iymnbcl,iymxbcl
real,intent(inout)::yldot(*)
real,intent(in)::yl(*)
integer,intent(in):: neq
!use Xpoint_indices      # ixpt1,ixpt2,iysptrx
integer j2l,j5l,i2l,i5l,ifld,igsp,iv,ix,iy,iv1
if (isbcwdt .eq. 0) then  ! omit b.c. eqns
!!!MER   NOTE: what about internal guard cells (for dnbot,dnull,limiter) ???
            j2l = 1
            j5l = ny
            i2l = 1
            i5l = nx
         else                      ! include b.c. eqns
            j2l = (1-iymnbcl)
            j5l = ny+1-(1-iymxbcl)
            i2l = (1-ixmnbcl)
            i5l = nx+1-(1-ixmxbcl)
         endif
         do iy = j2l, j5l    !if j2l=j2, etc., omit the boundary equations
            do ix = i2l, i5l
              do ifld = 1, nisp
                if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  yldot(iv) = (1.-fdtnixy(ix,iy,ifld))*yldot(iv)
                  if(zi(ifld).eq.0. .and. ineudif.eq.3) then
                    yldot(iv) = yldot(iv) - (1/n0(ifld))* (exp(yl(iv))-exp(ylodt(iv)))/dtuse(iv)
                  else
                    yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                  endif
                endif
              enddo
               if(ix.ne.nx+2*isbcwdt) then
                              ! nx test - for algebr. eq. unless isbcwdt=1
                  do ifld = 1, nusp
                    if(isuponxy(ix,iy,ifld).eq.1) then
                      iv = idxu(ix,iy,ifld)
                      yldot(iv) = (1.-fdtupxy(ix,iy,ifld))*yldot(iv)
                      yldot(iv) = yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                    endif
                  enddo
               endif
               if (isteonxy(ix,iy) == 1) then
                 iv =  idxte(ix,iy)
                 yldot(iv) = (1.-fdttexy(ix,iy))*yldot(iv)
                 yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtuse(iv)
               endif
               if (istionxy(ix,iy) == 1) then
                 iv1 = idxti(ix,iy)
                 yldot(iv1) = (1.-fdttixy(ix,iy))*yldot(iv1)
                 yldot(iv1)=yldot(iv1) - (yl(iv1)-ylodt(iv1))/dtuse(iv1)
               endif
               do igsp = 1, ngsp
                  if(isngonxy(ix,iy,igsp).eq.1) then
                     iv = idxg(ix,iy,igsp)
                     yldot(iv) = (1.-fdtngxy(ix,iy,igsp))*yldot(iv)
                     if(ineudif.eq.3) then
                       yldot(iv) = yldot(iv) - (1/n0g(igsp))*(exp(yl(iv))-exp(ylodt(iv)))/dtuse(iv)
                     else
                       yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                     endif
                  endif
               enddo
               do igsp = 1, ngsp
                  if(istgonxy(ix,iy,igsp).eq.1) then
                     iv = idxtg(ix,iy,igsp)
                     yldot(iv) = (1.-fdttgxy(ix,iy,igsp))*yldot(iv)
                     yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                  endif
               enddo
               if (isphionxy(ix,iy).eq.1 .and. isbcwdt.eq.1) then
                  iv = idxphi(ix,iy)
                  yldot(iv) = (1.-fdtphixy(ix,iy))*yldot(iv)
                  yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtuse(iv)
               endif

            enddo
         enddo

!...  Now do an additional relaxation of the potential equations with
!...  timestep dtphi
        if (dtphi < 1e10) then
          do iy = 0, ny+1
            do ix = 0, nx+1
              if (isphionxy(ix,iy) == 1) then
                iv = idxphi(ix,iy)
                yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtphi
              endif
            enddo
          enddo
        endif
end subroutine

        subroutine tick(t)
        implicit none
            integer, intent(OUT) :: t

            call system_clock(t)
        end subroutine tick

        real function tock(t)
         implicit none
            integer, intent(in) :: t
            integer :: now, clock_rate

            call system_clock(now,clock_rate)

            tock = real(now - t)/real(clock_rate)
        end function tock

subroutine PrintTimingPandf()
    use TimingPandf
    write(*,*) '----- Timing Pandf in evalrhs mode----'
    write(*,*) 'TimeParallelPandf',TimeParallelPandf
    write(*,*) 'TimeserialPandf',TimeSerialPandf
    write(*,*) '--------------------------------------'
    write(*,*)'TimeConvert0', TimeConvert0
    write(*,*)'TimeConvert1', TimeConvert1
    write(*,*) '--------------------------------------'
    write(*,*)'TimeCopy0', TimeCopy0
    write(*,*)'TimeCopy1', TimeCopy1
    write(*,*)'TimeCopy2', TimeCopy2
    write(*,*) '--------------------------------------'
write(*,*)'TimeBlock1', TimeBlock1
write(*,*)'TimeBlock2' ,TimeBlock2
write(*,*)'TimeBlock3', TimeBlock3
write(*,*)'TimeBlock4', TimeBlock4
write(*,*)'TimeBlock5', TimeBlock5
write(*,*)'TimeBlock6', TimeBlock6
write(*,*)'TimeBlock7', TimeBlock7
write(*,*)'TimeBlock8', TimeBlock8
write(*,*)'TimeBlock9', TimeBlock9
write(*,*)'TimeBlock10', TimeBlock10
write(*,*)'TimeBlock11' ,TimeBlock11
write(*,*)'TimeBlock12', TimeBlock12
write(*,*)'TimeBlock13' ,TimeBlock13
write(*,*)'TimeBlock14', TimeBlock14
write(*,*)'TimeBlock15', TimeBlock15
write(*,*)'TimeBlock16', TimeBlock16
write(*,*) '--------------------------------------'
write(*,*)'TimeParaConv1', TimeParaConv1
write(*,*)'TimeParaBlock5', TimeParaBlock5
write(*,*)'TimeParaBlock8', TimeParaBlock8
write(*,*) '--------------------------------------'
write(*,*)'TimeConv0', TimeConv0
write(*,*)'TimeConv1', TimeConv1
write(*,*)'TimeConv2' ,TimeConv2
write(*,*)'TimeConv3', TimeConv3
write(*,*)'TimeConv4', TimeConv4
write(*,*)'TimeConv5', TimeConv5
end subroutine PrintTimingPandf

subroutine OMPSplitIndexPandf(ieqmin,ieqmax,Nthreads,ic,linc,rinc,ivthread,neq)
    Use Indexes,only: igyl
    implicit none
    integer,intent(in) ::ieqmin,ieqmax,Nthreads,neq
    integer,intent(out)::ic(Nthreads),linc(Nthreads),rinc(Nthreads),ivthread(1:neq)
    integer :: ihigh(Nthreads),ilow(Nthreads),ixthread(ieqmin:ieqmax)
    integer:: Nsize(Nthreads),Msize,R,i,imax,iv,lpad=1,rpad=1

    if (ieqmax-ieqmin+1<2) call xerrab('Number of equations to solve <2')

        !        if (OMPLoadWeight.eq.1) then

        !do i=1,Nthreads
         !   if (weight(i)<=0) call xerrab('OMPSplitIndex: weight <0')
            !write(*,*) weight(i)
        !enddo

        ! Normalized weights
        !weight(1:Nthreads)=weight(1:Nthreads)/sum(weight(1:Nthreads))*real(Nthreads)
        do i=1,Nthreads
            !Nsize(i)=int(real((ieqmax-ieqmin+1)/Nthreads)*weight(i))
            Nsize(i)=int(real((ieqmax-ieqmin+1)/Nthreads))
        enddo
        write(*,*) Nsize

        do i=1,Nthreads
            if (Nsize(i)<1) call xerrab('Nsize<1')
            !if (Nsize(i)<2) Nsize(i)=Nsize(i)+1
        enddo

        if (ieqmax-ieqmin+1.ne.sum(Nsize)) then
            imax=1
            do i=2,Nthreads
                if (Nsize(i)>Nsize(i-1)) then
                    imax=i
                endif
            enddo
            Nsize(imax) = Nsize(imax) + ((ieqmax-ieqmin+1)-sum(Nsize))
        endif
        !write(*,*) Nsize,neq,sum(Nsize)
        if (ieqmax-ieqmin+1.ne.sum(Nsize)) call xerrab('Nsize .ne. ieqmax-ieqmin+1')

        ilow(1)=ieqmin
        ihigh(1)=ieqmin+Nsize(1)-1
        do i=2,Nthreads
            ilow(i)=ihigh(i-1)+1
            ihigh(i)=ilow(i)+Nsize(i)-1
        enddo
        write(*,*) ilow
        write(*,*) ihigh
        if (ihigh(Nthreads).ne.ieqmax) call xerrab('ihigh(Nthreads)!=ieqmax')

        ic=ilow+int((ihigh-ilow)/2)
        linc=ic-ilow+lpad
        rinc=ihigh-ic+rpad
        do iv=1,neq
            do i=1,Nthreads
                if (igyl(iv,1).le.ihigh(i) .and. igyl(iv,1).ge.ilow(i)) then
                    ivthread(iv)=i
                endif
            enddo
        enddo

end subroutine OMPSplitIndexPandf

subroutine OMPSplitIndexyPandf(ieqmin,ieqmax,Nthreads,ic,inc,ivthread,neq)
    Use Indexes,only: igyl
    use OMPPandf1Settings,only: OMPPandf1Padyinc
    implicit none
    integer,intent(in) ::ieqmin,ieqmax,Nthreads,neq
    integer,intent(out)::ic(Nthreads),inc(Nthreads),ivthread(1:neq)
    integer :: ihigh(Nthreads),ilow(Nthreads),ixthread(ieqmin:ieqmax)
    integer:: Nsize(Nthreads),Msize,R,i,imax,iv,imin,iend,istart
    if (ieqmax-ieqmin+1<2) call xerrab('Number of equations to solve <2')

        !        if (OMPLoadWeight.eq.1) then

        !do i=1,Nthreads
         !   if (weight(i)<=0) call xerrab('OMPSplitIndex: weight <0')
            !write(*,*) weight(i)
        !enddo

        ! Normalized weights
        !weight(1:Nthreads)=weight(1:Nthreads)/sum(weight(1:Nthreads))*real(Nthreads)
        do i=1,Nthreads
            !Nsize(i)=int(real((ieqmax-ieqmin+1)/Nthreads)*weight(i))
            Nsize(i)=int(real((ieqmax-ieqmin+1)/Nthreads))
        enddo

        do i=1,Nthreads
            if (Nsize(i)<1) call xerrab('Nsize<1')
        enddo
        if (Nthreads.ge.4) then
        istart=2
        iend=Nthreads-1
        else
        istart=1
        iend=Nthreads
        endif
imin=istart
        do while (ieqmax-ieqmin+1.gt.sum(Nsize))
                Nsize(imin)=Nsize(imin)+1
                if (imin.ge.iend) then
                    imin=istart
                else
                    imin=imin+1
                endif
        enddo

        if (ieqmax-ieqmin+1.ne.sum(Nsize)) then
           write(*,*) Nsize,sum(Nsize)
           call xerrab('OMPPandf1:Nsize .ne. ieqmax-ieqmin+1')
       endif
        ilow(1)=ieqmin
        ihigh(1)=ieqmin+Nsize(1)-1
        do i=2,Nthreads
            ilow(i)=ihigh(i-1)+1
            ihigh(i)=ilow(i)+Nsize(i)-1

        enddo
        if (ihigh(Nthreads).ne.ieqmax) call xerrab('ihigh(Nthreads)!=ieqmax')

        ic=ilow+int((ihigh-ilow)/2)
        do i=1,Nthreads
           !write(*,*) max(ihigh(i)-ic(i),ic(i)-ilow(i)),OMPPandf1Padyinc
           inc(i)=max(ihigh(i)-ic(i),ic(i)-ilow(i))+OMPPandf1Padyinc
        enddo
        do iv=1,neq
            do i=1,Nthreads
                if (igyl(iv,2).le.ihigh(i) .and. igyl(iv,2).ge.ilow(i)) then
                    ivthread(iv)=i
                endif
            enddo
        enddo

end subroutine OMPSplitIndexyPandf

!subroutine ParallelPandf1x(neq,time,yl,yldot)
!    use Output
!    use omp_lib
!    use OMPSettings,only: Nthreads
!    use OMPPandf1Settings,only: NthreadsPandf1,OMPPandf1FlagVerbose,OMPCheckParallelPandf1,&
!    OMPTimeParallelPandf1,OMPTimeSerialPandf1,OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
!    use OMPPandf1,only:OMPic,OMPlinc,OMPrinc,OMPivthread,OMPTimeCollectPandf1,OMPTimeLocalPandf1
!    use Dim,only:nx
!    use Selec, only: xlinc,xrinc
!    integer xlinc_bkp,xrinc_bkp,iv,tid
!    integer,intent(in)::neq
!    real,intent(in)::yl(*)
!    real,intent(out)::yldot(*)
!    real,intent(in)::time
!    real::yldotcopy(1:neq)
!    real yldotsave(1:neq),ylcopy(1:neq+2)
!     character*80 ::FileName
!
!    integer::ith
!    ylcopy(1:neq+2)=yl(1:neq+2)
!    NthreadsPandf1=min(Nthreads,nx)
!    !NthreadsPandf1=nx+2
!    if (OMPPandf1Verbose.gt.0.and. OMPPandf1FlagVerbose.eq.1) then
!     write(*,'(a,a,i3,a,i3)') OMPPandf1Stamp,' Number of threads for omp calculations:',NthreadsPandf1,'for nx=',nx
!     OMPPandf1FlagVerbose=0
!    endif
!    call OMPSplitIndexPandf(0,nx+1,NthreadsPandf1,OMPic,OMPlinc,OMPrinc,OMPivthread,neq)
!    if (OMPPandf1Verbose.gt.0) then
!        write(*,'(a,a,I3)') OMPPandf1Stamp,' nx=',nx
!        write(*,*)OMPPandf1Stamp, ' ic(ith), linc(ith), rinc(ith)'
!        do ith=1,NthreadsPandf1
!            write(*,'(a,I3,a,I7,I7,I7)') '  *    ithread ', ith,':',OMPic(ith),OMPlinc(ith),OMPrinc(ith)
!        enddo
!    endif
!    xlinc=1
!    xrinc=1
!    xlinc_bkp=xlinc
!    xrinc_bkp=xrinc
!
!    OMPTimeParallelPandf1=omp_get_wtime()
! write (*,*) 'v115'
! OMPTimeLocalPandf1(1)=omp_get_wtime()
!    !call convsr_vo (-1, -1, yl)
!
!    !call convsr_aux (-1, -1)
!
!OMPTimeLocalPandf1(1)=omp_get_wtime()-OMPTimeLocalPandf1(1)
!write(*,*) ' Time in convert:',OMPTimeLocalPandf1(1)
!   !!$omp parallel do default(shared)&
!    !!$omp& firstprivate(yldotcopy,ylcopy)&
!    !!$omp& private(iv,tid)
!
!
!
!    loopthread: do ith=1,NthreadsPandf1 !ith from 1 to Nthread, tid from 0 to Nthread-1
!        tid=omp_get_thread_num()
!        if (OMPPandf1Debug.gt.0) write(*,*) OMPPandf1Stamp,'Thread id:',tid,' <-> ith:',ith
!        ! we keep all these parameters as it is easier to debug LocalJacBuilder and deal with private/shared attributes
!        !xlinc=OMPlinc(ith)
!        !xrinc=OMPrinc(ith)
!        OMPTimeLocalPandf1(ith)=omp_get_wtime()
!        call pandf1 (OMPic(ith), -1, 0, neq, time, ylcopy, yldotcopy)
!        OMPTimeLocalPandf1(ith)=omp_get_wtime()-OMPTimeLocalPandf1(ith)
!        !write(FileName,'(a,I2.2)') "dump_parapandf_",ith
!        !call DebugHelper(trim(FileName))
!        OMPTimeCollectPandf1(ith)=omp_get_wtime()
!        !!!$omp  critical
!
!        do iv=1,neq
!        if (OMPivthread(iv)==ith) then
!
!        yldot(iv)=yldotcopy(iv)
!
!        endif
!        !!!$omp  end critical
!        enddo
!        OMPTimeCollectPandf1(ith)=omp_get_wtime()-OMPTimeCollectPandf1(ith)
!
!    if (OMPPandf1Verbose.gt.1) then
!            write(*,*) OMPPandf1Stamp,' Time in thread #', tid,':',OMPTimeLocalPandf1(ith),OMPTimeCollectPandf1(ith)
!    endif
!    OMPTimeLocalPandf1(ith)=OMPTimeCollectPandf1(ith)+OMPTimeLocalPandf1(ith)
!    enddo loopthread
!!!$omp  END PARALLEL DO
!
!    xlinc=xlinc_bkp
!    xrinc=xrinc_bkp
!
!      OMPTimeParallelPandf1=omp_get_wtime()-OMPTimeParallelPandf1
!      if (OMPCheckParallelPandf1.gt.0) then
!          OMPTimeSerialPandf1=omp_get_wtime()
!call pandf1 (-1, -1, 0, neq, time, ylcopy, yldotsave)
!!write(FileName,'(a)') "dump_serialpandf"
!        !call DebugHelper(trim(FileName))
!          OMPTimeSerialPandf1=omp_get_wtime()-OMPTimeSerialPandf1
!          write(*,*) "Timing Pandf1 serial:",OMPTimeSerialPandf1,"/parallel:",OMPTimeParallelPandf1
!          call Compare(yldot,yldotsave,neq)
!          write(*,*) "serial and parallel pandf are identical"
!        endif
!
!
!end subroutine ParallelPandf1x
#ifdef _OPENMP
subroutine ParallelPandf1(neq,time,yl,yldot)

    use omp_lib
    use OmpCopybbb
    use OMPSettings,only: Nthreads
    use OMPPandf1Settings,only: NthreadsPandf1,OMPPandf1FlagVerbose,OMPPandf1Check,&
    OMPTimeParallelPandf1,OMPTimeSerialPandf1,OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    use OMPPandf1,only:OMPic,OMPyinc,OMPivthread,OMPTimeCollectPandf1,OMPTimeLocalPandf1
    use OMPPandf1Settings,only:OMPPandf1RunPara,OMPPandf1FirstRun
    use Dim,only:ny
    use Selec, only:yinc
Use Grid,only:ijactot

    integer yinc_bkp,iv,tid,EffNthreadsPandf1
    integer,intent(in)::neq
    real,intent(in)::yl(*)
    real,intent(out)::yldot(*)
    real,intent(in)::time
    real::yldotcopy(1:neq)
    real yldotsave(1:neq),ylcopy(1:neq+2)
     character*80 ::FileName
    real time1,time2
    integer::ith
    !write(*,*)' ======================================== '
    ylcopy(1:neq+1)=yl(1:neq+1)




    if (ijactot.gt.0) then

     if (OMPPandf1FirstRun.eq.1) then

     EffNthreadsPandf1=min(Nthreads,min(NthreadsPandf1,ny))
     call gchange('OMPPandf1',0)

    call OMPSplitIndexyPandf(0,ny+1,EffNthreadsPandf1,OMPic,OMPyinc,OMPivthread,neq)
    if (OMPPandf1Verbose.gt.0) then
    write(*,'(a,a,I3,a,I3)') OMPPandf1Stamp,'EffNthreadsPandf1',EffNthreadsPandf1 ,' ny=',ny
endif
    if (OMPPandf1Verbose.gt.1) then
        write(*,*)OMPPandf1Stamp, ' ic(ith), linc(ith)'
        do ith=1,EffNthreadsPandf1
            write(*,'(a,I3,a,I7,I7)') '  *    ithread ', ith,':',OMPic(ith),OMPyinc(ith)
        enddo
    endif


     OMPPandf1FirstRun=0
    endif

    if (EffNthreadsPandf1.ne.min(Nthreads,min(NthreadsPandf1,ny))) then
    EffNthreadsPandf1=min(Nthreads,min(NthreadsPandf1,ny))
    call OMPSplitIndexyPandf(0,ny+1,EffNthreadsPandf1,OMPic,OMPyinc,OMPivthread,neq)
    endif
    if (OMPPandf1Verbose.gt.0.and. OMPPandf1FlagVerbose.eq.1) then
     write(*,'(a,a,i3,a,i3)') OMPPandf1Stamp,' Number of threads for omp calculations:',EffNthreadsPandf1,'for ny=',ny
     OMPPandf1FlagVerbose=0
    endif
    yinc_bkp=yinc
    Time1=omp_get_wtime()
    call OmpCopyPointerup
    !$omp parallel do default(shared) if(OMPPandf1RunPara.gt.0) &
    !$omp& private(iv,tid) firstprivate(ylcopy) private(yldotcopy)
    loopthread: do ith=1,EffNthreadsPandf1 !ith from 1 to Nthread, tid from 0 to Nthread-1
        tid=omp_get_thread_num()
        if (OMPPandf1Debug.gt.0) write(*,*) OMPPandf1Stamp,'Thread id:',tid,' <-> ith:',ith
        ! we keep all these parameters as it is easier to debug LocalJacBuilder and deal with private/shared attributes
        yinc=OMPyinc(ith)
        OMPTimeLocalPandf1(ith)=omp_get_wtime()

        call pandf1 (-1,OMPic(ith), 0, neq, time, ylcopy, yldotcopy)

        OMPTimeLocalPandf1(ith)=omp_get_wtime()-OMPTimeLocalPandf1(ith)
        OMPTimeCollectPandf1(ith)=omp_get_wtime()
        do iv=1,neq
        if (OMPivthread(iv)==ith) then
        yldot(iv)=yldotcopy(iv)
        endif
        enddo


        OMPTimeCollectPandf1(ith)=omp_get_wtime()-OMPTimeCollectPandf1(ith)
yinc=yinc_bkp
    if (OMPPandf1Verbose.gt.1) then
            write(*,*) OMPPandf1Stamp,' Time in thread #', tid,':',OMPTimeLocalPandf1(ith),OMPTimeCollectPandf1(ith)
    endif
    OMPTimeLocalPandf1(ith)=OMPTimeCollectPandf1(ith)+OMPTimeLocalPandf1(ith)
    enddo loopthread
!$omp  END PARALLEL DO
Time1=omp_get_wtime()-Time1

      OMPTimeParallelPandf1=Time1+OMPTimeParallelPandf1
        !do iv=1,neq
        !yldot(iv)=yldotcopy(iv,OMPivthread(iv))
        !enddo




      if (OMPPandf1Check.gt.0) then
          Time2=omp_get_wtime()
         call pandf1 (-1, -1, 0, neq, time, ylcopy, yldotsave)
         Time2=omp_get_wtime()-Time2
          OMPTimeSerialPandf1=Time2+OMPTimeSerialPandf1
          write(*,*) "Timing Pandf1 serial:",OMPTimeSerialPandf1,"(",Time2,")/parallel:",OMPTimeParallelPandf1,'(',Time1,')'
          write(*,*) 'ijactot ',ijactot
          call Compare(yldot,yldotsave,neq)
          write(*,*) "serial and parallel pandf are identical"
        endif
  else
       call pandf1 (-1,-1, 0, neq, time, yl, yldot)
endif
end subroutine ParallelPandf1

subroutine OmpCopyPointerbbb2
use OmpCopybbb
!call OmpCopyPointerctaue
!call OmpCopyPointerctaui
!call OmpCopyPointercfvcsx
!call OmpCopyPointercfvcsy
!call OmpCopyPointercfvgpx
!call OmpCopyPointercfvgpy
!call OmpCopyPointeriwalli
!call OmpCopyPointeriwallo
!call OmpCopyPointerfngysi
!call OmpCopyPointerfngyso
!call OmpCopyPointeryld_carbi
!call OmpCopyPointeryld_carbo
!call OmpCopyPointerfqpsatlb
!call OmpCopyPointerfqpsatrb
!call OmpCopyPointertotfeexl
!call OmpCopyPointertotfeexr
!call OmpCopyPointertotfeixl
!call OmpCopyPointertotfeixr
!call OmpCopyPointersputflxlb
!call OmpCopyPointersputflxrb
!call OmpCopyPointersputflxw
!call OmpCopyPointersputflxpf
call OmpCopyPointerparvis
call OmpCopyPointertravis
!call OmpCopyPointerdiffusivwrk
!call OmpCopyPointercoll_fe
!call OmpCopyPointercoll_fi
!call OmpCopyPointerni
!call OmpCopyPointernm
!call OmpCopyPointernz2
!call OmpCopyPointeruu
!call OmpCopyPointeruup


!call OmpCopyPointerup
!call OmpCopyPointerupi
!call OmpCopyPointeruz
!
!
!call OmpCopyPointerv2
!call OmpCopyPointerv2xgp
!call OmpCopyPointerv2ce
!call OmpCopyPointerv2cb
!call OmpCopyPointerve2cb
!call OmpCopyPointerv2cd
!call OmpCopyPointerve2cd
!call OmpCopyPointerq2cd
!call OmpCopyPointerv2rd
!call OmpCopyPointerv2dd
!call OmpCopyPointervy
!call OmpCopyPointervygp
!call OmpCopyPointervytan
!call OmpCopyPointervygtan
!call OmpCopyPointervyce
!call OmpCopyPointervycb
!call OmpCopyPointerveycb
!call OmpCopyPointervycp
!call OmpCopyPointerveycp
!call OmpCopyPointervyrd
!call OmpCopyPointervydd
!call OmpCopyPointervyavis
!call OmpCopyPointervex
!call OmpCopyPointerupe
!call OmpCopyPointervey
!
!!!!!!!!!!!!!!!!!!!!!!
!call OmpCopyPointervycf
!call OmpCopyPointervycr
!call OmpCopyPointerte
!call OmpCopyPointerti
!call OmpCopyPointerng
!call OmpCopyPointerlng
!call OmpCopyPointeruug
!call OmpCopyPointervyg
!call OmpCopyPointertg
!call OmpCopyPointertev
!call OmpCopyPointertiv
!call OmpCopyPointerniy0
!call OmpCopyPointerniy1
!call OmpCopyPointerniy0s
!call OmpCopyPointerniy1s
!call OmpCopyPointerney0
!call OmpCopyPointerney1
!call OmpCopyPointernity0
!call OmpCopyPointernity1
!call OmpCopyPointertey0
!call OmpCopyPointertey1
!call OmpCopyPointertiy0
!call OmpCopyPointertiy1
!call OmpCopyPointertiy0s
!call OmpCopyPointertiy1s
!call OmpCopyPointertgy0
!
!!call OmpCopyPointertgy1
!
!write(*,*) '#Master::: Calling routine to copy variable ngy0'
!
!call OmpCopyPointerngy0
!
!write(*,*) '#Master::: Calling routine to copy variable ngy1'
!
!call OmpCopyPointerngy1
!
!write(*,*) '#Master::: Calling routine to copy variable pgy0'
!
!call OmpCopyPointerpgy0
!
!write(*,*) '#Master::: Calling routine to copy variable pgy1'
!
!call OmpCopyPointerpgy1
!
!write(*,*) '#Master::: Calling routine to copy variable pg'
!
!call OmpCopyPointerpg
!
!write(*,*) '#Master::: Calling routine to copy variable phiy0'
!
!call OmpCopyPointerphiy0
!
!write(*,*) '#Master::: Calling routine to copy variable phiy1'
!
!call OmpCopyPointerphiy1
!
!write(*,*) '#Master::: Calling routine to copy variable phiy0s'
!
!call OmpCopyPointerphiy0s
!
!write(*,*) '#Master::: Calling routine to copy variable phiy1s'
!
!call OmpCopyPointerphiy1s
!
!write(*,*) '#Master::: Calling routine to copy variable pr'
!
!call OmpCopyPointerpr
!
!write(*,*) '#Master::: Calling routine to copy variable prev'
!
!call OmpCopyPointerprev
!
!write(*,*) '#Master::: Calling routine to copy variable prtv'
!
!call OmpCopyPointerprtv
!
!write(*,*) '#Master::: Calling routine to copy variable pri'
!
!call OmpCopyPointerpri
!
!write(*,*) '#Master::: Calling routine to copy variable priv'
!
!call OmpCopyPointerpriv
!
!write(*,*) '#Master::: Calling routine to copy variable priy0'
!
!call OmpCopyPointerpriy0
!
!write(*,*) '#Master::: Calling routine to copy variable priy1'
!
!call OmpCopyPointerpriy1
!
!write(*,*) '#Master::: Calling routine to copy variable pre'
!
!call OmpCopyPointerpre
!
!write(*,*) '#Master::: Calling routine to copy variable ne'
!
!call OmpCopyPointerne
!
!write(*,*) '#Master::: Calling routine to copy variable nit'
!
!call OmpCopyPointernit
!
!write(*,*) '#Master::: Calling routine to copy variable phi'
!
!call OmpCopyPointerphi
!
!write(*,*) '#Master::: Calling routine to copy variable phiv'
!
!call OmpCopyPointerphiv
!
!write(*,*) '#Master::: Calling routine to copy variable zeff'
!
!call OmpCopyPointerzeff
!
!write(*,*) '#Master::: Calling routine to copy variable loglambda'
!
!call OmpCopyPointerloglambda
!
!write(*,*) '#Master::: Calling routine to copy variable netap'
!
!call OmpCopyPointernetap
!call OmpCopyPointerznot
!call OmpCopyPointerupxpt
!call OmpCopyPointernixpt
!call OmpCopyPointervisyxpt
!call OmpCopyPointervyhxpt
!call OmpCopyPointervyvxpt
!call OmpCopyPointerfmihxpt
!call OmpCopyPointerfmivxpt
!call OmpCopyPointerrtaux
!call OmpCopyPointerrtauy
!call OmpCopyPointerrtau
!call OmpCopyPointerbetap
!call OmpCopyPointerfqp
!call OmpCopyPointerfq2
!call OmpCopyPointerfqx
!call OmpCopyPointerfqxb
!call OmpCopyPointerfdiaxlb
!call OmpCopyPointerfdiaxrb
!call OmpCopyPointerfloxebgt
!call OmpCopyPointerfloxibgt
!call OmpCopyPointerfqy
!
!write(*,*) '#Master::: Calling routine to copy variable fqyb'
!
!call OmpCopyPointerfqyb
!
!write(*,*) '#Master::: Calling routine to copy variable fqyn'
!
!call OmpCopyPointerfqyn
!
!write(*,*) '#Master::: Calling routine to copy variable fqym'
!
!call OmpCopyPointerfqym
!
!write(*,*) '#Master::: Calling routine to copy variable fqymi'
!
!call OmpCopyPointerfqymi
!
!write(*,*) '#Master::: Calling routine to copy variable fqya'
!
!call OmpCopyPointerfqya
!
!write(*,*) '#Master::: Calling routine to copy variable fqydt'
!
!call OmpCopyPointerfqydt
!
!write(*,*) '#Master::: Calling routine to copy variable fqydti'
!
!call OmpCopyPointerfqydti
!
!write(*,*) '#Master::: Calling routine to copy variable fqyao'
!
!call OmpCopyPointerfqyao
!
!write(*,*) '#Master::: Calling routine to copy variable fqyae'
!
!call OmpCopyPointerfqyae
!
!write(*,*) '#Master::: Calling routine to copy variable fqyai'
!
!call OmpCopyPointerfqyai
!
!write(*,*) '#Master::: Calling routine to copy variable fqyd'
!
!call OmpCopyPointerfqyd
!
!write(*,*) '#Master::: Calling routine to copy variable fqygp'
!
!call OmpCopyPointerfqygp
!
!write(*,*) '#Master::: Calling routine to copy variable fq2d'
!
!call OmpCopyPointerfq2d
!
!write(*,*) '#Master::: Calling routine to copy variable fnix'
!
!call OmpCopyPointerfnix
!call OmpCopyPointerfnixcb
!call OmpCopyPointerfniy
!call OmpCopyPointerfniy4ord
!call OmpCopyPointerfniycb
!call OmpCopyPointerfmix
!call OmpCopyPointerfmiy
!call OmpCopyPointerfmixy
!call OmpCopyPointerfmity
!call OmpCopyPointerfeex
!call OmpCopyPointerfeey
!call OmpCopyPointerfeexy
!call OmpCopyPointerfeey4ord
!call OmpCopyPointerfeix
!call OmpCopyPointerfeiy
!call OmpCopyPointerfegx
!call OmpCopyPointerfegy
!call OmpCopyPointerqipar
!call OmpCopyPointerfniycbo
!call OmpCopyPointerfeiycbo
!call OmpCopyPointerfeeycbo
!call OmpCopyPointerfeixy
!call OmpCopyPointerfeiy4ord
!call OmpCopyPointerfngx
!call OmpCopyPointerfngx4ord
!call OmpCopyPointerflngx
!call OmpCopyPointerfngy
!call OmpCopyPointerfngy4ord
!call OmpCopyPointerflngy
!call OmpCopyPointerfngxy
!call OmpCopyPointerflngxy
!call OmpCopyPointerbcel
!call OmpCopyPointerbcer
!call OmpCopyPointerbcil
!call OmpCopyPointerbcir
!call OmpCopyPointerkappal
!call OmpCopyPointerkappar
!call OmpCopyPointerdphi_iy1
!call OmpCopyPointerkincorlb
!call OmpCopyPointerkincorrb
!call OmpCopyPointerex
!call OmpCopyPointerey
!call OmpCopyPointergpix
!call OmpCopyPointergpiy
!call OmpCopyPointergpex
!call OmpCopyPointergpey
!call OmpCopyPointergprx
!call OmpCopyPointergpry
!call OmpCopyPointergtex
!call OmpCopyPointergtey
!!
!!write(*,*) '#Master::: Calling routine to copy variable gtix'
!!
!!call OmpCopyPointergtix
!!
!!write(*,*) '#Master::: Calling routine to copy variable gtiy'
!!
!!call OmpCopyPointergtiy
!!
!!write(*,*) '#Master::: Calling routine to copy variable frice'
!!
!!call OmpCopyPointerfrice
!!call OmpCopyPointerfrici
!!call OmpCopyPointerfricnrl
!!call OmpCopyPointeralfe
!!call OmpCopyPointerbetai
!!call OmpCopyPointerw
!!call OmpCopyPointerw0
!!call OmpCopyPointerw1
!!
!!write(*,*) '#Master::: Calling routine to copy variable w2'
!!
!!call OmpCopyPointerw2
!!
!!write(*,*) '#Master::: Calling routine to copy variable w3'
!!
!!call OmpCopyPointerw3
!!
!!write(*,*) '#Master::: Calling routine to copy variable wvh'
!!
!!call OmpCopyPointerwvh
!!
!!write(*,*) '#Master::: Calling routine to copy variable flox'
!!
!!call OmpCopyPointerflox
!!
!!write(*,*) '#Master::: Calling routine to copy variable floy'
!!
!!call OmpCopyPointerfloy
!!
!!write(*,*) '#Master::: Calling routine to copy variable conx'
!!
!!call OmpCopyPointerconx
!!
!!write(*,*) '#Master::: Calling routine to copy variable cony'
!!
!!call OmpCopyPointercony
!!
!!write(*,*) '#Master::: Calling routine to copy variable floxe'
!!
!!call OmpCopyPointerfloxe
!!
!!write(*,*) '#Master::: Calling routine to copy variable floye'
!!
!!call OmpCopyPointerfloye
!!
!!write(*,*) '#Master::: Calling routine to copy variable floxi'
!!
!!call OmpCopyPointerfloxi
!!
!!write(*,*) '#Master::: Calling routine to copy variable floyi'
!!
!!call OmpCopyPointerfloyi
!!
!!write(*,*) '#Master::: Calling routine to copy variable floxg'
!!
!!call OmpCopyPointerfloxg
!!
!!write(*,*) '#Master::: Calling routine to copy variable floyg'
!!
!!call OmpCopyPointerfloyg
!!
!!write(*,*) '#Master::: Calling routine to copy variable conxe'
!!
!!call OmpCopyPointerconxe
!!
!!write(*,*) '#Master::: Calling routine to copy variable conye'
!!
!!call OmpCopyPointerconye
!!
!!write(*,*) '#Master::: Calling routine to copy variable conxi'
!!
!!call OmpCopyPointerconxi
!!
!!write(*,*) '#Master::: Calling routine to copy variable conyi'
!!
!!call OmpCopyPointerconyi
!!
!!write(*,*) '#Master::: Calling routine to copy variable conxg'
!!
!!call OmpCopyPointerconxg
!!
!!write(*,*) '#Master::: Calling routine to copy variable conyg'
!!
!!call OmpCopyPointerconyg
!!
!!write(*,*) '#Master::: Calling routine to copy variable floxge'
!!
!!call OmpCopyPointerfloxge
!!
!!write(*,*) '#Master::: Calling routine to copy variable floyge'
!!
!!call OmpCopyPointerfloyge
!!
!!write(*,*) '#Master::: Calling routine to copy variable conxge'
!!
!!call OmpCopyPointerconxge
!!
!!write(*,*) '#Master::: Calling routine to copy variable conyge'
!!
!!call OmpCopyPointerconyge
!!
!!write(*,*) '#Master::: Calling routine to copy variable visx'
!!
!!call OmpCopyPointervisx
!!
!!write(*,*) '#Master::: Calling routine to copy variable visy'
!!
!!call OmpCopyPointervisy
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxe'
!!
!!call OmpCopyPointerhcxe
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcye'
!!
!!call OmpCopyPointerhcye
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxij'
!!
!!call OmpCopyPointerhcxij
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcyij'
!!
!!call OmpCopyPointerhcyij
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxg'
!!
!!call OmpCopyPointerhcxg
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcyg'
!!
!!call OmpCopyPointerhcyg
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxi'
!!
!!call OmpCopyPointerhcxi
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxineo'
!!
!!call OmpCopyPointerhcxineo
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcyi'
!!
!!call OmpCopyPointerhcyi
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxn'
!!
!!call OmpCopyPointerhcxn
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcyn'
!!
!!call OmpCopyPointerhcyn
!!
!!write(*,*) '#Master::: Calling routine to copy variable kxbohm'
!!
!!call OmpCopyPointerkxbohm
!!
!!write(*,*) '#Master::: Calling routine to copy variable kybohm'
!!
!!call OmpCopyPointerkybohm
!!
!!write(*,*) '#Master::: Calling routine to copy variable dif_use'
!!
!!call OmpCopyPointerdif_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable difp_use'
!!
!!call OmpCopyPointerdifp_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable dif2_use'
!!
!!call OmpCopyPointerdif2_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable tray_use'
!!
!!call OmpCopyPointertray_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable trax_use'
!!
!!call OmpCopyPointertrax_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable kye_use'
!!
!!call OmpCopyPointerkye_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable kyi_use'
!!
!!call OmpCopyPointerkyi_use
!!call OmpCopyPointerkxe_use
!!call OmpCopyPointerkxi_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable dutm_use'
!!
!!call OmpCopyPointerdutm_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable vy_use'
!!
!!call OmpCopyPointervy_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable vy_cft'
!!
!!call OmpCopyPointervy_cft
!!
!!write(*,*) '#Master::: Calling routine to copy variable vyte_cft'
!!
!!call OmpCopyPointervyte_cft
!!
!!write(*,*) '#Master::: Calling routine to copy variable vyti_cft'
!!
!!call OmpCopyPointervyti_cft
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuiz'
!!
!!call OmpCopyPointernuiz
!!
!!write(*,*) '#Master::: Calling routine to copy variable nucx'
!!
!!call OmpCopyPointernucx
!!
!!write(*,*) '#Master::: Calling routine to copy variable nucxi'
!!
!!call OmpCopyPointernucxi
!!
!!write(*,*) '#Master::: Calling routine to copy variable nueli'
!!
!!call OmpCopyPointernueli
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuelg'
!!
!!call OmpCopyPointernuelg
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuix'
!!
!!call OmpCopyPointernuix
!!
!!write(*,*) '#Master::: Calling routine to copy variable nurc'
!!
!!call OmpCopyPointernurc
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuvl'
!!
!!call OmpCopyPointernuvl
!!
!!write(*,*) '#Master::: Calling routine to copy variable eqp'
!!
!!call OmpCopyPointereqp
!!
!!write(*,*) '#Master::: Calling routine to copy variable eqpg'
!!
!!call OmpCopyPointereqpg
!!
!!write(*,*) '#Master::: Calling routine to copy variable eeli'
!!
!!call OmpCopyPointereeli
!!
!!write(*,*) '#Master::: Calling routine to copy variable pradhyd'
!!
!!call OmpCopyPointerpradhyd
!!
!!write(*,*) '#Master::: Calling routine to copy variable eta1'
!!
!!call OmpCopyPointereta1
!!
!!write(*,*) '#Master::: Calling routine to copy variable rtaue'
!!
!!call OmpCopyPointerrtaue
!!
!!write(*,*) '#Master::: Calling routine to copy variable dclass_e'
!!
!!call OmpCopyPointerdclass_e
!!
!!write(*,*) '#Master::: Calling routine to copy variable dclass_i'
!!
!!call OmpCopyPointerdclass_i
!!
!!write(*,*) '#Master::: Calling routine to copy variable visxneo'
!!
!!call OmpCopyPointervisxneo
!!
!!write(*,*) '#Master::: Calling routine to copy variable visvol_v'
!!
!!call OmpCopyPointervisvol_v
!!
!!write(*,*) '#Master::: Calling routine to copy variable visvol_q'
!!
!!call OmpCopyPointervisvol_q
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuii'
!!
!!call OmpCopyPointernuii
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuiistar'
!!
!!call OmpCopyPointernuiistar
!!
!!write(*,*) '#Master::: Calling routine to copy variable alfneo'
!!
!!call OmpCopyPointeralfneo
!!
!!write(*,*) '#Master::: Calling routine to copy variable k2neo'
!!
!!call OmpCopyPointerk2neo
!!
!!write(*,*) '#Master::: Calling routine to copy variable ktneo'
!!
!!call OmpCopyPointerktneo
!!
!!write(*,*) '#Master::: Calling routine to copy variable snic'
!!
!!call OmpCopyPointersnic
!!
!!write(*,*) '#Master::: Calling routine to copy variable sniv'
!!
!!call OmpCopyPointersniv
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorc'
!!
!!call OmpCopyPointerpsorc
!!call OmpCopyPointerpsor
!!call OmpCopyPointerpsorxrc
!!call OmpCopyPointerpsorxr
!!
!!write(*,*) '#Master::: Calling routine to copy variable psor_tmpov'
!!
!!call OmpCopyPointerpsor_tmpov
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorgc'
!!
!!call OmpCopyPointerpsorgc
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorg'
!!
!!call OmpCopyPointerpsorg
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorrgc'
!!
!!call OmpCopyPointerpsorrgc
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorrg'
!!
!!call OmpCopyPointerpsorrg
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorcxgc'
!!
!!call OmpCopyPointerpsorcxgc
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorcxg'
!!
!!call OmpCopyPointerpsorcxg
!!
!!write(*,*) '#Master::: Calling routine to copy variable psori'
!!
!!call OmpCopyPointerpsori
!!
!!write(*,*) '#Master::: Calling routine to copy variable psordis'
!!
!!call OmpCopyPointerpsordis
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorbgg'
!!
!!call OmpCopyPointerpsorbgg
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorbgz'
!!
!!call OmpCopyPointerpsorbgz
!!
!!write(*,*) '#Master::: Calling routine to copy variable erliz'
!!
!!call OmpCopyPointererliz
!!
!!write(*,*) '#Master::: Calling routine to copy variable erlrc'
!!
!!call OmpCopyPointererlrc
!!
!!write(*,*) '#Master::: Calling routine to copy variable vsoreec'
!!
!!call OmpCopyPointervsoreec
!!
!!write(*,*) '#Master::: Calling routine to copy variable vsoree'
!!
!!call OmpCopyPointervsoree
!!
!!write(*,*) '#Master::: Calling routine to copy variable pwrebkg'
!!
!!call OmpCopyPointerpwrebkg
!!
!!write(*,*) '#Master::: Calling routine to copy variable pwribkg'
!!
!!call OmpCopyPointerpwribkg
!!
!!write(*,*) '#Master::: Calling routine to copy variable wjdote'
!!
!!call OmpCopyPointerwjdote
!!
!!write(*,*) '#Master::: Calling routine to copy variable smoc'
!!
!!call OmpCopyPointersmoc
!!
!!write(*,*) '#Master::: Calling routine to copy variable smov'
!!
!!call OmpCopyPointersmov
!!
!!write(*,*) '#Master::: Calling routine to copy variable msor'
!!
!!call OmpCopyPointermsor
!!
!!write(*,*) '#Master::: Calling routine to copy variable msorxr'
!!
!!call OmpCopyPointermsorxr
!!
!!write(*,*) '#Master::: Calling routine to copy variable seec'
!!
!!call OmpCopyPointerseec
!!
!!write(*,*) '#Master::: Calling routine to copy variable seev'
!!
!!call OmpCopyPointerseev
!!
!!write(*,*) '#Master::: Calling routine to copy variable seic'
!!
!!call OmpCopyPointerseic
!!
!!write(*,*) '#Master::: Calling routine to copy variable seiv'
!!
!!call OmpCopyPointerseiv
!!
!!
!!
!!call OmpCopyPointerresco
!!
!! to copy variable resng'
!!
!!call OmpCopyPointerresng
!!
!!write(*,*) '#Master::: Calling routine to copy variable reseg'
!!
!!call OmpCopyPointerreseg
!!
!!write(*,*) '#Master::: Calling routine to copy variable resmo'
!!
!!call OmpCopyPointerresmo
!!
!!write(*,*) '#Master::: Calling routine to copy variable resee'
!!
!!call OmpCopyPointerresee
!!
!!write(*,*) '#Master::: Calling routine to copy variable resei'
!!
!!call OmpCopyPointerresei
!!
!!write(*,*) '#Master::: Calling routine to copy variable resphi'
!!
!!call OmpCopyPointerresphi
!!
!!write(*,*) '#Master::: Calling routine to copy variable sng_ue'
!!
!!call OmpCopyPointersng_ue
!!
!!write(*,*) '#Master::: Calling routine to copy variable seg_ue'
!!
!!call OmpCopyPointerseg_ue
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorold'
!!
!!call OmpCopyPointerpsorold
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorxrold'
!!
!!call OmpCopyPointerpsorxrold
!!
!!write(*,*) '#Master::: Calling routine to copy variable msorold'
!!
!!call OmpCopyPointermsorold
!!
!!write(*,*) '#Master::: Calling routine to copy variable msorxrold'
!!
!!call OmpCopyPointermsorxrold
!!
!!write(*,*) '#Master::: Calling routine to copy variable yldot_pert'
!!
!!call OmpCopyPointeryldot_pert
!!call OmpCopyPointernzloc
!!call OmpCopyPointerimpradloc
!!call OmpCopyPointerpwrzec
!!call OmpCopyPointerpwrze
!!call OmpCopyPointerpradc
!!call OmpCopyPointerpradcff
!!call OmpCopyPointerprad
!!call OmpCopyPointerpradzc
!!call OmpCopyPointerpradz
!!call OmpCopyPointerna
!!call OmpCopyPointerntau
!!call OmpCopyPointernratio
!!call OmpCopyPointertempa
!!call OmpCopyPointerden
!!call OmpCopyPointergradp
!!call OmpCopyPointergradt
!!call OmpCopyPointerdztot
end subroutine OmpCopyPointerbbb2


subroutine ompCopyConvert
      use OmpCopybbb
      implicit none
      call OmpCopyPointerne
      call OmpCopyPointernit
      call OmpCopyPointernm
      call OmpCopyPointernz2
      call OmpCopyPointerni
      call OmpCopyPointerng
      call OmpCopyPointerlng
      call OmpCopyPointerte
      call OmpCopyPointertg
      call OmpCopyPointerti
      call OmpCopyPointerphi
      call OmpCopyPointerup
      call OmpCopyPointerpr
      call OmpCopyPointerzeff
      call OmpCopyPointerpri
      call OmpCopyPointerpre
      call OmpCopyPointerznot
      call OmpCopyPointertg
      call OmpCopyPointerpg
      call OmpCopyPointergprx
      call OmpCopyPointerney0
      call OmpCopyPointerney1
      call OmpCopyPointernity0
      call OmpCopyPointernity1
      call OmpCopyPointergpry
      call OmpCopyPointergpix
      call OmpCopyPointerniy0
      call OmpCopyPointerniy1
      call OmpCopyPointerniy0s
      call OmpCopyPointerniy1s
      call OmpCopyPointerpriy0
      call OmpCopyPointerpriy1
      call OmpCopyPointergpiy
      call OmpCopyPointertey0
      call OmpCopyPointertey1
      call OmpCopyPointertiy0
      call OmpCopyPointertiy1
      call OmpCopyPointerphiy0
      call OmpCopyPointerphiy1
      call OmpCopyPointertiy0s
      call OmpCopyPointertiy1s
      call OmpCopyPointerphiy0s
      call OmpCopyPointerphiy1s
      call OmpCopyPointerngy0
      call OmpCopyPointerngy1
      call OmpCopyPointertgy0
      call OmpCopyPointertgy1
      call OmpCopyPointerpgy0
      call OmpCopyPointerpgy1
      call OmpCopyPointergpex
      call OmpCopyPointergtex
      call OmpCopyPointergtix
      call OmpCopyPointerex
      call OmpCopyPointergpey
      call OmpCopyPointergtey
      call OmpCopyPointergtiy
      call OmpCopyPointerey
      call OmpCopyPointerphiv
      call OmpCopyPointertiv
      call OmpCopyPointertev
      call OmpCopyPointerprev
      call OmpCopyPointerprtv
      call OmpCopyPointerpriv
      end subroutine ompCopyConvert
#endif
