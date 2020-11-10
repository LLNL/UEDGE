#ifdef MPIJAC
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
