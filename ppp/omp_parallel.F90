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
#ifndef _OPENMP
     call xerrab("UEDGE was not compiled with OpenMP. Cannot use OMP parallel features.")
#else
    ! prepare MPI/Hybrid stamps for output
    if (HybridOMPMPI>0) then
        write(MPIRankTag,'(I4)') MPIrank
        write(OMPJacStamp,'(a,a,a)') '[',trim(adjustl(trim(MPIRankTag))),'] OMPJac* '
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
#endif
end subroutine InitOMP
!-------------------------------------------------------------------------------------------------
#ifdef _OPENMP
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


!!!! >>>>>> Routine for pandf1 parallelization <<<<<<<<<<<<<
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

#ifdef _OPENMP
subroutine OMPPandf1Rhs(neq,time,yl,yldot)

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
end subroutine OMPPandf1Rhs
#endif
