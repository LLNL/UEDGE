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

#endif


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
