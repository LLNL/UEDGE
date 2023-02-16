#include "../bbb/bbb.h"
#include "../com/com.h"
!#include "../mppl.h"
!#include "../sptodp.h"
!-------------------------------------------------------------------------------------------------
subroutine InitParallel

    Use ParallelSettings,only:OMPParallelJac,MPIParallelJac,OMPParallelPandf1,ParallelWARNING
    Use ParallelEval,only:ParallelJac,ParallelPandf1
    Use HybridSettings,only:HybridOMPMPI
    implicit none


    if ((OMPParallelJac==1 .or. OMPParallelPandf1.gt.0) .and. MPIParallelJac==0) then
    if (ParallelWARNING.gt.0) then
    write(*,*) "<<< WARNING >>> : You are using openmp routines to evaluate the Jacobian or/and pandf1."
    write(*,*) "<<< WARNING >>> : Read first the documentation about how to use these options."
    write(*,*) "<<< WARNING >>> : If this is your first run with this configuration, you should verify whether&
    the parallel evaluation of pandf1 and jacobian are correct by setting ppp.CheckJac=1 and ppp.OMPCheckPandf1=1.&
     If not error is thrown, turn these settings off (=0)."
    write(*,*) "<<< WARNING >>> : You can turn off this warning by setting ppp.ParallelWarning=0."

    endif
        call InitOMP
        if (OMPParallelJac.gt.0) ParallelJac=1
        if (OMPParallelPandf1.gt.0) ParallelPandf1=1
        if (OMPParallelJac==0 .and. OMPParallelPandf1.gt.0) then
        call xerrab('Cannot run omp parallel evaluation of pandf1 without running jacobian omp evaluation.')
        endif
    elseif (OMPParallelJac==0 .and. MPIParallelJac==1) then
        write(*,'(a)') 'Parallelization of Jacobian evaluation with MPI not available yet'
        call InitMPI
        ParallelJac=1
    elseif (OMPParallelJac==1 .and. MPIParallelJac==1) then
        write(*,'(a)') 'Parallelization of Jacobian evaluation with MPI not available yet'
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
    integer i,iv
    real t_start
    real,external::tick,tock
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


      write(*,*) '--- Checking whether parallel and serial evaluations of jacobian are the same...'
       write(*,*) '----- Performing serial Evaluation of jacobian...'
      t_start=tick()
      call jac_calc (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jaccopy, jacopy, iacopy)
      write(*,*) '----- Serial Evaluation of jacobian done in: ',tock(t_start) ,' s'
      write(*,*) '----- Comparing jacobians (nzmax=',ia(neq)-1,')'
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

           write(*,*) '----- Writing jacobian into serialjac.dat and paralleljac.dat'
           call jac_write('paralleljac.dat',neq, jac, ja, ia)
           call jac_write('serialjac.dat',neq, jaccopy, jacopy, iacopy)

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
      write(*,*) '--- Parallel and serial evaluations of jacobian are identical...'

      endif

end subroutine jac_calc_parallel
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
    use Jacaux,only:ExtendedJacPhi
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
        if (ExtendedJacPhi.gt.0 .and. isphion*isnewpot.eq.1 .and. mod(iv,numvar).eq.0) then
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

!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------

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
integer::istop
real yldot(neq+2),yldotsave(neq+2)
istop=0
            do iv=1,neq
                if (abs(yldotsave(iv)-yldot(iv)).gt.1e-14) then
                    if (max(abs(yldot(iv)),abs(yldotsave(iv)))>0) then
                    if (abs(yldotsave(iv)-yldot(iv))/max(abs(yldot(iv)),abs(yldotsave(iv)))>1e-14) then
                    write(*,*) '>>>>',iv,yldotsave(iv),yldot(iv),abs(yldotsave(iv)-yldot(iv))/max(abs(yldot(iv)),abs(yldotsave(iv)))
                    istop=1
                    endif
                    endif

                endif
            enddo
           if (istop.gt.0) call xerrab('diff in rhsnk')


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

        real function tick()
        implicit none
            integer :: now, clock_rate
            call system_clock(now,clock_rate)
            tick=real(now)/real(clock_rate)
        end function tick

        real function tock(t)
         implicit none
            real, intent(in) :: t
            integer :: now, clock_rate
            call system_clock(now,clock_rate)

            tock = real(now)/real(clock_rate)-t
        end function tock






