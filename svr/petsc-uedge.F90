!
! $Id: petsc-uedge.F90,v 7.2 2019/10/14 22:59:29 meyer8 Exp $
!

!  Module used within a single call of PetscSnes()
! ------------------------------------------------
      module Petscmod
        implicit none
#include "finclude/petsc.h90"
!!#include "finclude/petscvec.h"
!!#include "finclude/petscvec.h90"
!!#include "finclude/petscmat.h"
!!#include "finclude/petscmat.h90"
!!#include "finclude/petscviewer.h"
!!#include "finclude/petscksp.h"
!!#include "finclude/petscpc.h"
!!#include "finclude/petscsnes.h"
!!#include "finclude/petscis.h"

      type Petscctx
        PetscScalar,pointer :: savfpoint(:),ylpoint(:)
        PetscScalar,pointer :: wkpoint(:)
        PetscInt,pointer :: icnstrpoint(:)
        PetscReal :: rlx
        PetscInt :: nnzmx    !maximum nom of nonzeros in Jacobian matrix, see Use(Jacobian)
        PetscInt :: incpset  !maximum num of Jacobian reevaluation, see Use(Lsode)
      end type Petscctx
      PetscLogStage :: stages(5)
      end module Petscmod

! ---------------------------------------------------------------------
      module Petscmodinterfaces
        use Petscmod
        implicit none

      Interface SNESSetApplicationContext
        Subroutine SNESSetApplicationContext(snes,ctx,ierr)
        use Petscmod
          SNES ::           snes
          type(Petscctx) :: ctx
          PetscErrorCode :: ierr
        End Subroutine
      End Interface SNESSetApplicationContext

      Interface SNESGetApplicationContext
        Subroutine SNESGetApplicationContext(snes,ctx,ierr)
        use Petscmod
          SNES ::           snes
          type(Petscctx), pointer :: ctx
          PetscErrorCode :: ierr
        End Subroutine
      End Interface SNESGetApplicationContext

      end module Petscmodinterfaces

! ---------------------------------------------------------------------
      module Petscmodinterface2
        use Petscmod
        implicit none

      Interface PetscSnesCreateObjs
        Subroutine PetscSnesCreateObjs(comm,n,neqg,nnzmx,ierr)
        use Petscmod
          integer4 ::       comm
          PetscInt ::       n,neqg,nnzmx
          PetscErrorCode :: ierr
        End Subroutine
      End Interface PetscSnesCreateObjs

      Interface PetscSnesApply
        Subroutine PetscSnesApply(n,yl,savf,suscal,sfscal,nnzmx,ireorder,incpset,icflag,rlx,icnstr,dtreal,iterm,ierr)
        use Petscmod
        PetscInt ::           n,nnzmx,ireorder,incpset
        PetscScalar,target :: yl(n),savf(n),suscal(n),sfscal(n)
        PetscInt,target ::    icnstr(n)
        PetscInt ::           icflag
        PetscReal ::          rlx,dtreal
        PetscInt ::           iterm
        PetscErrorCode ::     ierr
        End Subroutine
      End Interface PetscSnesApply
      end module Petscmodinterface2

! ---------------------------------------------------------------------
      module Snesmonitormod
      use Petscmod
      implicit none
      save
      type monctx
        PetscInt ::      buildJacobian,JacWasJustRebuilt
        PetscReal ::     thisdt
      end type monctx
      end module Snesmonitormod

! ---------------------------------------------------------------------
      subroutine PetscInitWrap()
      use Petscmod
      implicit none

      PetscErrorCode :: ierr
      PetscTruth     :: isInit
      call PetscInitialized(isInit, ierr)
      if (.not. isInit) then
        call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
        call PetscLogStageRegister("PetscSnesCreateObjs",stages(1),ierr)
        call PetscLogStageRegister("PetscSnesSetup",stages(2),ierr)
        call PetscLogStageRegister("PetscSnesApply",stages(3),ierr)
        call PetscLogStageRegister("PetscSnesFinalize",stages(4),ierr)
      end if

      end subroutine PetscInitWrap
! ---------------------------------------------------------------------
      subroutine PetscFinalizeWrap()
      use PETSc_Snes_Param
      implicit none

      PetscErrorCode :: ierr

!     Destroy global objects
      if (psp_F .ne. 0) then
        call VecDestroy(psp_F,ierr)
      endif
      if (psp_X .ne. 0) then
        call VecDestroy(psp_X,ierr)
      endif
      if (psp_Xtmp .ne. 0) then
        call VecDestroy(psp_Xtmp,ierr)
      endif
      if (psp_Jmat .ne. 0) then
        call MatDestroy(psp_Jmat,ierr)
      endif
      if (psp_Jmf .ne. 0) then
        call MatDestroy(psp_Jmf,ierr)
      endif
      if (psp_ISUVEC .ne. 0) then
        call VecDestroy(psp_ISUVEC,ierr)
      endif
      if (psp_SUVEC .ne. 0) then
        call VecDestroy(psp_SUVEC,ierr)
      endif
      if (psp_SFVEC .ne. 0) then
        call VecDestroy(psp_SFVEC,ierr)
      endif
      if (psp_usecolor .ne. 0) then
        !write(6,*) ' MatFDColoringDestroy'
        !call flush(6)
        call MatFDColoringDestroy(psp_fdcoloring,ierr)
      end if
      if (psp_snes .ne. 0) then
        call SNESDestroy(psp_snes,ierr)
      endif

      call PetscFinalize(ierr)
!      write(6,*) ' PetscFinalizeWrap'
!      call flush(6)
      end subroutine PetscFinalizeWrap

! ---------------------------------------------------------------------
      subroutine PetscInsertOpts()
!      use MpiVars
      use Npes_mpi
      use UEpar
      PetscErrorCode ::  ierr
      PetscInt ::           uedgeComm
      PetscTruth PETSC_TRUE,PETSC_FALSE
      parameter (PETSC_TRUE = .true.,PETSC_FALSE = .false.)


!     insert options from file
      if (petscoptfile .ne. '') then
        call PetscOptionsInsertFile(uedgeComm,petscoptfile,PETSC_TRUE,ierr)
      endif
!     insert options from the specified options string
      if (petscopts .ne. '') then
        call PetscOptionsInsertString(petscopts,ierr)
      endif
      end subroutine PetscInsertOpts

! ---------------------------------------------------------------------
      subroutine RhsCmp_Debug(comm,cmp_rhs,snes,pctx,n,yl)

      use Petscmod
      use PETSc_Snes_Param
      use Indices_loc_glob_map !ivloc2sdg,ivloc2mdg,ivl2gstnl - rm later!
      use Indices_domain_dcl   !ivloc2sdgl,ivloc2mdgl

      implicit none

      integer4 ::           comm
      PetscReal ::          cmp_rhs
      SNES ::               snes
      type(Petscctx) ::     pctx
      PetscScalar,target :: yl(n)

      PetscErrorCode ::     ierr
      integer4 ::           mype,npes
      Vec                   Vseq,Vmpi,Xseq,Xmpi
      VecScatter            scatter
      PetscScalar, pointer :: array(:),array_seq(:)
      PetscInt ::           id,i,j,n,neqg
      PetscScalar ::        error
      PetscReal ::          delta,norm
      PetscViewer ::        viewer
      PetscScalar,target :: yldot(n),yl_tmp(n)
      PetscTruth ::         flg
      PetscInt ::           uedgeComm

      external FormUedgeFunction

      call MPI_Comm_size(comm,npes,ierr)
      call MPI_Comm_rank(comm,mype,ierr)
!      write(6,*)" RhsCmp_Debug() npes,mype: ",npes,mype
!      call flush(6)

      call VecNorm(psp_X,NORM_2,norm,ierr)
      if (mype .eq. 0) then
        write(6,*)" initial Xnorm: ",norm
        call flush(6)
      end if
      call VecGetSize(psp_X,neqg,ierr)
      call VecGetLocalSize(psp_X,n,ierr)

      psp_isf = 0
      delta   = 0.0
      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-cmp_Frhs_delta',delta,flg,ierr)
      delta = 1.0 + delta
      if (flg .and. (mype .eq. 0)) then
        write(6,*)" scale X with",delta
        call flush(6)
      endif
      call VecScale(psp_X,delta,ierr)
      call FormUedgeFunction(snes,psp_X,psp_F,pctx,ierr)
      !write(6,*)" psp_F:"
      !call flush(6)
      !call VecView(psp_F,PETSC_VIEWER_STDOUT_WORLD,ierr)
      !write(6,*)" psp_SUVEC:"
      !call flush(6)
      !call VecView(psp_SUVEC,PETSC_VIEWER_STDOUT_WORLD,ierr)

      if (npes .eq. 1) then
        write(6,*)"write X in binary to datafile/X ..."
        call flush(6)
        call PetscViewerBinaryOpen(comm,"datafile/X",FILE_MODE_WRITE,viewer,ierr)
        call VecView(psp_X,viewer,ierr)
        call PetscViewerDestroy(viewer,ierr)

        write(6,*)"write F in binary to datafile/Frhs ..."
        write(6,*)" "
        call flush(6)
        call PetscViewerBinaryOpen(comm,"datafile/Frhs",FILE_MODE_WRITE,viewer,ierr)
        call VecView(psp_F,viewer,ierr)
        call PetscViewerDestroy(viewer,ierr)
      else
        ! compare parallel F with the sequential F
#if defined(DEBUG_rhs_1)
!     Only compare values in proc[0]
      if (mype .eq. 0) then
        write(6,*)"[",mype,"] load F from datafile/Frhs ..."
        write(6,*)" "
        call flush(6)
        call PetscViewerBinaryOpen(MPI_COMM_SELF,"datafile/Frhs",FILE_MODE_READ,viewer,ierr)
        call VecLoad(viewer, VECSEQ, Vseq,ierr)
        call PetscViewerDestroy(viewer,ierr)

        call VecGetArrayF90(psp_F,array,ierr)
        call VecGetArrayF90(Vseq,array_seq,ierr)
        if (mype .eq. 0) then
          write(6,*)"[",mype,"] Compare rhs function  F_mpi and F_seq.........."
          call flush(6)
        endif
        do i=1,n
          error = array(i) - array_seq(ivloc2sdgl(i))
          !error = (array(i) - array_seq(ivloc2sdgl(i)))/(array_seq(ivloc2sdgl(i)))
          if (error .lt. 0.0) then
            error = -error
          endif
          if (error .gt. cmp_rhs) then
            write(6,*) "[",mype,"]F_mpi(",i,")",array(i),"F_seq(",ivloc2sdgl(i),")",array_seq(ivloc2sdgl(i)),"diff",error
            call flush(6)
          end if
        enddo
          call flush(6)
          call VecRestoreArrayF90(Vseq,array_seq,ierr)
          call VecRestoreArrayF90(psp_F,array,ierr)
          call VecDestroy(Vseq,ierr)
      endif !(mype .eq. 0) then
#endif

!#if defined(DEBUG_rhs_2)
!       Compare all X and F values
!-------------------------------------
        call VecScatterCreateToZero(psp_X,scatter,Xmpi,ierr)
        call VecScatterBegin(scatter,psp_X,Xmpi,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(scatter,psp_X,Xmpi,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterDestroy(scatter,ierr)

        call VecScatterCreateToZero(psp_F,scatter,Vmpi,ierr)
        call VecScatterBegin(scatter,psp_F,Vmpi,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(scatter,psp_F,Vmpi,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterDestroy(scatter,ierr)
        !if (mype .eq. 0) then
        !  write(6,*)" Vmpi=psp_F:"
        !  call flush(6)
        !  call VecView(Vmpi,PETSC_VIEWER_STDOUT_SELF,ierr)
        !end if

        if (mype .eq. 0) then
          write(6,*)"[",mype,"] load X from datafile/X ..."
          call flush(6)
          call PetscViewerBinaryOpen(MPI_COMM_SELF,"datafile/X",FILE_MODE_READ,viewer,ierr)
          call VecLoad(viewer, VECSEQ, Xseq,ierr)
          call PetscViewerDestroy(viewer,ierr)
          !call VecView(Vseq,PETSC_VIEWEVR_STDOUT_SELF,ierr)

          call VecGetArrayF90(Xmpi,array,ierr)
          call VecGetArrayF90(Xseq,array_seq,ierr)
          write(6,*)"[",mype,"] Compare solution vector Xmpi and Xseq with tol",cmp_rhs
          call flush(6)
          do i=1,neqg
            id = (i-1)/n
            j  = i - id*n
            id = id+1
            !write(6,*) array(i),"X_np1(",ivloc2sdg(j,id),")",array_seq(ivloc2sdg(j,id))
            error = (array(i) - array_seq(ivloc2sdg(j,id)))
            if (error .lt. 0.0) then
              error = -error
            endif
            if (error .gt. cmp_rhs) then
              write(6,*) "Xmpi(",i,")",array(i),"Xseq(",ivloc2sdg(j,id),")",array_seq(ivloc2sdg(j,id)),"diff",error
            end if
          enddo
          call flush(6)
          call VecRestoreArrayF90(Xseq,array_seq,ierr)
          call VecRestoreArrayF90(Xmpi,array,ierr)

          write(6,*)" "
          write(6,*)"[",mype,"] load F from datafile/Frhs ..."
          call flush(6)
          call PetscViewerBinaryOpen(MPI_COMM_SELF,"datafile/Frhs",FILE_MODE_READ,viewer,ierr)
          call VecLoad(viewer, VECSEQ, Vseq,ierr)
          call PetscViewerDestroy(viewer,ierr)
          !call VecView(Vseq,PETSC_VIEWEVR_STDOUT_SELF,ierr)

          call VecGetArrayF90(Vmpi,array,ierr)
          call VecGetArrayF90(Vseq,array_seq,ierr)
          write(6,*)"[",mype,"] Compare rhs Fmpi and Fseq with tol",cmp_rhs
          call flush(6)
          do i=1,neqg
            id = (i-1)/n
            j  = i - id*n
            id = id+1
            !write(6,*) array(i),"F_np1(",ivloc2sdg(j,id),")",array_seq(ivloc2sdg(j,id))
            error = (array(i) - array_seq(ivloc2sdg(j,id)))/(array_seq(ivloc2sdg(j,id)))
            if (error .lt. 0.0) then
              error = -error
            endif
            if (error .gt. cmp_rhs) then
              write(6,*) "F_mpi(",i,")",array(i),"F_seq(",ivloc2sdg(j,id),")",array_seq(ivloc2sdg(j,id)),"rel_diff",error
              !write(6,*) "F_mpi(",i,")",array(i),"F_seq(",ivloc2sdg(j,id),")",array_seq(ivloc2sdg(j,id)),"diff",error
            end if
          enddo
          call flush(6)
          call VecRestoreArrayF90(Vseq,array_seq,ierr)
          call VecRestoreArrayF90(Vmpi,array,ierr)
          call VecDestroy(Vseq,ierr)
          call VecDestroy(Xseq,ierr)
        endif
        call VecDestroy(Xmpi,ierr)
        call VecDestroy(Vmpi,ierr)
!#endif
      endif

      psp_isf = 1
      call VecNorm(psp_F,NORM_2,norm,ierr)
      if (mype .eq. 0) then
        write(6,*)" initial Fnorm without sf scaling: ",norm
        call flush(6)
      end if
      call MPI_Barrier(comm,ierr)

      return
      end subroutine RhsCmp_Debug

! ---------------------------------------------------------------------
!  FormInitialGuess - Form initial approximation
!    Copy Uedge yl into X and scale X
!
!  Input Parameters:
!    neq  - local num of equations
!    yl   - initial approximate solution local array provided by Uedge
!    X    - vector
!
!  Output Parameter:
!    X - vector
! ---------------------------------------------------------------------
      subroutine FormInitialGuess(neq,yl,X,ierr)
      use Petscmod
      use PETSc_Snes_Param
      implicit none

      PetscInt,       INTENT(IN)    :: neq
      PetscScalar,    INTENT(IN)    :: yl(neq)
      Vec,            INTENT(INOUT) :: X
      PetscErrorCode, INTENT(OUT)   :: ierr
      PetscInt                      :: i
      PetscScalar,pointer           :: xx(:),su(:)

      call VecGetArrayF90(X,xx,ierr)
      if (psp_iscolnorm.gt.0) then
        call VecGetArrayF90(psp_SUVEC,su,ierr)
        do i=1,neq
           xx(i) = su(i)*yl(i)
        end do
        call VecRestoreArrayF90(psp_SUVEC,su,ierr)
      else
        do i=1,neq
           xx(i) = yl(i)
        end do
      endif
      call VecRestoreArrayF90(X,xx,ierr)

      end subroutine FormInitialGuess

! ---------------------------------------------------------------------
!  FormFinalSolution - Form final approximation for Uedge
!    Recover X with scaling suscal and convert it back to Uedge yl
!
!  Input Parameters:
!    n    - local num of equations
!    X    - vector
!
!  Output Parameter:
!    yl - Uedge solution local array without suscal scaling
! ---------------------------------------------------------------------
      subroutine FormFinalSolution(n,X,yl,ierr)
      use Petscmod
      use PETSc_Snes_Param
      implicit none

      PetscInt,       INTENT(IN)    :: n
      PetscScalar,    INTENT(OUT)   :: yl(n)
      Vec,            INTENT(IN)    :: X
      PetscErrorCode, INTENT(OUT)   :: ierr
      PetscInt                      :: i
      PetscScalar,pointer           :: xx(:),su(:)

      call VecGetArrayF90(X,xx,ierr)
      if (psp_iscolnorm.gt.0) then   
        call VecGetArrayF90(psp_SUVEC,su,ierr)     
        do i = 1,n
          yl(i) = xx(i)/su(i)
        end do   
        call VecRestoreArrayF90(psp_SUVEC,su,ierr)    
      else
        do i = 1,n
          yl(i) = xx(i)
        end do
      end if
      call VecRestoreArrayF90(X,xx,ierr)
      end subroutine FormFinalSolution

! ---------------------------------------------------------------------
!  FormUedgeFunction - Evaluates nonlinear function, F(x).
!
!    Note that psp_mon_ftime is never manually initialized.  It is automatically
!    initialized to 0 during the allocation.  If you're coming across a problem
!    it may be that the compiler doesn't do this.
!
!  Input Parameters:
!    snes - the SNES context
!    X    - input vector
!    pctx - defined in Petscmod to carry info to this function
!
!  Output Parameter:
!    F - function vector
! ---------------------------------------------------------------------
      subroutine FormUedgeFunction(snes,X,F,pctx,ierr)
      use Petscmod
      use PETSc_Snes_Param
      use PETSc_SNES_Param_Monitor
      implicit none

      SNES,       INTENT(IN)    :: snes
      Vec,            INTENT(IN)    :: X
      Vec,            INTENT(OUT)   :: F
      type(Petscctx), INTENT(INOUT) :: pctx
      PetscErrorCode, INTENT(OUT)   :: ierr

      integer             :: n,i
      PetscInt            :: nn
      PetscScalar,pointer :: xx(:),ff(:),sf(:),su(:),yl(:)
      PetscLogDouble      :: tstart, tend

      call PetscGetTime(tstart,ierr)

      yl=>pctx%ylpoint
      call VecGetLocalSize(X,nn,ierr)
      n = nn
      call VecGetArrayF90(X,xx,ierr)
      if (psp_iscolnorm.gt.0) then  ! scaling with ISUVEC
        call VecGetArrayF90(psp_SUVEC,su,ierr)
        do i=1,n
          yl(i) = xx(i)/su(i)  !scale with su^-1
        end do
        call VecRestoreArrayF90(psp_SUVEC,su,ierr)
      else
        do i=1,n
          yl(i) = xx(i)
        end do
      end if
      call VecRestoreArrayF90(X,xx,ierr)

      call VecGetArrayF90(F,ff,ierr)
      call rhsnk(n, yl, ff)

      if (psp_isf.eq.1) then
        call VecGetArrayF90(psp_SFVEC,sf,ierr)
        do i = 1,n
          ff(i) = sf(i)*ff(i)
        end do
        call VecRestoreArrayF90(psp_SFVEC,sf,ierr)
      end if
      call VecRestoreArrayF90(F,ff,ierr)

      call PetscGetTime(tend,ierr)
      if (psp_mon_on) then
        psp_mon_fcalls(psp_mon_arrlen) = psp_mon_fcalls(psp_mon_arrlen) + 1
        psp_mon_ftime(psp_mon_arrlen) = psp_mon_ftime(psp_mon_arrlen) + (tend-tstart)
      end if

      end subroutine FormUedgeFunction

!----------------------------------------------------------------------------------------
      subroutine PetscJac_calc(neq,t,yl,yldot0,ml,mu,wk,nnzmx,flag,Jmat,ierr)

!     Modified from jac_calc ()

! ... Common blocks:
            use Dim                      ! nx,ny,
                                         ! nusp[for fnorm not used here]
            use Timing                   ! istimingon,ttjstor,ttotjf,ttimpjf
            use Math_problem_size        ! neqmx,numvar
            use Grid                     ! ngrid,ig,ijac,ijactot
            use Indexes                  ! igyl,iseqalg
            use Variable_perturbation    ! del,dylconst
            use Jacobian_clipping        ! jaccliplim,istopjac,irstop,icstop
            use Jacobian_csc             ! rcsc,jcsc,icsc,yldot_pert
            use Ynorm                    ! suscal,sfscal
            use UEpar                    ! isphion,isnewpot,svrpkg,isbcwdt
            use Model_choice             ! iondenseqn
            use Imprad                   ! isimpon
            use Bcond                    ! isextrnpf,isextrtpf,isextrngc,
                                         ! isextrnw,isextrtw
            use Time_dep_nwt             ! nufak,dtreal,ylodt,dtuse
            use Selec                    ! yinc
            use Npes_mpi                 ! npes,mype,ismpion

            use Petscmod

! ... Functions:
      implicit none

! ... Input arguments:
      integer neq         ! local number of equations (all grid points)
      real :: t           ! physical time
      real :: yl(*)       ! dependent variables
      real :: yldot0(neq) ! right-hand sides evaluated at yl
      integer :: ml, mu   ! lower and upper bandwidths
      PetscInt :: nnzmx    ! maximum number of nonzeros in Jacobian

! ... Work-array argument:
      real :: wk(neq)     ! work space available to this subroutine

! ... Output arguments:
      real :: jac(nnzmx)     ! nonzero Jacobian elements
      integer :: ja(nnzmx)   ! col indices of nonzero Jacobian elements
      integer :: ia(neq+1)   ! pointers to beginning of each row in jac,ja

      logical tstguardc
      real(kind=4) :: gettime

! ... Local variables:
      integer :: nnz, ii, iv, ii1, ii2, xc, yc, ix, iy
      real :: yold, dyl, jacelem
      real(kind=4) :: sec4, tsjstor, tsimpjf, dtimpjf

      Mat ::            Jmat
      PetscErrorCode :: ierr
      PetscInt ::       npandf,rstart,rend
      PetscLogDouble :: ptstart,tfinal,etime
      MatStructure ::   flag

! ... Get initial value of system cpu timer.
      if (istimingon .eq. 1) tsjstor = gettime(sec4)

!      write(6,*) " PetscJac_calc: MatStructure flag =  ", flag,"nx,ny: ",nx,ny
!      call flush(6)
      call MatGetOwnershipRange(Jmat,rstart,rend,ierr)

! ... Set up diagnostic arrays for debugging
      do iv = 1, neq
        yldot_unpt(iv) = yldot0(iv)  ! for diagnostic only
        yldot_pert(iv) = 0.
      enddo

! ... Begin loop over dependent variables.
      npandf = 0
      etime  = 0.0
      nnz = 1
      do iv = 1, neq !for each local column

! ... Set beginning and ending indices of right-hand sides that might be perturbed.
         ii1 = max(iv-mu, 1)
         ii2 = min(iv+ml, neq)
! ... Reset range if this is a potential perturbation with isnewpot=1
         if (isphion*isnewpot.eq.1) then
            ii1 = max(iv-4*numvar*nx, 1)      ! 3*nx may be excessive
            ii2 = min(iv+4*numvar*nx, neq)    ! 3*nx may be excessive
         endif
! ... Reset range if extrapolation boundary conditions are used
         if (isextrnpf+isextrtpf+isextrngc+isextrnw+isextrtw.gt.0) then
            ii1 = max(iv-2*numvar*(nx+3), 1)      ! guess to include extrap. bc
            ii2 = min(iv+2*numvar*(nx+3), neq)    ! guess to include extrap. bc
         endif

! ... Initialize all of those right-hand sides to their unperturbed values.
         do ii = ii1, ii2   ! below wk is reset, but only over limited range
            wk(ii) = yldot0(ii)
         enddo

! ... Set spatial-location indices for this dependent variable.
         xc = igyl(iv,1)
         yc = igyl(iv,2)
!         write(6,*)"[",mype,"] iv,xc, yc: ",iv,xc,yc,"rstart,rend: ",rstart,rend
!         call flush(6)

! ... Save value of dependent variable, then perturb it.
!    The perturbation to the variable is proportional to parameter
!    del and to a measure of the size of the variable.  That measure
!    increases with the absolute value of the variable if it exceeds
!    the typical size given by dylconst/suscal but can never be less
!    than that typical size.
         yold = yl(iv)
         dyl = del * (abs(yold) + dylconst / suscal(iv))
         yl(iv) = yold + dyl

! ... Calculate right-hand sides near location of perturbed variable.

         call PetscGetTime(ptstart,ierr)
         call pandf1 (xc, yc, iv, neq, t, yl, wk)
         call PetscGetTime(tfinal,ierr)
         npandf = npandf + 1
         etime  = etime + (tfinal - ptstart)

! ... Calculate possibly nonzero Jacobian elements for this variable,
!c     and store nonzero elements in compressed sparse column format.
         jcsc(iv) = nnz      ! sets index for first Jac. elem. of var. iv
         do ii = ii1, ii2    !for each local row
            jacelem = (wk(ii) - yldot0(ii)) / dyl

! ...  Add diagonal 1/dt for nksol - PetscSnes()
            if (iv.eq.ii) then !diagonal entry
              if (iseqalg(iv)*(1-isbcwdt).eq.0) then
                jacelem = jacelem - 1/dtuse(iv)
              endif
              ix = igyl(iv,1)
              iy = igyl(iv,2)
              if (idxphi(ix,iy)==iv .and. dtphi<1e10) then !selects phi eqn
                jacelem = jacelem - 1/dtphi
              endif
            endif

! ...  Add a pseudo timestep to the diagonal !! if eqn is not algebraic
            if (nufak .gt. 0) then
               if (iv.eq.ii .and. yl(neq+1).eq.1) then
                 jacelem = jacelem - nufak  !omit .and. iseqalg(iv).eq.0)
               endif
            endif

          if (flag .eq. DIFFERENT_NONZERO_PATTERN) then
! ...       Add a nonzero entry into Jmat
            if (abs(jacelem*sfscal(iv)) .gt. jaccliplim) then
               if (nnz .gt. nnzmx) then
                  write(6,*) '*** jac_calc - Storage exceeded at (i,j)=',ii,iv,' nnz=',nnz,' nnzmx=',nnzmx
                 call xerrab("")
               endif
!               if (rdoff.ne.0.e0) jacelem=jacelem*(1.0e0+ranf()*rdoff)
               rcsc(nnz) = jacelem
               icsc(nnz) = ii
               nnz       = nnz + 1
#if defined(DEBUG_SNES)
               if (mype .eq. 1) then
                 write(6,*) "[",mype,"] MatSetValues at ",ii-1+rstart,iv-1+rstart
                 call flush(6)
               endif
#endif
               call MatSetValues(Jmat,1,ii-1+rstart,1,iv-1+rstart,jacelem,INSERT_VALUES,ierr)
            endif
          else
            write(6,*) "use same NONZERO_PATTERN"
          end if

            if (istopjac.gt.0 .and. ii.eq.irstop .and. iv.eq.icstop) then
               yldot_pert(ii) = wk(ii)      ! for diagnostic only
               if (istopjac == 2) then
                 yl(iv) = yold

                 call PetscGetTime(ptstart,ierr)
                 call pandf1 (xc, yc, iv, neq, t, yl, wk)
                 call PetscGetTime(tfinal,ierr)
                 etime  = etime + (tfinal - ptstart)
                 npandf = npandf + 1

               endif

               call remark("***** non-zero jac_elem at irstop,icstop")
               write(6,*) 'irstop = ', irstop, ', icstop = ', icstop
               call xerrab("")
            endif

         enddo  ! end of ii loop over equations

! ... Restore dependent variable and plasma variables near its location.
         yl(iv) = yold

         call PetscGetTime(ptstart,ierr)
         call pandf1 (xc, yc, iv, neq, t, yl, wk)
         call PetscGetTime(tfinal,ierr)
         etime  = etime + (tfinal - ptstart)
         npandf = npandf + 1

! ... End loop over dependent variables and finish Jacobian storage.
! 18   continue
!...  If this is the last variable before jumping to new cell, reset pandf
         if (mod(iv,numvar).eq.0 .and. isjacreset.ge.1) then
            call PetscGetTime(ptstart,ierr)
            call pandf1 (xc, yc, iv, neq, t, yl, wk)
            call PetscGetTime(tfinal,ierr)
            etime  = etime + (tfinal - ptstart)
            npandf = npandf + 1
         endif
      enddo             ! end of main loop over yl variables
      call MatAssemblyBegin(Jmat,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(Jmat,MAT_FINAL_ASSEMBLY,ierr)

      jcsc(neq+1) = nnz

! ... Accumulate cpu time spent here.
      if (istimingon .eq. 1) ttjstor = ttjstor + gettime(sec4) - tsjstor

      return
      end subroutine PetscJac_calc

!-----------------------------------------------------------------------
!  This is the subroutine to sequentially compute the Jacobian for PETSc purposes
!  It utilizes the jac_calc routine from the UEDGE package
!
!  Input Parameters:
!    snes     - SNES context for solving the problem
!    X        - Value at which the Jacobian will be evaluated
!    pctx     - Petscmod defined data structure to get us information
!
!  Output Parameters:
!    Jmat     - Jacobian matrix
!    jac_prec - Preconditioner matrix for the nonlinear system (may be same as Jmat)
!    flag     - Nonzero Pattern flag
! ---------------------------------------------------------------------
      subroutine FormUedgeJacobian(snes,X,Jmat,jac_prec,flag,pctx,ierr)
      use Petscmod
      use Decomp        ! lbw,ubw
      use Time_dep_nwt  ! nufak,ydt_max,ydt_max0,alfnuf,expnuf,nufak0
                        ! inufaknk,dtoptx,dtoptv
      use PETSc_Snes_Param
      implicit none

      integer*8 :: temp1, temp2

      PetscInt ::        snes
      Vec ::             X,F
      Mat ::             Jmat,jac_prec
      MatStructure ::    flag
      type(Petscctx) ::  pctx
      PetscErrorCode ::  ierr
      PetscLogDouble ::  tstart,tfinal,etime

! ... Local variables:
      PetscReal        tp
      PetscInt         i,nn,nnzmx,incpset
      integer          n
      PetscScalar,pointer :: xx(:),f0(:),wk(:),su(:),ff(:),yl(:)

      write(6,*) "    build UedgeJacobian ..."
      call flush(6)
! ... Now we need to allocate space so that the UEDGE functions can do their thing
! ... Also we need to move the variables to Fortran variables
      wk=>pctx%wkpoint
      yl=>pctx%ylpoint
      nnzmx   = pctx%nnzmx
      incpset = pctx%incpset

      call VecGetLocalSize(X,nn,ierr)
      n = nn

!  Copy input vector array into yl (for n+2 considerations)
      call VecGetArrayF90(X,xx,ierr)
      if (psp_iscolnorm .ge. 0) then
        call VecGetArrayF90(psp_SUVEC,su,ierr)
        do i=1,n
          yl(i) = xx(i)/su(i)  !scale with su^-1
        end do
        call VecRestoreArrayF90(psp_SUVEC,su,ierr)
      else
        do i=1,n
          yl(i) = xx(i)
        end do
      end if
      call VecRestoreArrayF90(X,xx,ierr)

! ... Calculate ydt_max = max(f*sf) to control yl(n+2) = nufak (was in psetnk)
      call VecNorm(psp_F,NORM_INFINITY,ydt_max,ierr)
      if (ydt_max < 1.e-16) ydt_max = 1.e-16
      if (ydt_max0 .eq. 0) ydt_max0 = ydt_max
      nufak = min(nufak*alfnuf*(ydt_max/ydt_max0)**expnuf, nufak0)
      if (inufaknk .eq. 1) then   ! deter. if nufak is used in Krylov step
        yl(n+2) = nufak
      else
        yl(n+2) = 0.
      endif
      if (expnuf.ne.0.) write(6,*) ' nufak = ', nufak
      ydt_max0 = ydt_max

! ... Call pandf to set terms with yl(n+1) flag on for Jacobian
      f0=>pctx%savfpoint
      call rhsnk(n, yl, f0)

! ... Calculate Jacobian matrix.
      tp = 0.
!      call PetscGetTime(tstart,ierr)
      flag=DIFFERENT_NONZERO_PATTERN
      call Petscjac_calc (n,tp,yl,f0,lbw,ubw,wk,nnzmx,flag,jac_prec,ierr)
      if (Jmat .ne. jac_prec) then
        call MatAssemblyBegin(Jmat,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(Jmat,MAT_FINAL_ASSEMBLY,ierr)
      endif
!      call PetscGetTime(tfinal,ierr)
!      write(6,*) "    Time for Petscjac_calc: ",tfinal-tstart

!  Do the sf scaling
      if (psp_isf.eq.1) then
        call MatDiagonalScale(jac_prec,psp_SFVEC,PETSC_NULL_OBJECT,ierr)
      end if
!  Do ISUVEC=su^-1 scaling.
      if (psp_iscolnorm .gt. 0) then
        call MatDiagonalScale(jac_prec,PETSC_NULL_OBJECT,psp_ISUVEC,ierr)
      end if
      end subroutine FormUedgeJacobian

! ---------------------------------------------------------------------
!  This function imposes the constraints required by UEDGE on the variables
!  These constraints are found in icnstr, icflag, rlx
!  If icflag==1, icnstr is imposed: icnstr(i)>0 -> yl(i) = x(i) - minchange*y(i) >0
!                                   icnstr(i)<0 -> yl(i) = x(i) - minchange*y(i) <0
!                                   icnstr(i)==0 -> no constraint
!
!  rlx defines the inequality: abs(u_step/u)<rlx
!  If the rlx constraint is violated, the step size is scaled back to satisfy it
!
!  Input Parameters:
!    snes - SNES context for solving the problem
!    X    - Current Newton iterate
!    Y    - Current Newton step
!    pctx - Petscmod defined structure to pass data to this function
!
!  Output Parameters:
!    changed_y      - Boolean saying that the length of search has been changed
! ---------------------------------------------------------------------
      subroutine FormPreCheck(snes,X,Y,pctx,changed_y,ierr)
      use Petscmod
      use Npes_mpi         ! npes,mype,ismpion
      implicit none

      SNES,           INTENT(IN)    :: snes
      Vec,            INTENT(IN)    :: X
      Vec,            INTENT(INOUT) :: Y
      type(Petscctx), INTENT(INOUT) :: pctx
      PetscErrorCode, INTENT(OUT)   :: ierr
      PetscScalar,    INTENT(OUT)   :: changed_y

      PetscInt ::       i,neq
      PetscReal ::      rlx
      PetscInt,pointer :: icnstr(:)
      PetscScalar ::    minchange,w,tiny,sbuf(3),rbuf(3)
      integer4 ::       comm
      PetscScalar,pointer :: xp(:),yp(:)

      minchange = 1.0
      tiny      = 1.0e-16

      icnstr=>pctx%icnstrpoint
      rlx = pctx%rlx

      call PetscObjectGetComm(snes,comm,ierr)

!  As of right now, this function will never changed the search direction, just the length
      changed_y = 0.0

      call VecGetLocalSize(X,neq,ierr)
      call VecGetArrayF90(X,xp,ierr)
      call VecGetArrayF90(Y,yp,ierr)
      do i=1,neq
        if (icnstr(i).ne.0) then
          w = xp(i) - yp(i)
          if ((icnstr(i).gt.0 .and. w.le.0) .or. (icnstr(i).lt.0 .and. w.ge.0)) then
            changed_y = 1.0
            minchange = min(minchange,xp(i)/yp(i))    ! yp(i) cannot be 0
          end if
        end if
      end do
      call VecRestoreArrayF90(X,xp,ierr)
      call VecRestoreArrayF90(Y,yp,ierr)

      call MPI_Comm_size(comm,npes,ierr)
      if (npes .gt. 1) then
        sbuf(1) = changed_y
        sbuf(2) = 1.0/minchange
        call MPI_ALLREDUCE(sbuf,rbuf,2,MPIU_SCALAR,MPI_MAX,comm,ierr)
        changed_y = rbuf(1)
        minchange = 1.0/rbuf(2)
      end if

!...  scale direction Y by minchange
      if (changed_y .eq. 1.0) then
        minchange = rlx*minchange + tiny
#if defined(DEBUG_SNES)
        if (mype .eq. 0) then
          write(6,*)  "  FormPreCheck: scale direction Y with minchange",minchange,"rlx",rlx
          call flush(6)
        end if
#endif
        call VecScale(Y,minchange,ierr)
      end if

!      if(pctx%countoutloud.eq.1) write(6,*) "FormPreCheck called for ",pctx%icniter,"time"
!      pctx%icniter=pctx%icniter+1
      end subroutine FormPreCheck

! ---------------------------------------------------------------------
!  subroutine psp_monAdjustSize(sizeIncrease,ierr)
!
!  This subroutine allows for an adjustment in the size of the monitor
!  arrays.  Because we don't know ahead of time how many SNES steps will
!  be needed, these arrays may need to be lengthened over time.  Of course
!  if we were using a real programming language there would be an object
!  that would do this for us.  But whatever...
!
!  Also, now passing sizeIncrease=0 will delete any data in 
!
!  Note that I think values are always initialized to 0 by Fortran after
!  allocating.  Even though this is my belief I still initialize ftime
!  and fcalls to 0 in the first array location.  Then future initializations
!  take place in FormMonitor.
! ---------------------------------------------------------------------
      subroutine psp_monAdjustSize(sizeIncrease,ierr)
      use Petscmod
      use PETSc_Snes_Param
      use PETSc_Snes_Param_Monitor
      implicit none

      PetscInt,       INTENT(IN)  :: sizeIncrease
      PetscErrorCode, INTENT(OUT) :: ierr
      
      integer i,allocInc,arrLen ! I think these are needed cuz Fortran's not cool with Petsc vars

      PetscReal,      ALLOCATABLE, DIMENSION(:) :: psp_mon_fnrm_tmp
      PetscReal,      ALLOCATABLE, DIMENSION(:) :: psp_mon_dtreal_tmp
      PetscLogDouble, ALLOCATABLE, DIMENSION(:) :: psp_mon_time_tmp
      PetscLogDouble, ALLOCATABLE, DIMENSION(:) :: psp_mon_ftime_tmp
      PetscInt,       ALLOCATABLE, DIMENSION(:) :: psp_mon_kspits_tmp
      PetscInt,       ALLOCATABLE, DIMENSION(:) :: psp_mon_fevals_tmp
      PetscInt,       ALLOCATABLE, DIMENSION(:) :: psp_mon_fcalls_tmp
      PetscInt,       ALLOCATABLE, DIMENSION(:) :: psp_mon_snesit_tmp

      arrLen = SIZE(psp_mon_fnrm) ! All arrays have the same length
      allocInc = sizeIncrease

      ! If we are decreasing the size of the arrays we need to consider that
      ! by only allocating the space we need
      arrLen = MIN(arrLen,arrLen+allocInc)

      ! If the user passes a negative value smaller than the length of the data
      ! this makes sure the function just deletes everything rather than a
      ! negative allocation which may have unstable results
      if (arrLen .le. 0) then
        allocInc = 0
      end if

      if (allocInc .eq. 0) then ! Delete the existing data
        if (ALLOCATED(psp_mon_fnrm)) then ! Make sure there is data to be deleted
          DEALLOCATE(psp_mon_fnrm)
          DEALLOCATE(psp_mon_dtreal)
          DEALLOCATE(psp_mon_time)
          DEALLOCATE(psp_mon_ftime)
          DEALLOCATE(psp_mon_kspits)
          DEALLOCATE(psp_mon_fevals)
          DEALLOCATE(psp_mon_fcalls)
          DEALLOCATE(psp_mon_snesit)
        end if
      else if (ALLOCATED(psp_mon_fnrm)) then ! Anything except the first SNES solve
        ALLOCATE(psp_mon_fnrm_tmp(1:arrLen))
        ALLOCATE(psp_mon_dtreal_tmp(1:arrLen))
        ALLOCATE(psp_mon_time_tmp(1:arrLen))
        ALLOCATE(psp_mon_ftime_tmp(1:arrLen))
        ALLOCATE(psp_mon_kspits_tmp(1:arrLen))
        ALLOCATE(psp_mon_fevals_tmp(1:arrLen))
        ALLOCATE(psp_mon_fcalls_tmp(1:arrLen))
        ALLOCATE(psp_mon_snesit_tmp(1:arrLen))
        do i=1,arrLen
          psp_mon_fnrm_tmp(i)   = psp_mon_fnrm(i)
          psp_mon_dtreal_tmp(i) = psp_mon_dtreal(i)
          psp_mon_time_tmp(i)   = psp_mon_time(i)
          psp_mon_ftime_tmp(i)  = psp_mon_ftime(i)
          psp_mon_kspits_tmp(i) = psp_mon_kspits(i)
          psp_mon_fevals_tmp(i) = psp_mon_fevals(i)
          psp_mon_fcalls_tmp(i) = psp_mon_fcalls(i)
          psp_mon_snesit_tmp(i) = psp_mon_snesit(i)
        end do
        DEALLOCATE(psp_mon_fnrm)
        DEALLOCATE(psp_mon_dtreal)
        DEALLOCATE(psp_mon_time)
        DEALLOCATE(psp_mon_ftime)
        DEALLOCATE(psp_mon_kspits)
        DEALLOCATE(psp_mon_fevals)
        DEALLOCATE(psp_mon_fcalls)
        DEALLOCATE(psp_mon_snesit)
        
        ! Now we need to move the data back from the temporary arrays
        
        ALLOCATE(psp_mon_fnrm(1:arrLen+allocInc))
        ALLOCATE(psp_mon_dtreal(1:arrLen+allocInc))
        ALLOCATE(psp_mon_time(1:arrLen+allocInc))
        ALLOCATE(psp_mon_ftime(1:arrLen+allocInc))
        ALLOCATE(psp_mon_kspits(1:arrLen+allocInc))
        ALLOCATE(psp_mon_fevals(1:arrLen+allocInc))
        ALLOCATE(psp_mon_fcalls(1:arrLen+allocInc))
        ALLOCATE(psp_mon_snesit(1:arrLen+allocInc))
        do i=1,arrLen ! Still use the old
          psp_mon_fnrm(i)   = psp_mon_fnrm_tmp(i)
          psp_mon_dtreal(i) = psp_mon_dtreal_tmp(i)
          psp_mon_time(i)   = psp_mon_time_tmp(i)
          psp_mon_ftime(i)  = psp_mon_ftime_tmp(i)
          psp_mon_kspits(i) = psp_mon_kspits_tmp(i)
          psp_mon_fevals(i) = psp_mon_fevals_tmp(i)
          psp_mon_fcalls(i) = psp_mon_fcalls_tmp(i)
          psp_mon_snesit(i) = psp_mon_snesit_tmp(i)
        end do
        DEALLOCATE(psp_mon_fnrm_tmp)
        DEALLOCATE(psp_mon_dtreal_tmp)
        DEALLOCATE(psp_mon_time_tmp)
        DEALLOCATE(psp_mon_ftime_tmp)
        DEALLOCATE(psp_mon_kspits_tmp)
        DEALLOCATE(psp_mon_fevals_tmp)
        DEALLOCATE(psp_mon_fcalls_tmp)
        DEALLOCATE(psp_mon_snesit_tmp)
      else ! If this is the first SNES solve, the arrays aren't allocated
        ALLOCATE(psp_mon_fnrm(1:allocInc))
        ALLOCATE(psp_mon_dtreal(1:allocInc))
        ALLOCATE(psp_mon_time(1:allocInc))
        ALLOCATE(psp_mon_ftime(1:allocInc))
        ALLOCATE(psp_mon_kspits(1:allocInc))
        ALLOCATE(psp_mon_fevals(1:allocInc))
        ALLOCATE(psp_mon_fcalls(1:allocInc))
        ALLOCATE(psp_mon_snesit(1:allocInc))
        psp_mon_fcalls(1) = 0
        psp_mon_ftime(1) = 0
      end if

      end subroutine psp_monAdjustSize
      
! ---------------------------------------------------------------------
!  subroutine DetermineJacRebuild
!  This subroutine looks at the evolution of the time stepping and the
!  SNES solve to determine if the Jacobian needs to be rebuilt
!  It will internally set the function SNESSetLagJacobian, and it will
!  also update values in snesmonitor.  It returns a string output
!  to the screen.
! ---------------------------------------------------------------------
      subroutine DetermineJacRebuild(snes,snesmonitor,its,rebuild,jac_rebuild_output,ierr)
      use Snesmonitormod
      use PETSc_Snes_Param ! Anything with a psp_ header is stored in here
      implicit none
      
      SNES,           INTENT(INOUT) :: snes
      type(monctx),   INTENT(INOUT) :: snesmonitor
      PetscInt,       INTENT(IN)    :: its
      PetscInt,       INTENT(OUT)   :: rebuild
      character*80,   INTENT(OUT)   :: jac_rebuild_output
      PetscErrorCode, INTENT(OUT)   :: ierr
      
      PetscReal :: dt_factor,max_dtchange
      PetscInt  :: newdt
      
      max_dtchange = psp_jacMaxdtChange
      if (max_dtchange .le. 1e-15) then
        newdt = 0 ! max_dt=0 means never reevaluate for time step
      else
        dt_factor = max(snesmonitor%thisdt/psp_olddt,psp_olddt/snesmonitor%thisdt)
        if (dt_factor .le. max_dtchange) then
          newdt = 0
        else
          newdt = 1
        end if
      end if

      rebuild = -1
!      if (its.gt.1) then
!        if (mod(its-1,psp_jaclag).eq.0) then
!          rebuild = 1
!        end if
!      end if
      if ((its.gt.0) .and. (mod(its,psp_jaclag).eq.0)) then
        rebuild = 1
      end if

      if ((psp_skipJacobian .eq. 1) .and. (its.eq.0)) then
        if (psp_nsnes .eq. 0) then
          rebuild = 2
        elseif (newdt .eq. 1) then
          rebuild = 5
        end if
      elseif (its.eq.0) then
        rebuild = 3
      end if
      if (snesmonitor%buildJacobian .eq. 1) then
        snesmonitor%buildJacobian = 0 ! Reset after linear solve failure
        rebuild = 4
      end if

      if (rebuild .gt. 0) then
        call SNESSetLagJacobian(snes,1,ierr)
        snesmonitor%JacWasJustRebuilt = 1
      else
        call SNESSetLagJacobian(snes,-1,ierr)
        snesmonitor%JacWasJustRebuilt = 0
      end if

      if (rebuild .eq. 1) then
        write(jac_rebuild_output,*) "Rebuilding Jacobian, Lag =",psp_jaclag
      elseif (rebuild .eq. 2) then
        write(jac_rebuild_output,*) "Building Jacobian for first time step"
      elseif (rebuild .eq. 3) then
        if (psp_nsnes .eq. 0) then
         print *,jac_rebuild_output
         write(jac_rebuild_output,*) "Building Jacobian for first time step"
        else
          write(jac_rebuild_output,*) "Rebuilding Jacobian, Required at each SNES step"
        endif
      elseif (rebuild .eq. 4) then
        write(jac_rebuild_output,*) "Rebuilding Jacobian, Previous Linear Solve Failed"
      elseif (rebuild .eq. 5) then
        write(jac_rebuild_output,1001) snesmonitor%thisdt,psp_olddt,max_dtchange
1001   format("Rebuilding Jacobian: newdt=",ES8.2E2,", olddt=",ES8.2E2,", max_dtfactor=",F6.2)
      end if
      
      end subroutine DetermineJacRebuild

! ---------------------------------------------------------------------
!  Subroutine FormMonitor
!  This function lets up keep track of the SNES progress at each step
!  Eventually all the data recorded by this will be output to a file
!  This function is only used if the option
!     -print_stats <string>
!  is used at runtime from the file petscopt
!  string should be the name of the variable that you want to use in Matlab
!  The Matlab file will be f<string>.m
!
!  This function calls DetermineJacRebuild and passes that infor along to
!  SNES to tell it when the Jacobian needs rebuilt.  No PETSc specific
!  Jacobian rebuilding options are valid; any options are handled manually
!  in PetscSnes and then dealt with here.
!
!  At the end of each SNES step, this function is called by SNES.  The
!  psp_mon_ftime and psp_mon_fcalls values for the next SNES solve are
!  initialized to 0.
!
!  Input Parameters:
!    snes    - SNES nonlinear solver context
!    its     - current nonlinear iteration
!    norm    - 2-norm of current residual (may be approximate)
!    snesmonitor - monctx designed module (included in Snesmonitormod)
! ---------------------------------------------------------------------
      subroutine FormMonitor(snes,its,norm,snesmonitor,ierr)
      use Npes_mpi          ! npes,mype,ismpion - rm later!
      use Snesmonitormod
      use PETSc_Snes_Param
      use PETSc_Snes_Param_Monitor
      implicit none

      SNES,       INTENT(INOUT) :: snes
      PetscInt,       INTENT(IN)    :: its
      PetscScalar,    INTENT(IN)    :: norm
      type(monctx),   INTENT(INOUT) :: snesmonitor
      PetscErrorCode, INTENT(OUT)   :: ierr
      
      PetscInt     :: rebuild
      KSP          :: ksp
      character*80 :: jac_rebuild_output

      if (psp_mon_on) then
        call SNESGetKSP(snes,ksp,ierr)
        psp_mon_fnrm(psp_mon_arrlen) = norm
        psp_mon_dtreal(psp_mon_arrlen) = snesmonitor%thisdt
        psp_mon_snesit(psp_mon_arrlen) = its
        call PetscGetTime(psp_mon_time(psp_mon_arrlen),ierr)
        if (its .eq. 0) then
          psp_mon_kspits(psp_mon_arrlen) = 0
        else
          call KSPGetIterationNumber(ksp,psp_mon_kspits(psp_mon_arrlen),ierr)
        end if
        call SNESGetNumberFunctionEvals(snes,psp_mon_fevals(psp_mon_arrlen),ierr)

        psp_mon_arrlen = psp_mon_arrlen + 1
        if (psp_mon_arrlen .gt. SIZE(psp_mon_fnrm)) then
          call psp_monAdjustSize(psp_mon_allocInc,ierr)
        end if
        psp_mon_fcalls(psp_mon_arrlen) = 0
        psp_mon_ftime(psp_mon_arrlen) = 0
      end if

      call DetermineJacRebuild(snes,snesmonitor,its,rebuild,jac_rebuild_output,ierr)
      if (rebuild.gt.0 .and. mype.eq.0) then
        write(6,*) jac_rebuild_output
        call flush(6)
      end if

      end subroutine FormMonitor

! ---------------------------------------------------------------------
!  Subroutine PetscSnesSetSf
!     Calculate the scaling vector of function values.
!     The Jacobian at X is calculated, and sf(i) is set
!     to the reciprocal of the max_norm of row i.
!     This routine is written based on bbb/oderhs.m/sfsetnk()
!
!  Input Parameters:
!    snes - SNES nonlinear solver context
!    X    - vector at which the Jacobian is evaluated
!  Output Parameter:
!    SF   - scaling vector for the rhs function
! ---------------------------------------------------------------------
      subroutine PetscSnesSetSf(snes,X,SF,ierr)
      use Petscmod
      use PETSc_Snes_Param
      use Npes_mpi         ! npes,mype,ismpion - rm later!
      implicit none

      SNES,           INTENT(INOUT) :: snes
      Vec,            INTENT(INOUT) :: X,SF
      PetscErrorCode, INTENT(OUT)   :: ierr

      MatStructure                        :: Jstruct
      PetscInt                            :: i,n,nt
      PetscInt, ALLOCATABLE, DIMENSION(:) :: idx
      PetscReal                           :: norm
      PetscScalar,pointer                 :: xx(:),xxtmp(:)

      call SNESComputeJacobian(snes,X,psp_Jmat,psp_Jmat,Jstruct,ierr)
#if defined(DEBUG_SF)
      write(6,*) "PetscSnesSetSf- Jmat: "
      call flush(6)
      call MatView(psp_Jmat,PETSC_VIEWER_STDOUT_SELF,ierr)
#endif
      call VecGetSize(SF,nt,ierr)
      ALLOCATE(idx(1:nt))
      call MatGetDiagonal(psp_Jmat,psp_Xtmp,ierr)
      call MatGetRowMaxAbs(psp_Jmat,SF,idx,ierr)
      DEALLOCATE(idx)
#if defined(DEBUG_SNES)
      if (mype .eq. 0) then 
        write(6,*) "SF: "
        flush(6)
      endif
      call VecView(SF,PETSC_VIEWER_STDOUT_WORLD,ierr)
      write(6,*) "[",mype,"] PetscSnesSetSf() SNESComputeJacobian is done "
      call VecNorm(SF,NORM_2,norm,ierr)
      write(6,*) "[",mype,"] ...MaxRowNorm",norm
      call flush(6)
#endif

      call VecGetLocalSize(SF,n,ierr)
      call VecGetArrayF90(SF,xx,ierr)
        do i=1,n
          if (xx(i) .ne. 0.0) then
            xx(i) = 1.0/xx(i)
            if (xx(i) .lt. 0.0) xx(i) = - xx(i)
          else
            write(6,*) "[",mype,"] ***Error: Jrowmax(",i,")=",xx(i)
            call flush(6)
            SETERRQ(1,' ',ierr)
          endif
        end do
      call VecRestoreArrayF90(SF,xx,ierr)
      end subroutine PetscSnesSetSf

!************************************************************
!  subroutine PetscSnesFillSf()
!
! -------------------------------------------------------------------
      subroutine PetscSnesFillSf(neq,sfscal,SF,ierr)
      use Petscmod
      use PETSc_Snes_Param
      implicit none

      PetscInt,       INTENT(IN)    :: neq
      PetscScalar,    INTENT(IN)    :: sfscal(neq)
      Vec,            INTENT(INOUT) :: SF
      PetscErrorCode, INTENT(OUT)   :: ierr
      
      PetscInt            :: i
      PetscScalar,pointer :: xx(:)

      call VecGetArrayF90(SF,xx,ierr)
      do i=1,neq
        xx(i) = sfscal(i)
      end do
      call VecRestoreArrayF90(SF,xx,ierr)
      end subroutine PetscSnesFillSf

!************************************************************
!  subroutine PetscSnesFillSu()
!
! -------------------------------------------------------------------
      subroutine PetscSnesFillSu(neq,suscal,SU,ISU,ierr)
      use Petscmod
      use PETSc_Snes_Param
      implicit none

      PetscInt,       INTENT(IN)    :: neq
      PetscScalar,    INTENT(INOUT) :: suscal(neq)
      Vec,            INTENT(INOUT) :: SU,ISU
      PetscErrorCode, INTENT(OUT)   :: ierr
      
      PetscInt            :: i
      PetscScalar,pointer :: xx(:),xx1(:)

      call VecGetArrayF90(ISU,xx,ierr)
      call VecGetArrayF90(SU,xx1,ierr)
      do i=1,neq
        if (suscal(i) .gt. -0.000001) then
          if (suscal(i) .lt. 0.000001) then
            write(6,*)"su",i,suscal(i)," reset su=1.e-6!"
            call flush(6)
            suscal(i) = 0.000001
          end if
        end if
        xx1(i) = suscal(i)
        xx(i) = 1.0/suscal(i)
      end do
      call VecRestoreArrayF90(SU,xx1,ierr)
      call VecRestoreArrayF90(ISU,xx,ierr)
!      write(6,*)"SUVEC: "
!      call flush(6)
!      call VecView(SU,PETSC_VIEWER_STDOUT_WORLD,ierr)
      end subroutine PetscSnesFillSu

!************************************************************
!  subroutine PetscSnesApply()
!
! -------------------------------------------------------------------
      subroutine PetscSnesApply(n,yl,savf,suscal,sfscal,nnzmx,ireorder,incpset,icflag,rlx,icnstr,dtreal,iterm,ierr)
      use Petscmod
      use Petscmodinterfaces
      use Snesmonitormod
      use PETSc_Snes_Param
      implicit none

      PetscInt ::           n,nnzmx,ireorder,incpset
      PetscScalar,target :: yl(n),savf(n),suscal(n),sfscal(n)
      PetscInt,target ::    icnstr(n)
      PetscInt ::           icflag
      PetscReal ::          rlx,dtreal
      PetscInt ::           iterm

      PetscInt ::           i,its
      MatStructure ::       flag
      PetscErrorCode ::     ierr
      PetscInt ::           ksp
      PC ::                 pc
      Vec ::                X1
      PetscScalar,target :: wk(n+2)
      PetscViewer ::        output

      type(Petscctx), pointer:: pctx

      type(monctx) ::       snesmonitor
      SNESConvergedReason :: snesreason

      PetscScalar ::        neg_one
      PetscReal ::          norm
      PetscScalar, pointer ::  xx(:),xx1(:)
      real ::               gettime, sec4

      external FormPreCheck,FormMonitor

!      write(6,*)"*********************************"
!      write(6,*)"Apply: neq=",n
!      write(6,*)"................................."
!      call flush(6)

!      call SNESGetApplicationContext(psp_snes,pctx,ierr)

      return
      end subroutine PetscSnesApply

!************************************************************
!  subroutine PetscSnesCreateObjs()
!
! -------------------------------------------------------------------
      subroutine PetscSnesCreateObjs(comm,n,neqg,nnzmx,ierr)

      use Petscmod
      use PETSc_Snes_Param
      use Npes_mpi         ! npes,mype,ismpion
      use Dim              ! nx,ny    - rm later!
      use Math_problem_size  ! neqmx,numvar - rm later!
      Use Indices_loc_glob_map ! ivcum,ivloc2sdg,ivloc2mdg,ivl2gstnl - rm later!

      implicit none

      integer4 ::           comm
      PetscInt ::           n,neqg,nnzmx
      PetscErrorCode ::     ierr

      PetscInt ::           i,bs,rstart,rend,j,k,id
      PetscViewer ::        viewer
      PetscReal  ::         norm
      PetscInt    ::        mycols(neqg),nmycols,band,row,col,rows(numvar),cols(9*numvar),bcols(9),bj
      PetscScalar ::        zeros(9*numvar*numvar),zero
      PetscTruth ::         flg

!      write(6,*)"[",mype,"] Create: neq, neqg=",n,neqg," numvar",numvar
!      write(6,*)"................................."
!      call flush(6)

!     Create SNES objects 
!----------------------------------------------------------
      call SNESCreate(comm,psp_snes,ierr)

!     Create vectors for holding solutions and rhsfunc
!----------------------------------------------------
      ! write(6,*) "PetscSnesCreateObjs: comm =", comm
      ! call flush(6)
      call VecCreate(comm,psp_F,ierr)
      call VecSetSizes(psp_F,n,neqg,ierr)
      call VecSetFromOptions(psp_F,ierr)
      call VecDuplicate(psp_F,psp_X,ierr)
      call VecDuplicate(psp_F,psp_Xtmp,ierr)

!     Create scaling vectors
!--------------------------------------------
      if (psp_iscolnorm .gt. 0) then
        call VecDuplicate(psp_X,psp_ISUVEC,ierr)
        call VecDuplicate(psp_X,psp_SUVEC,ierr)
      endif

      if (psp_isf.eq.1) then
        call VecDuplicate(psp_F,psp_SFVEC,ierr)
      endif

!     Create Jacobian matrix Jmat used as preconditioner
!-------------------------------------------------------
      call MatCreate(comm,psp_Jmat,ierr)
      call MatSetSizes(psp_Jmat,n,n,PETSC_DECIDE,PETSC_DECIDE,ierr)

!     default JmatType:
      if (npes .eq. 1) then
        psp_JmatType = 0
      else
        psp_JmatType = 2
      endif
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-Jmat_type',psp_JmatType,PETSC_NULL_INTEGER,ierr)

!     default usecolor
      psp_usecolor = 0
      call PetscOptionsHasName(PETSC_NULL_CHARACTER,"-Jmat_coloring",flg,ierr)
      if (flg .or. (psp_JmatType .eq. 2)) then !or Boxed sparse pattern is provided
        psp_usecolor = 1
      endif

      call MatSetType(psp_Jmat,MATAIJ,ierr)
      call MatSetFromOptions(psp_Jmat,ierr)
      i = 9*numvar
      call MatSeqAIJSetPreallocation(psp_Jmat,i,PETSC_NULL_INTEGER,ierr)
      call MatMPIAIJSetPreallocation(psp_Jmat,i,PETSC_NULL_INTEGER,i,PETSC_NULL_INTEGER,ierr)
     
      if (psp_JmatType .eq. 2) then
!         set sparse pattern of Jmat using info from the domain decomposition
!         -------------------------------------------------------------------
          if (mype == 0) then
            zero = 0.0
            do i=1,numvar*numvar*9
               zeros(i) = 0.0
            end do
            do id=1,npes
              do i=1,n,numvar         ! loop over local equations/rows
                row = (id-1)*n + i-1  ! global row index - assume all proc have same num of local rows
                do k=1,numvar
                  rows(k) = row + (k-1)
                end do

                k = 0
                do j=1,9*numvar       ! loop over possible neighbors -tmp!
                  if (ivl2gstnl(i,j,id) .eq. 0) then
                    exit
                  endif
                  cols(j) = ivl2gstnl(i,j,id) - 1
                  k = k+1
                end do !j=1,9*numvar
                call PetscSortInt(k,cols,ierr)

                bj = 1
                do j=1,k,numvar   ! loop over possible neighbors -tmp!
                  bcols(bj) = cols(j)/numvar
                  !write(6,*) "row=",row,cols(j), bcols(bj)
                  !call flush(6)
                  bj = bj+1
                end do !j=1,k
                !k = k/numvar  for MatSetValuesBlocked()
                !write(6,*) "k",k,"bj-1",bj-1
                !call flush(6)

                if (k .gt. 0) then
                  call MatSetValues(psp_Jmat,numvar,rows,k,cols,zeros,INSERT_VALUES,ierr)
                  !call MatSetValuesBlocked(psp_Jmat,1,row/numvar,k,bcols,zeros,INSERT_VALUES,ierr)
                endif

              end do      !i=1,n,numvar
            end do        !id=1,npes
          end if !(mype == 0)
          call MatAssemblyBegin(psp_Jmat,MAT_FINAL_ASSEMBLY,ierr)
          call MatAssemblyEnd(psp_Jmat,MAT_FINAL_ASSEMBLY,ierr)
      end if ! (psp_JmatType .eq. 2)

      return
      end subroutine PetscSnesCreateObjs

! ----------------------------------------------------------------------------------------
!     subroutine AnalyzePostSNESSolve(snes,snesmonitor,nsnessolve,snes_conv_output,iterm)
!
! This is just a function to take the results of a SNES solve and determine what the
! appropriate output message and return value should be.  In reality, all this could
! be in the function PetscSNES, but I like encapsulation.  If you don't like this,
! complain to Mike McCourt, we'll talk.
! ----------------------------------------------------------------------------------------
      subroutine AnalyzePostSNESSolve(snes,snesmonitor,nsnessolve,snes_conv_output,iterm,ierr)
      use Petscmod
      use Snesmonitormod
      use PETSc_Snes_Param
      implicit none

      SNES,           INTENT(IN)    :: snes
      type(monctx),   INTENT(INOUT) :: snesmonitor
      PetscInt,       INTENT(INOUT) :: nsnessolve
      character*80,   INTENT(OUT)   :: snes_conv_output
      PetscInt,       INTENT(OUT)   :: iterm
      PetscErrorCode, INTENT(OUT)   :: ierr

      KSP            :: ksp
      PetscInt       :: snesreason,kspreason
      iterm = 4

      call SNESGetConvergedReason(snes,snesreason,ierr)

      if (snesreason .eq. -3) then
        call SNESGetKSP(snes,ksp,ierr) ! Not called earlier if not in time step 0
        call KSPGetConvergedReason(ksp,kspreason,ierr)
        if (psp_nsnes.eq.0) then ! Don't reevaluate Jacobian at time step 0
          snes_conv_output = "    Linear Solve failed in first SNES solve , reducing time step"
        else
          nsnessolve = nsnessolve - 1
          if (kspreason .eq. -3) then
            if (snesmonitor%JacWasJustRebuilt .eq. 0) then
              snes_conv_output =  "    Linear Solve max iterations exceeded, Jacobian will be rebuilt"
            else
              snes_conv_output =  "    Linear Solve max iterations exceeded, reducing time step"
              nsnessolve = nsnessolve + 1 ! Failure after a rebuild means cut time step
            end if
          elseif (kspreason .eq. -9) then
            snes_conv_output =  "    Linear Solve NAN encountered, Jacobian will be rebuilt"
          else
            snes_conv_output =  "    Unknown Linear Solve failure encountered, Jacobian will be rebuilt"
          end if
        end if
        iterm = 4
        snesmonitor%buildJacobian = 1
      elseif (snesreason .eq. -5) then
        snes_conv_output = "    Max SNES steps exceeded, time step will be reduced"
        iterm = 4
      elseif (snesreason .eq. -6) then
        if (snesmonitor%JacWasJustRebuilt .eq. 0) then
          snes_conv_output = "    Line Search failed, Jacobian will be rebuilt"
          snesmonitor%buildJacobian = 1
          nsnessolve = nsnessolve - 1 ! Try a rebuild before cutting the time step
        else
          snes_conv_output = "    Line Search failed, time step will be reduced"
        end if
        iterm = 6
      elseif (snesreason .eq. -8) then
        if (snesmonitor%JacWasJustRebuilt .eq. 0) then
          snes_conv_output = "    Local Minimum found, Jacobian will be rebuilt"
          snesmonitor%buildJacobian = 1
          nsnessolve = nsnessolve - 1 ! Try a rebuild before cutting the time step
        else
          snes_conv_output = "    Local Minimum found, time step will be reduced"
        end if
        iterm = 7
      else
!           Do not know what value iterm should be set???
        snes_conv_output = "    Unknown SNES failure, time step will be reduced"
        iterm = 4
      endif

      end subroutine AnalyzePostSNESSolve

!************************************************************
!  subroutine PetscSnes()
!  The petscmod module is declared above
!  Include files for petsc routines in fortran
!
!  Input Parameters:
!    n      - num of local equations
!    neqg   - num of global equations
!    yl     - local array for field variable
!    savf   - f(u), function values at u
!    suscal - local scale factors for yl, see Use(Ynorm)
!    sfscal - local scale factors for f, see Use(Ynorm)
!    nnzmx  - maximum no. of nonzeros in Jacobian matrix, see Use(Jacobian)
!    ireorder - flag for reordering Jacobain, see Use(Jacreorder)
!    incpset - maximum num of Jacobian reevaluation, see Use(Lsode)
!    icflag - flag for icnstr
!    icnstr - local interger array. If icflag==1, icnstr is imposed:
!                           icnstr(i)>0 -> yl(i)>0
!                           icnstr(i)<0 -> yl(i)<0
!                           icnstr(i)==0 -> no constraint
!    rlx    - defines the inequality: abs(u_step/u)<rlx
!    dtreal - see Use(Time_dep_nwt)
!  Output Parameter:
!    iterm  - output flag for nksol, see Use(Lsode)
! ----------------------------------------------------------------------------------------------
      subroutine PetscSnes(comm,n,neqg,yl,savf,iscolnorm,suscal,sfscal,nnzmx,ireorder,incpset,icflag,rlx,icnstr,dtreal,mmaxu,iterm)
      use Petscmod
      use Petscmodinterfaces
      use Petscmodinterface2
      use Snesmonitormod
      use PETSc_Snes_Param
      use PETSc_Snes_Param_Monitor

      use Math_problem_size    ! neqmx,numvar
      use Indices_loc_glob_map !ivloc2sdg,ivloc2mdg,ivl2gstnl - rm later!
      use Indices_domain_dcl   !ivloc2sdgl,ivloc2mdgl

      implicit none

      integer4 ::           comm
      PetscInt ::           n,neqg,nnzmx,ireorder,incpset,iscolnorm,mmaxu
      PetscScalar,target :: yl(n+2),savf(n),suscal(n),sfscal(n)
      PetscInt,target ::    icnstr(n)
      PetscInt ::           icflag,jmat_file_val
      PetscReal ::          rlx,dtreal
      PetscInt ::           iterm,isf

      integer4 ::           mype,npes
      SNES ::               snes
      PetscInt ::           ksprestart,snesmaxit,snesfcts
      PetscReal ::          snesatol,snesrtol,snesstol
      PetscInt ::           i,its,j,row,col,ntime
      MatStructure ::       flag
      PetscErrorCode ::     ierr
      KSP ::                ksp
      PC ::                 pc
      Vec ::                X1,Vseq
      PetscScalar,target :: wk(n+2),yldot(n),yl_tmp(n)
      PetscViewer ::        output
      type(Petscctx) ::     pctx

      type(monctx) ::       snesmonitor
      PetscTruth ::         writeJmat,flg
      SNESConvergedReason :: snesreason
      KSPConvergedReason :: kspreason

      ISColoring ::         iscoloring

      PetscScalar ::        neg_one,zero
      PetscReal ::          norm,dtholder,cmp_rhs,delta
      PetscScalar, pointer ::  xx(:),xx1(:),array(:)
      real ::               gettime, sec4
      PetscViewer ::        viewer
      MatStructure ::       Jstruct

      PetscInt ::           id
      PetscLogDouble ::     time_loc(30),tstart,tend
      PetscInt ::           idx(3000),nsnessolve,nsnessolve_max,nsnesiter
      character*4           npes_str
      character*80          snes_conv_output
      character*10          snes_jmat_file
      PetscReal, ALLOCATABLE, DIMENSION(:) :: time_diff

!  All pctx-defined Fortran routines are declared external
      external FormUedgeFunction,FormUedgeJacobian,FormPreCheck
      external PetscSnesSetSf,PetscSnesFillSf,PetscSnesFillSu
      external FormMonitor
      external rhsnk

      Interface FormInitialGuess
        subroutine FormInitialGuess(neq,yl,X,ierr)
        use Ynorm                    ! suscal,sfscal,iscolnorm
        use Petscmod
        use PETSc_Snes_Param

        implicit none
        PetscInt ::           neq
        PetscScalar,target :: yl(neq)
        Vec ::                X
        PetscErrorCode ::     ierr
        end subroutine
      End Interface
      Interface FormFinalSolution
        subroutine FormFinalSolution(n,X,yl,ierr)
        use Ynorm                    ! suscal,sfscal,iscolnorm
        use Petscmod
        use PETSc_Snes_Param
        implicit none
        PetscInt ::           n
        PetscScalar,target :: yl(n)
        Vec ::                X
        PetscErrorCode ::     ierr
        end subroutine
      End Interface
      Interface  RhsCmp_Debug
        subroutine RhsCmp_Debug(comm,cmp_rhs,snes,pctx,n,yl)
        use Petscmod
        use PETSc_Snes_Param
        use Indices_loc_glob_map !ivloc2sdg,ivloc2mdg,ivl2gstnl - rm later!
        use Indices_domain_dcl   !ivloc2sdgl,ivloc2mdgl
        implicit none
        integer4 ::           comm
        PetscReal ::          cmp_rhs
        SNES ::               snes
        type(Petscctx) ::     pctx
        PetscScalar,target :: yl(n)
        PetscErrorCode ::     ierr
        PetscInt ::           n
        end subroutine
      End Interface

      call MPI_Barrier(comm,ierr)
      call PetscGetTime(time_loc(1),ierr)

      call MPI_Comm_size(comm,npes,ierr)
      call MPI_Comm_rank(comm,mype,ierr)
!      write(6,*) "comm = ", comm, ", npes = ", npes, ", mype = ", mype
!      call flush(6)

      if (npes .eq. 1) neqg = n
      if (mype .eq. 0) then
        write(6,*)"*****************************************************************"
        write(6,1000)psp_nsnes,n,neqg,numvar,dtreal
1000    format(I2,"-th PetscSnes: neq neqg=",I6,I7,",  numvar=",I2,", dtreal=",EN10.2)
        write(6,*)"................................................................."
        call flush(6)
      endif

      if (psp_nsnes .eq. 0) then !************************************

!     Create snes, vectors and Jacobian matrix
!     Note that iscolnorm is needed in PetscSnesCreateObjs
!---------------------------------------------------
        psp_iscolnorm = iscolnorm
        call PetscLogStagePush(stages(1),ierr)
        call PetscSnesCreateObjs(comm,n,neqg,nnzmx,ierr)
        call PetscLogStagePop(ierr)

!       Fill scaling vectors
!----------------------------------------------
!        call PetscLogStagePush(stages(2),ierr)
        if (psp_iscolnorm .gt. 0) then
          call PetscSnesFillSu(n,suscal,psp_SUVEC,psp_ISUVEC,ierr)
        endif

!       sfscal is provided by Uedge when using Uedge Jacobian
        if ((psp_isf.eq.1) .and. (psp_JmatType.eq.0)) then
          call PetscSnesFillSf(n,sfscal,psp_SFVEC,ierr)
        endif

!       Form initial guess
!------------------------------
        call FormInitialGuess(n,yl,psp_X,ierr)

!       Setup SNES context
!----------------------------------------------------------
        call SNESGetKSP(psp_snes,ksp,ierr)
        ksprestart  = mmaxu
        snesatol    = 1.e-10
        snesrtol    = 0.0
        snesstol    = 0.0
        snesmaxit   = 40
        snesfcts    = 10000
        call KSPGMRESSetRestart(ksp,ksprestart,ierr)
        call SNESSetTolerances(psp_snes,snesatol,snesrtol,snesstol,snesmaxit,snesfcts,ierr)

!  Determine if we need to use a reordering on the Jacobian
!---------------------------------------------------------------
        if (ireorder.eq.1) then
          call KSPGetPC(ksp,pc,ierr)
          call PCFactorSetMatOrderingType(pc,MATORDERING_RCM,ierr)
        endif

!  Determine what the Jacobian lag should be
!  This is passed with the option -snes_lag_jac <int>
!  Note this is different than the standard Jacobian option -snes_lag_lacobian,
!    which prevent PETSc from throwing an error that we already handle here
!  Default is likely 5, which is set in com/petsc-mod.F90
!-------------------------------------------------------------------
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-snes_lag_jac',psp_jaclag,flg,ierr)
      if (psp_jaclag .le. 0) then
        if (mype .eq. 0) then
1006  format("**NOTE** Nonsensical Jacobian lag, ",I5," passed, being reset to ",I2)
          write(6,1006)psp_jaclag,psp_jaclag_DEFAULT
        end if
        psp_jaclag = psp_jaclag_DEFAULT
      end if

!  Prepare the stats if needed 
!  To activate the stats collection and write to a file <str>, pass -stats_print <str>
!  The psp_mon_file default is likely 'trialuedge.m', defined in com/petcsMod.F90
!
!  To write the SNESView() from the final SNES solve to <str>, pass -snes_print <str>
!  The psp_mon_sfile default is likely 'snesout.txt', defined in com/petscMod.F90
!
!  To change the size of the dynamic memory allocation needed for statistics, pass -dynamic_alloc_increment <int> 
!  The psp_mon_allocInc default is likely 10, defined in com/petscMod.F90
!--------------------------------------------------------
        call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-dyanmic_alloc_increment',psp_mon_allocInc,flg,ierr)

        psp_mon_on = PETSC_FALSE
        call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-stats_print',psp_mon_on,ierr)
        if (psp_mon_on) then
          if (psp_mon_allocInc .lt. 1) then
            psp_mon_allocInc = psp_mon_allocInc_DEFAULT ! Found in com/petscMod.F90
          end if
          call psp_monAdjustSize(psp_mon_allocInc,ierr)

          call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-stats_print',psp_mon_file,psp_mon_on,ierr)
          if (psp_mon_file .eq. '') then
            psp_mon_file = psp_mon_file_DEFAULT
          end if
          psp_mon_on = PETSC_TRUE ! This may not be needed, but was once so I'm leaving it in
        end if

        psp_mon_print_snes = PETSC_FALSE
        call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-snes_print',psp_mon_print_snes,ierr)
        if (psp_mon_print_snes) then
          call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-snes_print',psp_mon_sfile,psp_mon_print_snes,ierr)
          if (psp_mon_sfile .eq. '') then
            psp_mon_sfile = psp_mon_sfile_DEFAULT ! Found in com/petscMod.F90
          end if
          psp_mon_print_snes = PETSC_TRUE ! In case the user wants to use the default
        end if

      endif              !************************************

      snes = psp_snes
     
      call SNESSetApplicationContext(snes,pctx,ierr)

!  Moving the data to the workspace for FormUedgeJacobian
!--------------------------------------------------------
      pctx%savfpoint=>savf
      pctx%wkpoint=>wk
      pctx%ylpoint=>yl
      pctx%nnzmx   = nnzmx
      pctx%incpset = incpset
      pctx%rlx     = rlx
      pctx%icnstrpoint=>icnstr
!      call MPI_Barrier(comm,ierr)

!     Set the rhs function 
!-----------------------------------------
      call SNESSetFunction(snes,psp_F,FormUedgeFunction,pctx,ierr)

!     Set the monitor object
!----------------------------------------
      call SNESMonitorSet(snes,FormMonitor,snesmonitor,PETSC_NULL_FUNCTION,ierr)

!     Create matrix-free matrix Jmf used for Jacobian*vector
!-----------------------------------------------------------
      if (psp_nsnes .eq. 0) then
        call MatCreateSNESMF(snes,psp_Jmf,ierr)
        call MatMFFDSetFromOptions(psp_Jmf,ierr)
      end if

!     Set the function to compute Jacobian
!-----------------------------------------
      if (psp_JmatType .eq. 0) then
        if ((mype .eq. 0) .and. (psp_nsnes .eq. 0)) write(6,*) " Use Uedge jacobian..."
        call SNESSetJacobian(snes,psp_Jmf,psp_Jmat,FormUedgeJacobian,pctx,ierr)
      else if (psp_JmatType .eq. 1) then
        if ((mype .eq. 0) .and. (psp_nsnes .eq. 0)) write(6,*) " Use PETSc slow fd jacobian ..."
        call SNESSetJacobian(snes,psp_Jmf,psp_Jmat,SNESDefaultComputeJacobian,pctx,ierr)
      else
        if ((mype .eq. 0) .and. (psp_nsnes .eq. 0)) write(6,*) " Use Uedge sparse patten with PETSc coloring ..."
      endif
      call flush(6)

      if (psp_nsnes .eq. 0) then
      end if

!  Handle changes in the time step if necessary
!----------------------------------------------
!  In call set_dt(n,yl,wk) we incur a function evaluation which gets
!  stored in wk, but it's not needed which is why wk is just a workspace
!
!  Also, these values get used in FormMonitor to determine when a
!  Jacobian update needs to occur because of a time step change

      snesmonitor%thisdt = dtreal
      if (snesmonitor%thisdt .ne. psp_olddt) then
        call set_dt(n,yl,wk)
      end if

!  Debugging: compare Frhs at the beginning of ntime time steps
!--------------------------------------------------------------
      ntime = 0
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-cmp_Frhs_ntime',ntime,flg,ierr)
      cmp_rhs = 1.e-16
      flg     = PETSC_FALSE
      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-cmp_Frhs',cmp_rhs,flg,ierr)
      if (flg .and. (psp_nsnes .eq. ntime)) then !************************************
        !write(6,*) "nsnes",psp_nsnes," ntime",ntime
        !call flush(6)
        call RhsCmp_Debug(comm,cmp_rhs,snes,pctx,n,yl)
      endif

!   Create coloring context. Set the function to evaluate Jacobian
!----------------------------------------------------------------
      if (psp_usecolor .eq. 1) then
        if ((mype .eq. 0) .and. (psp_nsnes .eq. 0)) then
          write(6,*) " use coloring, JmatType", psp_JmatType," isf",psp_isf
          call flush(6)
        endif
        if (psp_nsnes .eq. 0) then !************************************
          call PetscGetTime(tstart,ierr)
          if (psp_JmatType .lt. 2) then
!           Assemble Jmat
            isf = psp_isf
            psp_isf = 0 !if psp_isf=1, set it as 0 because sf is not available until Jacobian is ready
            call SNESComputeJacobian(snes,psp_X,psp_Jmat,psp_Jmat,Jstruct,ierr)
            psp_isf = isf
          endif

          call MatGetColoring(psp_Jmat,MATCOLORING_SL,iscoloring,ierr)
          call MatFDColoringCreate(psp_Jmat,iscoloring,psp_fdcoloring,ierr)
          call ISColoringDestroy(iscoloring,ierr)
          call MatFDColoringSetFromOptions(psp_fdcoloring,ierr)
          call PetscGetTime(tend,ierr)
          psp_mon_colortime = tend - tstart
        endif                      !************************************

        call MatFDColoringSetFunction(psp_fdcoloring,FormUedgeFunction,pctx,ierr)
        call SNESSetJacobian(snes,psp_Jmf,psp_Jmat,SNESDefaultComputeJacobianColor,psp_fdcoloring,ierr)
      endif

!      write(6,*) "psp_usecolor = 1 block concluded. npes = ", npes
!      call flush(6)

!  Enable runtime options
!---------------------------------------------------------
      if (psp_nsnes .eq. 0) then !************************************
        call SNESSetFromOptions(snes,ierr)  ! must be called prior to SNESLineSearchSetPreCheck()
      endif                      !************************************

      if (icflag .eq. 1) then !icflag=1:  flag to use constraint that ni, etc. not < 0
        call SNESLineSearchSetPreCheck(snes,FormPreCheck,pctx,ierr)
      endif

!  Set sfscal (sfscal is not provided by uedge when psp_Jmat_type>0)
!--------------------------------------------------------------------
      if ((psp_JmatType .gt. 0) .and. (psp_nsnes .eq. 0) .and. (psp_isf .eq. 1)) then !************************************
        psp_isf = 0
        call PetscSnesSetSf(snes,psp_X,psp_SFVEC,ierr)
        psp_isf = 1
      
#if defined(SetSf_Debugging)
        write(6,*)" sfscal, SFVEC:"
        write(6,*)" -----------------"
        call VecGetArrayF90(psp_SFVEC,xx,ierr)
        do i=1,n
          write(6,*) i,sfscal(i),xx(i),(sfscal(i) - xx(i))
        end do
        call flush(6)
        call VecRestoreArrayF90(psp_SFVEC,xx,ierr)
#endif
      endif                                                    !************************************

!  Determine whether the first Jacobian needs to be skipped
!  Get frequency of rebuilding Jacobian; then set it to the monitoring parameter
!  Also store the max time step change allowed before a Jac reevaluation is needed
!
!    Note: Passing -jac_max_dt_change 0 means never reevaluate for time step reasons
!----------------------------------------------------------------------------
      psp_skipJacobian  = 0 ! default
      snesmonitor%buildJacobian = 0 ! default
      call PetscOptionsHasName(PETSC_NULL_CHARACTER,"-skip_1st_jacobian",flg,ierr)
      if (flg) then
        psp_skipJacobian = 1
      endif
      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,"-jac_max_dt_change",dtholder,flg,ierr)
      if (flg) then
        if (dtholder .lt. 0) then
          dtholder = -dtholder
        end if
        if ((dtholder .lt. 1) .and. (dtholder .ne. 0)) then
          dtholder = 1/dtholder
        end if
        psp_jacMaxdtChange = dtholder+.001 ! Allowing for some roundoff
      else
        psp_jacMaxdtChange = 1.001 ! Allowing for some roundoff
      end if
!      call PetscLogStagePop(ierr)

!  Solve the system
!  Note that snesmonitor%buildJacobian = 1 for a linear solve fail
!  Other SNES failures also incur a Jac reevaluation, but that is handled in FormMonitor
!---------------------------------------------------------
      call PetscLogStagePush(stages(3),ierr)
      nsnessolve = 0
      do while (nsnessolve .le. 1)
        nsnessolve = nsnessolve + 2
        call SNESSolve(snes,PETSC_NULL_OBJECT,psp_X,ierr)

        call SNESGetConvergedReason(snes,snesreason,ierr)
        if (snesreason .gt. 0) then
          iterm = 1
!          exit
        else ! not converge
          call AnalyzePostSNESSolve(snes,snesmonitor,nsnessolve,snes_conv_output,iterm,ierr)
          if (mype .eq. 0) then
            write(6,*) snes_conv_output
            call flush(6)
          endif
        endif
      end do
      call PetscLogStagePop(ierr)

!  Recover X with scaling suscal and convert it back to Uedge yl
!---------------------------------------------------------------
      call PetscLogStagePush(stages(4),ierr)
!      write(6,*) "Calling FormFinalSolution."
!      call flush(6)
      call FormFinalSolution(n,psp_X,yl,ierr)

!      write(6,*) "Calling SNESGetIterationNumber."
!      call flush(6)
      call SNESGetIterationNumber(snes,its,ierr);
      call SNESGetFunctionNorm(snes,norm,ierr);
      if (mype .eq. 0) then
        write(6,*)"SNES its ",its," fnorm ",norm," convergedReason ",snesreason, "iterm ",iterm
        call flush(6)
      endif

!  Store the current time step to compare to the next time step
!  This value will be needed to determine if the Jacobian needs updated
!----------------------------------------------------------------------
      psp_olddt = snesmonitor%thisdt

!  Write the Jacobian matrix at time ntime
!
!  You may have to manually make a directory jacobians before you use this option the first time
!  That directory must be in the same location as the driver script
!
!  The format for the output of the Jacobian is
!    Jv<int>c<int>n<int>.s<int>
!      v : number of variables (1-6)
!      c : coloring choice (0-3)
!          0 : UEDGE Jacobian
!          1 : Dense Jacobian Only
!          2 : Tom's Sparse Boxed Pattern
!          3 : First Dense Evaluation and Future PETSc Colorings
!      n : number of processors (1,2,4,8)
!      s : snes solve number 
!  If you ran in parallel there would also be a permutation file:
!    Pv<int>c<int>n<int>.s<int>
!  So for example, Jv4c2n1.s0 has 4 variables, Tom's coloring, 1 processor and the first Jacobian
!
!  The s term is the term passed in uedgepetscopt with -write_Jmat <int>
!  The default value is zero, so -write_Jmat 0 is equivalent to -write_Jmat
!  If you pass a value which never gets reached, then no Jacobian is printed
!  It's tempting to call this term the time step, but it isn't because the time step can be split
!--------------------------------------------------------------------------
      writeJmat = PETSC_FALSE
      call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-write_Jmat',writeJmat,ierr)
      if (writeJmat) then
        ntime = 0
        call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-write_Jmat',ntime,writeJmat,ierr)
        writeJmat = PETSC_TRUE ! May be set to false if no integer is passed
      end if
      if (writeJmat .and. (psp_nsnes .eq. ntime)) then
        if (psp_JmatType .eq. 1) then
          if (psp_usecolor .eq. 1) then
            jmat_file_val = 3
          else
            jmat_file_val = 1
          endif
        elseif (psp_JmatType .eq.2) then
          jmat_file_val = 2 ! This is the boxed sparse pattern
        else
          jmat_file_val = 0 ! I guess we could still run with Uedge Jacobian
        endif
 1002   format('v',I1,'c',I1,'n',I1,'.s',I1)
 1003   format('v',I1,'c',I1,'n',I1,'.s',I2)
        if (ntime .le. 9) then
          write(snes_jmat_file,1002) numvar,jmat_file_val,npes,ntime
        else
          write(snes_jmat_file,1003) numvar,jmat_file_val,npes,ntime
        endif
        if (mype .eq. 0) then
          write(6,*)' write Jmat in binary to jacobians/J'//snes_jmat_file
          call flush(6)
        endif

        call PetscViewerBinaryOpen(comm,'jacobians/J'//snes_jmat_file,FILE_MODE_WRITE,viewer,ierr)
        call MatView(psp_Jmat,viewer,ierr)
        call PetscViewerDestroy(viewer,ierr)

        if ((mype .eq. 0) .and. (npes .gt. 1)) then ! Write Permutation to file
          write(6,*)" write permutation ivloc2sdg to jacobians/P"//snes_jmat_file
          call flush(6)
          call VecCreateSeq(MPI_COMM_SELF,neqg,Vseq,ierr)
          call VecGetArrayF90(Vseq,array,ierr)
          do i=1,neqg
            id = (i-1)/n
            j  = i - id*n
            id = id+1
            array(i) = ivloc2sdg(j,id)
          end do
          call VecRestoreArrayF90(Vseq,array,ierr)

          call PetscViewerBinaryOpen(MPI_COMM_SELF,"jacobians/P"//snes_jmat_file,FILE_MODE_WRITE,viewer,ierr)
          call VecView(Vseq,viewer,ierr)
          call PetscViewerDestroy(viewer,ierr)
          call VecDestroy(Vseq,ierr)
        end if

      endif

      if (mype .eq. 0) then
        call PetscGetTime(time_loc(2),ierr)
        write(6,2000) psp_nsnes,time_loc(2)-time_loc(1)
2000    format("*****",i3,"-th PetscSnes is done, Time:",f6.2,"  *****")
        call flush(6)
      endif
      psp_nsnes = psp_nsnes + 1

      call PetscLogStagePop(ierr)

      return
      end subroutine PetscSnes

! ---------------------------------------------------------------------------------
!  subroutine outputstatspetsc
!
!  The info for this function is activated at runtime by -print_stats <int>
!  These statistics will be output to the file snesout.m
!  Eventually I'll add an output to say what file you want the stats written to
!  All the values will be updated in the function FormMonitor
!  Note that the default value for psp_mon_allocInc is 10, defined in com/petscMod.F90
!
!  Since we sometimes do a final run of dt 1e-9 to output the variable values to an hd5
!  file, this function sets psp_mon_on=PETSC_FALSE at the end.  This prevents SNES from
!  trying to write to arrays that are no longer allocated.  If this needs to be
!  changed at some points, we'll figure it out, but for right now this seems like the
!  easiest way to make sure the memory is freed.
! ----------------------------------------------------------------------------------
      subroutine outputstatspetsc
      use Petscmod
      use PETSc_Snes_Param
      use PETSc_Snes_Param_Monitor
      use Npes_mpi         ! npes,mype,ismpion
      use Math_problem_size    ! neqmx,numvar
      implicit none

      PetscInt       :: i,n
      PetscReal      :: time0,time00,ftime0,ftime00
      PetscErrorCode :: ierr

      PetscReal      :: atol,rtol,stol,alpha,maxstep
      PetscInt       :: maxit,maxf
      SNES           :: snes
      KSP            :: ksp
      PetscViewer    :: snesview_output
      character*10   :: name

      snes = psp_snes
      call VecGetSize(psp_X,n,ierr)

      if (.not. psp_mon_on) then ! Can't do anything if you haven't collected any stats
        if (mype .eq. 0) then
          write(6,*) "bbb.outputstats() called but -stats_print <str> not activated"
          write(6,*) "If -stats_print <str> is not activated, data is not collected during the solve"
          write(6,*) "The <str> value gives the file you want the statistics written to"
          write(6,*) "The default value of ",psp_mon_file," will be actived if -stats_print is passed"
          call flush(6)
        end if
      else if (mype .eq. 0) then ! Only have one core do the output
        open(unit=5,file=psp_mon_file,action="write",status="replace")
        write(5,*) "% This data was output from UEDGE"
        write(5,*) "%"
        write(5,*) "pspm_np=",npes,"; % number of processors"
        write(5,*) "pspm_numvar=",numvar,"; % number of variables"
        write(5,*) "pspm_neq=",n,"; % number of equations"
        write(5,*) "%"

        call SNESGetTolerances(snes,atol,rtol,stol,maxit,maxf,ierr)
        write(5,*) "% SNES tolerances and such"
        write(5,*) "pspm_atol=",atol,"; % SNES absolute convergence"
        write(5,*) "pspm_rtol=",rtol,"; % SNES relative convergence"
        write(5,*) "pspm_stol=",stol,"; % SNES step change convergence"
        write(5,*) "pspm_maxit=",maxit,"; % max SNES its before divergence"
        write(5,*) "pspm_maxf=",maxf,"; % max f evals before divergence"
        write(5,*) "%"

        call SNESLineSearchGetParams(snes,alpha,maxstep,ierr)
        call SNESGetType(snes,name,ierr)
        write(5,*) "% SNESLS tolerances and such"
        write(5,*) "pspm_snestype='",name,"'; % Should always be Line Search"
        write(5,*) "pspm_alpha=",alpha,"; % minimum LS reduction factor"
        write(5,*) "pspm_maxstep=",maxstep,"; % maximum norm of Newton step"
        write(5,*) "%"

        call SNESGetKSP(snes,ksp,ierr)
        call KSPGetType(ksp,name,ierr)
        call KSPGetTolerances(ksp,rtol,atol,stol,maxit,ierr)
        write(5,*) "% KSP tolerances and such"
        write(5,*) "pspm_ksp_type='",name,"'; % Should be GMRES"
        write(5,*) "pspm_ksp_rtol=",rtol,"; % KSP relative convergence"
        write(5,*) "pspm_ksp_atol=",atol,"; % KSP absolute convergence"
        write(5,*) "pspm_ksp_dtol=",stol,"; % KSP divergence"
        write(5,*) "pspm_ksp_maxit=",maxit,"; % max KSP its before divergence"
        write(5,*) "%"

        write(5,*) "% Jacobian parameters (specified from uedgepetscopt)"
        write(5,*) "pspm_jac_lag=",psp_jaclag,"; % Jacobian rebuild lag"
        write(5,*) "pspm_jac_type=",psp_JmatType,"; % Method for computing Jacobian"
        write(5,*) "pspm_jac_color=",psp_usecolor,"; % Jacobian coloring activated"
        write(5,*) "pspm_jac_skip1=",psp_skipJacobian,"; % Don't always reevaluated at 1st SNES step"
        write(5,*) "pspm_jac_maxdt=",psp_jacMaxDtChange,"; % dt factor for reevaluating Jacobian"
        write(5,*) "%"

! The last value never gets used, so don't count it
        psp_mon_arrlen = psp_mon_arrlen - 1

! Output the time required to create the coloring, which could be significant for -Jmat_type 1
        write(5,*) "pspm_colortime=",psp_mon_colortime,"; % Time for creating the coloring context"

        write(5,*) "% Col  1: index"
        write(5,*) "% Col  2: Time step size (dtreal)"
        write(5,*) "% Col  3: SNES step"
        write(5,*) "% Col  4: Function norm (residual)"
        write(5,*) "% Col  5: KSP Its in this SNES solve"
        write(5,*) "% Col  6: Cumulative function evals in the time step"
        write(5,*) "% Col  7: Calls to FormUedgeFunction this SNES step"
        write(5,*) "% Col  8: Time spent in FormUedgeFunction this SNES step"
        write(5,*) "% Col  9: Cumulative time in the time step"
        write(5,*) "% Col 10: Cumulative time overall"
        write(5,*) "uedgepetscData=["
        time00 = psp_mon_time(1)
        ftime00 = psp_mon_ftime(1)
1004  format("% Col",I3,", ",I8  ,", ",I2,", ",I10  ,", ",I5,", ",I6,", ",I7,", ",I9  ,", ",I9  ,", ",I9)
1005  format("     ",I3,", ",E8.2,", ",I2,", ",E10.3,", ",I5,", ",I6,", ",I7,", ",F9.3,", ",F9.3,", ",F9.3)
        write(5,1004) 1,2,3,4,5,6,7,8,9,10
        do i=1,psp_mon_arrlen
          if (psp_mon_snesit(i) .eq. 0) then
            time0 = psp_mon_time(i)
            ftime0 = psp_mon_ftime(i)
          end if
          write(5,1005) i,psp_mon_dtreal(i),psp_mon_snesit(i),psp_mon_fnrm(i),psp_mon_kspits(i), &
                        psp_mon_fevals(i),psp_mon_fcalls(i),psp_mon_ftime(i),psp_mon_time(i)-time0+ftime0, &
                        psp_mon_time(i)-time00+ftime00
        end do      
        write(5,*) "];"

        call flush(5)
        close(unit=5)

        write(6,*) "Writing stats to ",psp_mon_file
        call flush(6)
        call psp_monAdjustSize(int(0,4),ierr)
        psp_mon_on = PETSC_FALSE ! So that hd5 runs don't write to psp_mon arrays
      end if

      if (psp_mon_print_snes) then
        if (mype .eq. 0) then
          write(6,*) "Writing SNESView() to ",psp_mon_sfile
          call flush(6)
        end if
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD,psp_mon_sfile,snesview_output,ierr)
        call SNESView(snes,snesview_output,ierr)
        call PetscViewerDestroy(snesview_output,ierr)
      end if

      end subroutine outputstatspetsc
