MODULE PETSc_Snes_Param
  implicit none
#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
  SAVE
  PetscInt :: psp_nsnes=0
  Vec :: psp_snes=0
  Mat :: psp_Jmat=0
  Mat :: psp_Jmf=0
  PetscInt :: psp_JmatType=1
  PetscInt :: psp_fdcoloring=0
  PetscInt :: psp_usecolor=0
  Vec :: psp_F=0
  Vec :: psp_Floc=0
  Vec :: psp_X=0
  Vec :: psp_Xtmp=0
  Vec :: psp_Xloc=0
  Vec :: psp_SFVEC=0
  Vec :: psp_SUVEC=0
  Vec :: psp_ISUVEC=0
  PetscInt :: psp_isf=1
  PetscInt :: psp_ijnr=0
  PetscInt :: psp_iscolnorm=0
  PetscInt :: psp_skipJacobian=0
  PetscReal :: psp_jacMaxdtChange
  PetscReal :: psp_olddt=0.0
  PetscInt :: snesdebug=0 ! Does this do anything?
  PetscInt, PARAMETER :: psp_snesits=30
  PetscReal, DIMENSION(psp_snesits+1) :: nksolfnrm
  PetscReal, DIMENSION(psp_snesits+1) :: nksoltime
  PetscInt,  DIMENSION(psp_snesits+1) :: nksollinits
  PetscInt,  DIMENSION(psp_snesits+1) :: nksolfeval

  PetscInt, PARAMETER :: psp_jaclag_DEFAULT = 5
  PetscInt            :: psp_jaclag         = psp_jaclag_DEFAULT
END MODULE PETSc_Snes_Param

MODULE PETSc_Snes_Param_Monitor
  implicit none
#include "finclude/petscdef.h"
  SAVE
  PetscTruth :: psp_mon_on
  PetscTruth :: psp_mon_print_snes
  PetscInt   :: psp_mon_arrlen=1

  character*30, PARAMETER :: psp_mon_file_DEFAULT     = "trialuedge.m"
  character*30, PARAMETER :: psp_mon_sfile_DEFAULT    = "snesout.txt"
  PetscInt,     PARAMETER :: psp_mon_allocInc_DEFAULT = 10
  character*30            :: psp_mon_file             = psp_mon_file_DEFAULT
  character*30            :: psp_mon_sfile            = psp_mon_sfile_DEFAULT
  PetscInt                :: psp_mon_allocInc         = psp_mon_allocInc_DEFAULT

  PetscLogDouble :: psp_mon_colortime

  PetscReal,      ALLOCATABLE, DIMENSION(:) :: psp_mon_fnrm
  PetscReal,      ALLOCATABLE, DIMENSION(:) :: psp_mon_dtreal
  PetscLogDouble, ALLOCATABLE, DIMENSION(:) :: psp_mon_time
  PetscLogDouble, ALLOCATABLE, DIMENSION(:) :: psp_mon_ftime
  PetscInt,       ALLOCATABLE, DIMENSION(:) :: psp_mon_kspits
  PetscInt,       ALLOCATABLE, DIMENSION(:) :: psp_mon_fevals
  PetscInt,       ALLOCATABLE, DIMENSION(:) :: psp_mon_fcalls
  PetscInt,       ALLOCATABLE, DIMENSION(:) :: psp_mon_snesit
END MODULE PETSc_Snes_Param_Monitor
