*  ---------------------------------------------------------------------
* @file: dvshdf5.m
*
* @purpose: dump data
*
* @version $Id: dvshdf5.m,v 7.0 2018/02/28 18:32:53 meyer8 Exp $
*
*  ---------------------------------------------------------------------

*  ---------------------------------------------------------------------
*     Initialize hdf5 by calling initialization methods in vshdf5.
*  ---------------------------------------------------------------------
      subroutine init_hdf5()
      Use(hdf5_api)      # hdf5 types
      CALL vshdf5_fcinit()
      return
      end

*  ---------------------------------------------------------------------
*     Setup and initialization
*     For parallel runs, we are only going to do it on node 0.
*  ---------------------------------------------------------------------
      subroutine dump_vshdf5(filename,intime,doTransposeFlag)
      Use(Dim)           # nisp,nusp,nfsp,ngsp, nhgsp
      Use(Interp)        # nis,ups,tes,tis,ngs,phis
      Use(Jacreorder)    # ireorder
      Use(Jacobian)      # nnzmx
      Use(Xpoint_indices)  # ixpt2
      Use(RZ_grid_info)  # rm, zm
      Use(hdf5_api)      # nnzmx
      Use(Neqdskv)       # zmagx
      Use(UEint)         # mhdgeo
      Use(Share)         # nxleg,nxcore
c_mpi      Use(MpiVars)
c_mpi      Use(Npes_mpi)
      Use(Conduc)        # dif_use, tra_use, ...
      Use(Comtra)        # isbohmcalc facbni, facbup, facbni2, facbee
      implicit none
      character*(*), intent(in) :: filename
      integer, intent(in) :: doTransposeFlag
      real, intent(in) :: intime
      TYPE(hdf5ErrorType) :: h5err
      TYPE(hdf5InOpts) :: h5in
      character(60) :: description
      real, dimension(:,:,:), allocatable :: vs3d
      real, dimension(:,:), allocatable :: vs2d
      integer i,ifld,lid, lenx, leny, ixs, ixe, iys, iye
      real(Size4) gettime, sec4
      INTEGER(HID_T) :: root_gid,bbbid,fid, runInfoId
      logical :: dumpPriv, dumpCore, dumpSOL, dumpPrivRite
      logical :: dumpPrivLeft,dumpCoreLeftCell,dumpCoreRiteCell
      integer ::  lxs, lxe, lys, lye, rxs, rxe, rys, rye
*  ---------------------------------------------------------------------
*     Setup and initialization
*     For parallel runs, we are only going to do it on node 0.
*  ---------------------------------------------------------------------
      ! write(*,*) "dump_vshdf5: mype = ", mype
      IF(mype/=0) THEN
        ! write(*,*) "dump_vshdf5: returning as mype = ", mype
        return
      ENDIF
      call vshdf5_inith5vars(h5in, h5err)
      h5in%wrd_type=H5T_NATIVE_DOUBLE
      h5in%doTranspose=.true.
      h5in%vsCentering="zonal"
      h5in%vsTime=intime
      h5in%wrVsTime=.true.
      h5in%verbose=.false.
      h5in%vstimegroup = "ue_time"
c      h5in%logUnit=99   ! Need to fix

c The below are done with -Dc__MPI= , which is the case for truly parallel
c__MPI      h5in%comm=MPI_COMM_SELF
c__MPI      h5in%info=MPI_INFO_NULL

*  ---------------------------------------------------------------------
*     Write the variables to an hdf5 file
*  ---------------------------------------------------------------------
      description=" "
      write(*,*) "dump_vshdf5: opening ", filename,
     &    ", with h5in%comm = ", h5in%comm
c JRC, xlf:(E) Null literal string is not permitted.  A single blank is assumed.
      CALL open_newh5file(filename,fid,description,root_gid,h5in,h5err)
      IF(h5err%errBool) THEN
        WRITE(*, *) "dump_vshdf5: open_newh5file failed with '",
     &  h5err%errorMsg, "'"
      ENDIF
      ! Put time group at the top level.  The references to this group
      ! get added automatically
      call make_time_group(root_gid,h5in,h5err)
      ! Put everything into the bbb subdirectory
      CALL make_group(root_gid,"bbb",bbbid,h5in,h5err)
      IF(h5err%errBool) THEN
        WRITE(*, *) "dump_vshdf5: make_group for bbb failed with '",
     &  h5err%errorMsg, "'"
      ENDIF

*  ---------------------------------------------------------------------
*     Write the run information
*  ---------------------------------------------------------------------
      CALL make_group(root_gid, "runInfo", runInfoId, h5in, h5err)
      IF(h5err%errBool) THEN
        WRITE(*, *) "dump_vshdf5: make_group for runInfo failed with '",
     &  h5err%errorMsg, "'"
      ENDIF
      CALL write_attribute(runInfoId, "vsSoftware", "UEDGE", h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
      CALL write_attribute(runInfoId, "vsVsSoftware", "2.9.9", h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
      CALL close_group("runInfo", runInfoId, h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg

*  ---------------------------------------------------------------------
*     Basic dumping of raw arrays
*  ---------------------------------------------------------------------
      h5in%units="m^-3"
      CALL dump_h5(bbbid,"nis",nis,h5in,h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
      h5in%units="m^-3"
      CALL dump_h5(bbbid,"ngs",ngs,h5in,h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
      h5in%units=" "
c JRC, xlf:(E) Null literal string is not permitted.  A single blank is assumed.
      CALL dump_h5(bbbid,"tes",tes,h5in,h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
      CALL dump_h5(bbbid,"tis",tis,h5in,h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
      CALL dump_h5(bbbid,"ups",ups,h5in,h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
      CALL dump_h5(bbbid,"psi",psi,h5in,h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
      CALL dump_h5(bbbid,"phis",phis,h5in,h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
      ! diffusivities if need !-------------------------------------
      ! It is hard to tell what to use and what not to use
      IF(read_diffs==1) THEN
      h5in%units="m^2/s"
      CALL dump_h5(bbbid,"dif_use", dif_use,h5in,h5err)
      CALL dump_h5(bbbid,"kye_use", kye_use,h5in,h5err)
      CALL dump_h5(bbbid,"kyi_use", kyi_use,h5in,h5err)
!      CALL dump_h5(bbbid,"tra_use", tra_use,h5in,h5err)
!      CALL dump_h5(bbbid,"dutm_use",dutm_use,h5in,h5err)
!      CALL dump_h5(bbbid,"vy_use",  vy_use,h5in,h5err)
!      CALL dump_h5(bbbid,"vyup_use",vyup_use,h5in,h5err)
!      CALL dump_h5(bbbid,"vyte_use",vyte_use,h5in,h5err)
!      CALL dump_h5(bbbid,"vyti_use",vyti_use,h5in,h5err)
!      CALL dump_h5(bbbid,"fniyos_use", fniyos_use,h5in,h5err)
!      CALL dump_h5(bbbid,"feeyosn_use",feeyosn_use,h5in,h5err)
!      CALL dump_h5(bbbid,"feiyosn_use",feiyosn_use,h5in,h5err)
      ENDIF

*  ---------------------------------------------------------------------
*  Figure out what to dump
*  ---------------------------------------------------------------------
      dumpPriv=.true.; dumpCore=.true.; dumpSOL=.true.
      dumpPrivLeft=.true.;     dumpCoreLeftCell=.true.
      dumpPrivRite=.true.;     dumpCoreRiteCell=.true.
      if (mhdgeo <= 0 .and. isgrdsym .ne. 1) then
*       omit dumping the poloidal left leg and left core cell data - because there is none
        dumpPrivLeft=.false.;  dumpCoreLeftCell=.false.
      elseif (nxomit>=nxleg(1,1).and.nxomit<nxleg(1,1)+nxcore(1,1)) then
*       omit dumping only the poloidal left leg region
        dumpPrivLeft=.false.
      elseif (nxomit >= nxleg(1,1)+nxcore(1,1)+nxcore(1,2)) then
*       omit dumping left leg, left core, and right core; only dump right leg
        dumpPrivLeft=.false.;  dumpCoreLeftCell=.false.
        dumpPrivRite=.false.;  dumpCoreRiteCell=.false.
      endif
*  ---------------------------------------------------------------------
*     Core grid
*     painful because we have to find the nodes of a cell-based grid
*  ---------------------------------------------------------------------
*    Write out the core grid
      ixs=ixpt1(1)+1;            ixe=ixpt2(1)
      iys=0;                     iye=iysptrx1(1)
      lenx=ixe-ixs+1;            leny=iye-iys+1
      ALLOCATE(vs3d(0:lenx, 0:leny, 2)) ! 1 = R, 2 = Z
*    Left lower block
      vs3d(0:lenx-1,0:leny-1,1)=rm(ixs:ixe,0:iye, 1)
      vs3d(0:lenx-1,0:leny-1,2)=zm(ixs:ixe,0:iye, 1) - zmagx
*    Right line (same as left from periodicity)
      vs3d(lenx,0:leny-1,1)=rm(ixe,0:iye, 2)
      vs3d(lenx,0:leny-1,2)=zm(ixe,0:iye, 2) - zmagx
*    Inner line
      vs3d(0:lenx-1,leny,1)=rm(ixs:ixe,iye, 3)
      vs3d(0:lenx-1,leny,2)=zm(ixs:ixe,iye, 3) - zmagx
*    Right Inner node
      vs3d(lenx,leny,1)    =rm(ixe,iye, 4)
      vs3d(lenx,leny,2)    =zm(ixe,iye, 4) - zmagx

*    Dump
      h5in%units="m"; h5in%mesh="mesh-structured"; h5in%vsMD="mesh"
      CALL dump_h5(bbbid,"coreMesh",vs3d,h5in,h5err)
      DEALLOCATE(vs3d)
*  ---------------------------------------------------------------------
*     Core variables
*  ---------------------------------------------------------------------
      h5in%mesh="coreMesh"       ! Everybody as the same mesh here
      !-----------------------------------------------------------
      ! psi (dump all cell data
      !-----------------------------------------------------------
      ALLOCATE(vs2d(0:lenx, 0:leny))
*    Left lower block
      vs2d(0:lenx-1,0:leny-1)=psi(ixs:ixe,0:iye, 1)
*    Right line (same as left from periodicity)
      vs2d(lenx,0:leny-1)=psi(ixe,0:iye, 2)
*    Inner line
      vs2d(0:lenx-1,leny)=psi(ixs:ixe,iye, 3)
*    Right Inner node
      vs2d(lenx,leny)    =psi(ixe,iye, 4)
      h5in%units="Wb";    h5in%vsMD="psi";  h5in%vsCentering="nodal"
      CALL dump_h5(bbbid,"psiCore",vs2d,h5in,h5err)
      DEALLOCATE(vs2d)
      lenx=lenx-1;             leny = leny-1
      !-----------------------------------------------------------
      ! nis and ups are both for 1:nisp
      !-----------------------------------------------------------
      h5in%vsCentering="zonal"
      ALLOCATE(vs3d(0:lenx, 0:leny, nisp))
      vs3d(:,:,:)=nis(ixs:ixe,0:iye,:)
      h5in%units="m^-3";    h5in%vsMD="nis"
      CALL dump_h5(bbbid,"nisCore",vs3d,h5in,h5err)
      vs3d(:,:,:)=ups(ixs:ixe,0:iye,:)
      h5in%units="m/s";    h5in%vsMD="ups"
      CALL dump_h5(bbbid,"upsCore",vs3d,h5in,h5err)
      IF(read_diffs==1) THEN
      h5in%units="m^2/s"
      vs3d(:,:,:)=dif_use(ixs:ixe,0:iye,:); h5in%vsMD="ndiff"
      CALL dump_h5(bbbid,"dif_useCore", vs3d,h5in,h5err)
      ENDIF
      DEALLOCATE(vs3d)
      !-----------------------------------------------------------
      ! ngs  - dimension 1:ngsp
      !-----------------------------------------------------------
      ALLOCATE(vs3d(0:lenx, 0:leny, ngsp))
      vs3d(:,:,:)=ngs(ixs:ixe,0:iye,:)
      h5in%units="m^-3";    h5in%vsMD="ngs"
      CALL dump_h5(bbbid,"ngsCore",vs3d,h5in,h5err)
      DEALLOCATE(vs3d)
      !-----------------------------------------------------------
      ! tes, tis, phis
      !-----------------------------------------------------------
      ALLOCATE(vs2d(0:lenx, 0:leny))
      vs2d(:,:)=tes(ixs:ixe,0:iye)/1.602176487e-19
      h5in%units="eV";    h5in%vsMD="tes"
      CALL dump_h5(bbbid,"tesCore",vs2d,h5in,h5err)
      vs2d(:,:)=tis(ixs:ixe,0:iye)/1.602176487e-19
      h5in%units="eV";    h5in%vsMD="tis"
      CALL dump_h5(bbbid,"tisCore",vs2d,h5in,h5err)
      vs2d(:,:)=phis(ixs:ixe,0:iye)
      h5in%units="V";     h5in%vsMD="phis"
      CALL dump_h5(bbbid,"phisCore",vs2d,h5in,h5err)
      IF(read_diffs==1) THEN
      h5in%units="m^2/s"
      vs2d(:,:)=kye_use(ixs:ixe,0:iye); h5in%vsMD="kye"
      CALL dump_h5(bbbid,"kye_useCore", vs2d,h5in,h5err)
      vs2d(:,:)=kyi_use(ixs:ixe,0:iye); h5in%vsMD="kyi"
      CALL dump_h5(bbbid,"kyi_useCore", vs2d,h5in,h5err)
      ENDIF
      DEALLOCATE(vs2d)

*  ---------------------------------------------------------------------
*     Write out the sol grid
*  ---------------------------------------------------------------------
      lenx=nx+2;       leny=ny+2-iysptrx1(1)-1
      ixs=0;            ixe=nx+1  !=lenx-2
      iys=iysptrx1(1)+1;    iye=ny+1  !=leny+iysptrx1(1)
      ALLOCATE(vs3d(0:lenx, 0:leny, 2)) ! 1 = R, 2 = Z
*    Left lower block
      vs3d(0:lenx-1,0:leny-1,1)=rm(:,iys:iye, 1)
      vs3d(0:lenx-1,0:leny-1,2)=zm(:,iys:iye, 1) - zmagx
*    Right line (same as left from periodicity)
      vs3d(lenx,0:leny-1,1)=rm(ixe,iys:iye, 2)
      vs3d(lenx,0:leny-1,2)=zm(ixe,iys:iye, 2) - zmagx
*    Inner line
      vs3d(0:lenx-1,leny,1)=rm(ixs:ixe,iye, 3)
      vs3d(0:lenx-1,leny,2)=zm(ixs:ixe,iye, 3) - zmagx
*    Right Inner node
      vs3d(lenx,leny,1)=rm(ixe,iye, 4)
      vs3d(lenx,leny,2)=zm(ixe,iye, 4) - zmagx

*    Dump
      h5in%units="m";     h5in%mesh="mesh-structured";  h5in%vsMD="mesh"
      CALL dump_h5(bbbid,"solMesh",vs3d,h5in,h5err)
      DEALLOCATE(vs3d)
*  ---------------------------------------------------------------------
*     Sol variables
*  ---------------------------------------------------------------------
      h5in%mesh="solMesh"

      ! psi (dump entire cell data) !------------------------------
      ALLOCATE(vs2d(0:lenx, 0:leny))
*    Left lower block
      vs2d(0:lenx-1,0:leny-1)=psi(:,iys:iye, 1)
*    Right line (same as left from periodicity)
      vs2d(lenx,0:leny-1)=psi(ixe,iys:iye, 2)
*    Inner line
      vs2d(0:lenx-1,leny)=psi(ixs:ixe,iye, 3)
*    Right Inner node
      vs2d(lenx,leny)=psi(ixe,iye, 4)
      h5in%units="Wb";    h5in%vsMD="psi";  h5in%vsCentering="nodal"
      CALL dump_h5(bbbid,"psiSol",vs2d,h5in,h5err)
      DEALLOCATE(vs2d)

      lenx=lenx-1;               leny=leny-1

      h5in%vsCentering="zonal"
      ! nis and ups are both for 1:nisp !--------------------------
      ALLOCATE(vs3d(0:lenx, 0:leny, nisp))
      vs3d(:,:,:)=nis(ixs:ixe,iys:iye,:)
      h5in%units="m^-3";    h5in%vsMD="nis"
      CALL dump_h5(bbbid,"nisSol",vs3d,h5in,h5err)
      vs3d(:,:,:)=ups(ixs:ixe,iys:iye,:)
      h5in%units="m/s";    h5in%vsMD="ups"
      CALL dump_h5(bbbid,"upsSol",vs3d,h5in,h5err)
      IF(read_diffs==1) THEN
      h5in%units="m^2/s"; h5in%vsMD="ndiff"
      vs3d(:,:,:)=dif_use(ixs:ixe,iys:iye,:)
      CALL dump_h5(bbbid,"dif_useSol", vs3d,h5in,h5err)
      ENDIF
      DEALLOCATE(vs3d)

      ! ngs !------------------------------------------------------
      ALLOCATE(vs3d(0:lenx, 0:leny, ngsp))
      vs3d(:,:,:)=ngs(ixs:ixe,iys:iye,:)
      h5in%units="m^-3";    h5in%vsMD="ngs"
      CALL dump_h5(bbbid,"ngsSol",vs3d,h5in,h5err)
      DEALLOCATE(vs3d)
      ! WRITE(*,*) 'v---------------------------'
      ! WRITE(*,*) "tesSol"
      ! WRITE(*,*) ixs,ixe
      ! WRITE(*,*) iys,iye
      ! WRITE(*,*) 'v---------------------------'

      ! tes, tis, phis !-------------------------------------------
      ALLOCATE(vs2d(0:lenx, 0:leny))
      vs2d(:,:)=tes(ixs:ixe,iys:iye)/1.602176487e-19
      h5in%units="eV";    h5in%vsMD="tes"
      CALL dump_h5(bbbid,"tesSol",vs2d,h5in,h5err)
      vs2d(:,:)=tis(ixs:ixe,iys:iye)/1.602176487e-19
      h5in%units="eV";    h5in%vsMD="tis"
      CALL dump_h5(bbbid,"tisSol",vs2d,h5in,h5err)
      vs2d(:,:)=phis(ixs:ixe,iys:iye)
      h5in%units="V";     h5in%vsMD="phis"
      CALL dump_h5(bbbid,"phisSol",vs2d,h5in,h5err)
      IF(read_diffs==1) THEN
      h5in%units="m^2/s"
      vs2d(:,:)=kye_use(ixs:ixe,iys:iye); h5in%vsMD="kye"
      CALL dump_h5(bbbid,"kye_useSol", vs2d,h5in,h5err)
      vs2d(:,:)=kyi_use(ixs:ixe,iys:iye); h5in%vsMD="kyi"
      CALL dump_h5(bbbid,"kyi_useSol", vs2d,h5in,h5err)
      ENDIF
      DEALLOCATE(vs2d)

*  ---------------------------------------------------------------------
*    Set up dimensional quantities for getting the array indices right
*    Getting all of the geometries (slab, cylinder, etc.) is very tricky
*  ---------------------------------------------------------------------
      IF (dumpPrivLeft .AND. dumpPrivRite) THEN
        lenx=nx+2 +ixpt1(1)-ixpt2(1);       leny=iysptrx1(1)+1
        lxs=0;             lxe=ixpt1(1);    lys=0;     lye=leny-1
        rxs=ixpt1(1)+1;    rxe=lenx-1;      rys=0;     rye=leny-1
      ELSEIF (dumpPrivLeft .AND. .NOT. dumpPrivRite) THEN
        !SEK - Just return at this point since I don't have it working
        CALL close_group("bbb",bbbid,h5err)
        CALL close_h5file(fid,root_gid,h5err)
        RETURN
        !SEK
        lenx=ixpt1(1)+1;                    leny=iysptrx1(1)+1
        lxs=0;             lxe=ixpt1(1);    lys=0;     lye=leny-1
      ELSEIF (.NOT. dumpPrivLeft .AND. dumpPrivRite) THEN
        !SEK - Just return at this point since I don't have it working
        CALL close_group("bbb",bbbid,h5err)
        CALL close_h5file(fid,root_gid,h5err)
        RETURN
        !SEK
        !ixs=ixpt2(1)+1;            ixe=nx+1
        lenx=nx+1-ixpt2(1);                 leny=iysptrx1(1)+1
        rxs=0;             rxe=lenx-1;      rys=0;    rye=leny-1
      ENDIF
*  ---------------------------------------------------------------------
*    Write out the Priv grid
*  ---------------------------------------------------------------------
      ALLOCATE(vs3d(0:lenx, 0:leny, 2)) ! 1 = R, 2 = Z
      iys=0;                              iye =iysptrx1(1)
      IF (dumpPrivLeft) THEN
*    Left lower block
        ixs=0;            ixe=ixpt1(1)
        vs3d(lxs:lxe,lys:lye,1)=rm(ixs:ixe,iys:iye, 1)
        vs3d(lxs:lxe,lys:lye,2)=zm(ixs:ixe,iys:iye, 1) - zmagx
*    Inner line
        vs3d(lxs:lxe,lye+1,1)=rm(ixs:ixe,iye, 3)
        vs3d(lxs:lxe,lye+1,2)=zm(ixs:ixe,iye, 3) - zmagx
*    Right line (same as left from periodicity)
        vs3d(lxe+1,lys:lye,1)=rm(ixe,iys:iye, 2)
        vs3d(lxe+1,lys:lye,2)=zm(ixe,iys:iye, 2) - zmagx
*    Right Inner node
        vs3d(lxe+1,lye+1,1)=rm(ixe,iye, 4)
        vs3d(lxe+1,lye+1,2)=zm(ixe,iye, 4) - zmagx
      ENDIF
      IF (dumpPrivRite) THEN
        ixs=ixpt2(1)+1;            ixe=nx+1
*    Right lower block -- skip redundant line
        vs3d(rxs:rxe,rys:rye,1)=rm(ixs:ixe,iys:iye, 1)
        vs3d(rxs:rxe,rys:rye,2)=zm(ixs:ixe,iys:iye, 1)-zmagx
*    Right line
        vs3d(rxe+1,rys:rye,1)=rm(ixe,iys:iye, 2)
        vs3d(rxe+1,rys:rye,2)=zm(ixe,iys:iye, 2) - zmagx
*    Inner line
        vs3d(rxs:rxe,rye+1,1)=rm(ixs:ixe,iye, 3)
        vs3d(rxs:rxe,rye+1,2)=zm(ixs:ixe,iye, 3) - zmagx
*    Upper right corner of lower right block
        vs3d(rxe+1,rye+1,1)=rm(ixe,iye, 4)
        vs3d(rxe+1,rye+1,2)=zm(ixe,iye, 4) - zmagx
      ENDIF
*    Dump
      h5in%units="m";     h5in%mesh="mesh-structured";  h5in%vsMD="mesh"
      CALL dump_h5(bbbid,"privMesh",vs3d,h5in,h5err)
      DEALLOCATE(vs3d)
*  ---------------------------------------------------------------------
*     Priv variables
*  ---------------------------------------------------------------------
      h5in%mesh="privMesh"
      !-----------------------------------------------------------
      ! psi (dump entire cell data)
      !-----------------------------------------------------------
      ALLOCATE(vs2d(0:lenx, 0:leny))
      iys=0;                              iye =iysptrx1(1)
      IF (dumpPrivLeft) THEN
*    Left lower block
        ixs=0;            ixe=ixpt1(1)
        vs2d(lxs:lxe,lys:lye)=psi(ixs:ixe,iys:iye, 1)
*    Inner line
        vs2d(lxs:lxe,lye+1)=psi(ixs:ixe,iye, 3)
*    Right line (same as left from periodicity)
        vs2d(lxe+1,lys:lye)=psi(ixe,iys:iye, 2)
*    Right Inner node
        vs2d(lxe+1,lye+1)=psi(ixe,iye, 4)
      ENDIF
      IF (dumpPrivRite) THEN
        ixs=ixpt2(1)+1;            ixe=nx+1
*    Right lower block -- skip redundant line
        vs2d(rxs:rxe,rys:rye)=psi(ixs:ixe,iys:iye, 1)
*    Right line
        vs2d(rxe+1,rys:rye)=psi(ixe,iys:iye, 2)
*    Inner line
        vs2d(rxs:rxe,rye+1)=psi(ixs:ixe,iye, 3)
*    Upper right corner of lower right block
        vs2d(rxe+1,rye+1)=psi(ixe,iye, 4)
      ENDIF
      h5in%units="Wb";    h5in%vsMD="psi"; h5in%vsCentering="nodal"
      CALL dump_h5(bbbid,"psiPriv",vs2d,h5in,h5err)
      DEALLOCATE(vs2d)

      lenx=lenx-1;               leny = leny-1
      h5in%vsCentering="zonal"
*    Left lower block
      !-----------------------------------------------------------
      ! nis and ups are both for 1:nisp
      !-----------------------------------------------------------
      ALLOCATE(vs3d(0:lenx, 0:leny, nisp))
      ixs=0;            ixe=ixpt1(1)
      IF (dumpPrivLeft) vs3d(lxs:lxe,:,:)=nis(ixs:ixe,iys:iye,:)
      ixs=ixpt2(1)+1;            ixe=nx+1
      IF (dumpPrivRite) vs3d(rxs:rxe,:,:)=nis(ixs:ixe,iys:iye,:)
      h5in%units="m^-3";    h5in%vsMD="nis"
      CALL dump_h5(bbbid,"nisPriv",vs3d,h5in,h5err)
      ixs=0; ixe=ixpt1(1)
      IF (dumpPrivLeft) vs3d(lxs:lxe,:,:)=ups(ixs:ixe,iys:iye,:)
      ixs=ixpt2(1)+1;            ixe=nx+1
      IF (dumpPrivRite) vs3d(rxs:rxe,:,:)=ups(ixs:ixe,iys:iye,:)
      h5in%units="m/s";    h5in%vsMD="ups"
      CALL dump_h5(bbbid,"upsPriv",vs3d,h5in,h5err)
      IF(read_diffs==1) THEN
      h5in%units="m^2/s"; h5in%vsMD="ndiff"
      ixs=0; ixe=ixpt1(1)
      IF (dumpPrivLeft) vs3d(lxs:lxe,:,:)=dif_use(ixs:ixe,iys:iye,:)
      ixs=ixpt2(1)+1;            ixe=nx+1
      IF (dumpPrivRite) vs3d(rxs:rxe,:,:)=dif_use(ixs:ixe,iys:iye,:)
      CALL dump_h5(bbbid,"dif_usePriv",vs3d,h5in,h5err)
      ENDIF
      DEALLOCATE(vs3d)
      !-----------------------------------------------------------
      ! ngs
      !-----------------------------------------------------------
      ALLOCATE(vs3d(0:lenx, 0:leny, ngsp))
      ixs=0; ixe=ixpt1(1)
      IF (dumpPrivLeft) vs3d(lxs:lxe,:,:)=ngs(ixs:ixe,iys:iye,:)
      ixs=ixpt2(1)+1;            ixe=nx+1
      IF (dumpPrivRite) vs3d(rxs:rxe,:,:)=ngs(ixs:ixe,iys:iye,:)
      h5in%units="m^-3";    h5in%vsMD="ngs"
      CALL dump_h5(bbbid,"ngsPriv",vs3d,h5in,h5err)
      DEALLOCATE(vs3d)
      !-----------------------------------------------------------
      ! tes, tis, phis, psi
      ! WRITE(*,*) 'v---------------------------'
      ! WRITE(*,*) "tesPriv"
      ! WRITE(*,*) ixs,ixe
      ! WRITE(*,*) iys,iye
      ! WRITE(*,*) 'v---------------------------'

      !-----------------------------------------------------------
      ALLOCATE(vs2d(0:lenx, 0:leny))
      ixs=0; ixe=ixpt1(1)
      IF (dumpPrivLeft) vs2d(lxs:lxe,:)=tes(ixs:ixe,iys:iye)
      ixs=ixpt2(1)+1;            ixe=nx+1
      IF (dumpPrivRite) vs2d(rxs:rxe,:)=tes(ixs:ixe,iys:iye)
      vs2d=vs2d/1.602176487e-19
      h5in%units="eV";    h5in%vsMD="tes"
      CALL dump_h5(bbbid,"tesPriv",vs2d,h5in,h5err)
      ixs=0; ixe=ixpt1(1)
      IF (dumpPrivLeft) vs2d(lxs:lxe,:)=tis(ixs:ixe,iys:iye)
      ixs=ixpt2(1)+1;            ixe=nx+1
      IF (dumpPrivRite) vs2d(rxs:rxe,:)=tis(ixs:ixe,iys:iye)
      vs2d=vs2d/1.602176487e-19
      h5in%units="eV";    h5in%vsMD="tis"
      CALL dump_h5(bbbid,"tisPriv",vs2d,h5in,h5err)
      ixs=0; ixe=ixpt1(1)
      IF (dumpPrivLeft) vs2d(lxs:lxe,:)=phis(ixs:ixe,iys:iye)
      ixs=ixpt2(1)+1;            ixe=nx+1
      IF (dumpPrivRite) vs2d(rxs:rxe,:)=phis(ixs:ixe,iys:iye)
      h5in%units="V";     h5in%vsMD="phis"
      CALL dump_h5(bbbid,"phisPriv",vs2d,h5in,h5err)

      IF(read_diffs==1) THEN
      h5in%units="m^2/s"
      h5in%vsMD="kye"
      ixs=0; ixe=ixpt1(1)
      IF (dumpPrivLeft) vs2d(lxs:lxe,:)=kye_use(ixs:ixe,iys:iye)
      ixs=ixpt2(1)+1;            ixe=nx+1
      IF (dumpPrivRite) vs2d(rxs:rxe,:)=kye_use(ixs:ixe,iys:iye)
      CALL dump_h5(bbbid,"kye_usePriv",vs2d,h5in,h5err)
      h5in%vsMD="kyi"
      ixs=0; ixe=ixpt1(1)
      IF (dumpPrivLeft) vs2d(lxs:lxe,:)=kyi_use(ixs:ixe,iys:iye)
      ixs=ixpt2(1)+1;            ixe=nx+1
      IF (dumpPrivRite) vs2d(rxs:rxe,:)=kyi_use(ixs:ixe,iys:iye)
      CALL dump_h5(bbbid,"kyi_usePriv",vs2d,h5in,h5err)
      ENDIF
      DEALLOCATE(vs2d)
      !-----------------------------------------------------------
      ! Done!
      !-----------------------------------------------------------
      CALL close_group("bbb",bbbid,h5err)
      CALL close_h5file(fid,root_gid,h5err)
      WRITE(*,*) 'Dump of ', filename, ' concluded.'
      call flush(6)

      return
      end
c****** end of subroutine dump_vshdf5 *********************

c ------------------------------------------------------------------------
      subroutine restore_vshdf5(filename,doTransposeFlag,errval)
      Use(Dim)           # nisp,nusp,nfsp,ngsp, nhgsp
      Use(Interp)        # nis,ups,tes,tis,ngs,phis
      Use(Jacreorder)    # ireorder
      Use(Jacobian)      # nnzmx
      Use(Xpoint_indices)  # ixpt2
      Use(RZ_grid_info)  # rm, zm
      Use(hdf5_api)
      Use(hdf5)
      Use(Neqdskv)       # zmagx
c_mpi      Use(MpiVars)
c_mpi      Use(Npes_mpi)
      Use(Conduc)        # dif_use, tra_use, ...
      Use(Comtra)        # isbohmcalc facbni, facbup, facbni2, facbee
      Use(Ueint)         # diffread
      implicit none
      character*(*), intent(in) :: filename
      integer, intent(in) :: doTransposeFlag
      integer, intent(out) :: errval
      TYPE(hdf5InOpts) :: h5in
      TYPE(hdf5ErrorType) :: h5err
      character(60) :: description
      real, dimension(:,:,:), allocatable :: vs3d
      real, dimension(:,:), allocatable :: vs2d
      integer i,ifld,lid
      integer lenx, leny, ixs, ixe, iys, iye
      real(Size4) gettime, sec4
      integer(int4) :: ierr
      INTEGER(HID_T) :: root_gid,bbbid,fid
*  ---------------------------------------------------------------------
*     Setup and initialization
*  ---------------------------------------------------------------------
c_mpi      call mpi_barrier(uedgeComm, ierr)
      # WRITE(*,*) " Calling vshdf5_fcinit."
      # CALL vshdf5_fcinit(h5in,h5err)  ! IMPORTANT
      # WRITE(*,*) " vshdf5_fcinit returned."
      IF(mype/=0) return
      call vshdf5_inith5vars(h5in, h5err)
      ! WRITE(*,*) " vshdf5_inith5vars returned."
      h5in%wrd_type=H5T_NATIVE_DOUBLE
      h5in%doTranspose=.true.
      h5in%vsCentering="zonal"
      h5in%wrVsTime=.true.
      ! h5in%logUnit=99
      h5in%verbose=.false.
! c_mpi      WRITE(*,*) " mpi_comm_world = ", mpi_comm_world
! c_mpi      WRITE(*,*) " uedgeComm = ", uedgeComm
! c_mpi      WRITE(*,*) " mpi_comm_self = ", mpi_comm_self
c_mpi      h5in%comm = mpi_comm_self
! c__MPI      WRITE(*,*) " mpi_info_null = ", mpi_info_null
      h5in%info = 0
c__MPI      h5in%info=mpi_info_null
! c__MPI      WRITE(*,*) " h5in%info = ", h5in%info
*  ---------------------------------------------------------------------
*     Open the file and get to the right group
*  ---------------------------------------------------------------------
      ! WRITE(*, fmt='(a)I3') "restore_vshdf5: opening file '", filename,
      WRITE(*, *) "restore_vshdf5: opening file '", filename,
     &    "' with h5in%comm = ", h5in%comm
      CALL open_oldh5file(filename,fid,root_gid,h5in,h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
      ! everything is in the bbb subdirectory
      CALL open_group(root_gid,"bbb",bbbid,h5err)
      IF(h5err%errBool) WRITE(*,fmt='(/,a,/)') h5err%errorMsg
* jrc 12may09: How can one continue if there is an error opening the group?

*  ---------------------------------------------------------------------
*     Fill the fundamental arrays
*  ---------------------------------------------------------------------
* fmt='(a)' removes leading blank
      errval = 0
      WRITE(*, fmt='(a)') "Restoring the fundamental arrays."
      h5in%verbose=.true.
      CALL read_h5(bbbid,"nis",nis,h5in,h5err)
      IF(h5err%errBool) CALL dvshdf5_restore_stop("nis", errval)
      WRITE(*, fmt='(a)') "nis restored."
      CALL read_h5(bbbid,"ngs",ngs,h5in,h5err)
      IF(h5err%errBool) CALL dvshdf5_restore_stop("ngs", errval)
      WRITE(*, fmt='(a)') "ngs restored."
      CALL read_h5(bbbid,"tes",tes,h5in,h5err)
      IF(h5err%errBool) CALL dvshdf5_restore_stop("tes", errval)
      WRITE(*, fmt='(a)') "tes restored."
      CALL read_h5(bbbid,"tis",tis,h5in,h5err)
      IF(h5err%errBool) CALL dvshdf5_restore_stop("tis", errval)
      WRITE(*, fmt='(a)') "tis restored."
      CALL read_h5(bbbid,"ups",ups,h5in,h5err)
      IF(h5err%errBool) CALL dvshdf5_restore_stop("ups", errval)
      WRITE(*, fmt='(a)') "ups restored."
      CALL read_h5(bbbid,"phis",phis,h5in,h5err)
      IF(h5err%errBool) CALL dvshdf5_restore_stop("phis", errval)
      WRITE(*, fmt='(a)') "phis restored."

      IF ( errval .NE. 0 )  THEN
            CALL close_group("bbb",bbbid,h5err)
            CALL close_h5file(fid,root_gid,h5err)
            RETURN
      ENDIF

*  ---------------------------------------------------------------------
      CALL dvshdf5_diff_restore3d("dif_use",dif_use)
      CALL dvshdf5_diff_restore2d("kye_use",kye_use)
      CALL dvshdf5_diff_restore2d("kyi_use",kyi_use)
!      CALL dvshdf5_diff_restore3d("tra_use",tra_use)
!      CALL dvshdf5_diff_restore3d("vy_use",vy_use)
!      CALL dvshdf5_diff_restore2d("vyup_use",vyup_use)
!      CALL dvshdf5_diff_restore2d("vyte_use",vyte_use)
!      CALL dvshdf5_diff_restore2d("vyti_use",vyti_use)
!      CALL dvshdf5_diff_restore3d("fniyos_use",fniyos_use)
!      CALL dvshdf5_diff_restore2d("feeyosn_use",feeyosn_use)
!      CALL dvshdf5_diff_restore2d("feiyosn_use",feiyosn_use)
      !-----------------------------------------------------------
      ! Done!
      !-----------------------------------------------------------
      CALL close_group("bbb",bbbid,h5err)
      CALL close_h5file(fid,root_gid,h5err)

      CALL common_setup(1)

      return
      contains
c ------------------------------------------------------------------------
c      All of the diffusivities have the same structure.
c ------------------------------------------------------------------------
       SUBROUTINE dvshdf5_diff_restore2d(varName,array)
       CHARACTER(*), INTENT(IN) :: varName
       REAL, DIMENSION(:,:), INTENT(INOUT) :: array

       IF(read_diffs==1) THEN
       CALL read_h5(bbbid,TRIM(varName),array,h5in,h5err)
       IF(h5err%errBool) THEN

        WRITE(*, fmt='(a)') "Warning: Problem restoring "//TRIM(varNAME)
        WRITE(*, fmt='(a)') "         Possible inconsistent with inputs"
        WRITE(*, fmt='(a)') TRIM(h5err%errorMsg)
       ENDIF
       ENDIF
       RETURN
       END SUBROUTINE dvshdf5_diff_restore2d
       SUBROUTINE dvshdf5_diff_restore3d(varName,array)
       CHARACTER(*), INTENT(IN) :: varName
       REAL, DIMENSION(:,:,:), INTENT(INOUT) :: array

       IF(read_diffs==1) THEN
        CALL read_h5(bbbid,TRIM(varName),array,h5in,h5err)
        IF(h5err%errBool) THEN
         WRITE(*,fmt='(a)') "Warning: Problem restoring "//TRIM(varNAME)
         WRITE(*,fmt='(a)') "         Possible inconsistent with inputs"
         WRITE(*,fmt='(a)') TRIM(h5err%errorMsg)
        ENDIF
       ENDIF
       RETURN
       END SUBROUTINE dvshdf5_diff_restore3d


c ------------------------------------------------------------------------
c     If I fail on a restore then I want to stop
c ------------------------------------------------------------------------
       SUBROUTINE dvshdf5_restore_stop(varName, errval)
       CHARACTER(*), INTENT(IN) :: varName
       INTEGER, INTENT(OUT) :: errval

       WRITE(*,fmt='(a)') "ERROR: Problem restoring "//TRIM(varNAME)
       WRITE(*,fmt='(a)') TRIM(h5err%errorMsg)
       WRITE(*,fmt='(a)') "EXITING"
       errval = 1
       RETURN
       END SUBROUTINE dvshdf5_restore_stop
      end
c****** end of subroutine restore_vshdf5 *********************
c ------------------------------------------------------------------------
c     This is from uefacets.py.
c ------------------------------------------------------------------------
      subroutine common_setup(in_restart)
      Use(Ueint)       # zmagx
 c_mpi     Use(Npes_mpi)       # zmagx
      implicit none
      integer, intent(in) :: in_restart
*  ---------------------------------------------------------------------
*     Setup and initialization
*  ---------------------------------------------------------------------

*     Set up mass and charge arrays and flux surface average arrays needed
*     for the variable setting methods.  Must have arrays allocated and so must
*     be called after buildData.  Must know whether this is a restart,
*     to set bbb.restart in order to be able to call this from
*     initialize() or from restore()
      restart=in_restart

      if (ismpion==1 .and. restart == 1) then
        write(*,*)  "WARNING: For parallel run, averages "
        write(*,*)  "are NOT computed by commonsetup."
cc        return
      endif

* If bbb.ueinit is here, then uedge runs
      if (ismpion == 0) then
        call ueinit
      elseif (mype == 0) then
        call globalmesh
      endif

* SEK: There was a lot of stuff here that is confusing me.
* I don't think it was needed, but who knows

* In above, making use of thin-ness of guard cell so  needn't worry
*  abut using centered cell value of r.

      return
      end
c****** end of subroutine common_setup *********************
c ------------------------------------------------------------------------
