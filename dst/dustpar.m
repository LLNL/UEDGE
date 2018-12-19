      subroutine dstprr
c
c--------------------------------------
c subroutine prepares PSI parameter
c--------------------------------------
c
      implicit none

      Use (dustinp)
      Use (std_cns)
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      Use (psiparloc)

      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)

      integer  i,k
      integer  psierr
      integer  psitabr
      external psitabr
      real*8   psimassa,psi_Aew,psi_izpot
      real*8   psi_Tmelt,psi_Hmelt_ev,psi_Tsub
      integer  psi_melt
      external psimassa,psi_Aew,psi_izpot
      external psi_Tmelt,psi_Hmelt_ev,psi_Tsub
      external psi_melt
c
c---input interpolation mesh parameters 
c
      call psimatinp(psierr)
      if (psierr .ne. 0) then
       print *, 'Error reading meterial parameters =', psierr
c       stop
       return
      endif
      call psitabicro(psierr)
      if (psierr .ne. 0) then
       print *, 'Error initializing of material parameters =', psierr
c       stop
       return
      endif
      
c---prepare psi parameters
      call psi_tmesh_ini
      call psi_tmesh_par
      psierr=psitabr(Nmat_ptm)	 ! number of target materials
      if (psierr .ne. 0) then
       print *, 'error of reading the interpolation tables =', psierr
c       stop
       return
      endif

      k=0
      do 20, i=1,Nmat_pt
       if (IZmatter(i) .eq. zmat) k=i
20    continue
      if (k .gt. 0) then
       imater=k
      else
       print *, 'No PSI data for dust material!'
c       stop
       return
      endif
      
      call psi_tmesh_imat(imater)
      md=psimassa(zmat)
      cwfunc  =psi_Aew     (zmat)
      canmelt =psi_melt    (zmat)
      Tmm0    =psi_Tmelt   (zmat)+273.15d0
      tsub    =psi_Tsub    (zmat)+273.15d0
      xmeltene=psi_Hmelt_ev(zmat)*9.6472d4/md

      
      k=0
      do 10, i=1,Nprj_pt
       if (IZprjctl(i) .eq. zip) k=i
10    continue
      if (k .gt. 0) then
       ip=k
       mp=AZprjctl(ip)
       mpnrm=mp*mol1nrm
       invmpnrm=1.d0/mpnrm
      else
       print *, 'No PSI data for main plasma component!'
c       stop
       return
      endif

      ipot=psi_izpot(zip)

      end
