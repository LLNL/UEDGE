      subroutine dstrji

      implicit none

      Use (dusttrj)

      integer i

      kaa=0
      do 1 i=1,naa
         ra(i)=0.d0
         za(i)=0.d0
        tha(i)=0.d0
        vra(i)=0.d0
        vza(i)=0.d0
        vta(i)=0.d0
        rda(i)=0.d0
        tda(i)=0.d0
       msda(i)=0.d0
        wta(i)=0.d0
        tta(i)=0.d0
        spa(i)=0.d0
        sta(i)=0.d0
        fpa(i)=0.d0
        fma(i)=0.d0
        qpa(i)=0.d0
        qma(i)=0.d0
       fida(i)=0.d0
       xliquda(i)=0.d0
       lghta(i)=0.d0
1     continue

      kjm=0
      kjms=0

      end  ! dstrji
c=======================================================================


      subroutine dstrjs(ixb,iyb,r0,z0,th0,avr,avz,avt,wes10,nn)

      implicit none

c      Use (dustdim)
      Use (dustinp)
      Use (dustcom)
      Use (dustcgm)
      Use (dusttrj)
      Use (dustcntrl)
      Use (std_cns)
      Use (psitab)
      Use (psitab_s)
      Use (psiparloc)

      integer ixb,iyb,i,nn
      real*8  r0,z0,th0,avr,avz,avt,wes10
      real*8  ran,logrd
      real*8  ro

      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     ,            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)


      rb=r0
      zb=z0
      thb=th0

      if (jstt .eq. 0) then
       rdb=rd0
       tdb=td0
       xliqudb=xliqud0
       call psi_tdmesh_int(tdb)
       call psi_tdmesh_vall
cc       call psi_tdmesh_rocp
       ro=rotd
      else
       tdb=twall(ixb,iyb)
       xliqudb=0.d0
       call psi_tdmesh_int(tdb)
       call psi_tdmesh_vall
cc       call psi_tdmesh_rocp
       ro=rotd
c --- cast radius
       rdb=rdist(nn)
c --- cast speed
       call std_rndd3x(ran)
       i=1
55     if (ran .gt. fvdist(i)) then
        i=i+1
        goto 55
       endif
       vp0=fdvd*(i-2+(ran-fvdist(i-1))/(fvdist(i)-fvdist(i-1)))
      endif
      vrb=avr*vp0
      vzb=avz*vp0
      vtb=avt*vp0
      msdb=pi4*onethird*rdb*rdb*rdb*ro*1.d-12
      
      tb      = 0.d0
      spb     = 0.d0
      stb     = 0.d0
      fpb     = 0.d0
      fmb     = 0.d0
      qpb     = 0.d0
      qmb     = 0.d0
      fidb    = 0.d0
      lghtb   = 0.d0
      wtb     = wes10

      kjm=1
      call dstrjif1

      return
      end  ! dstrjs
c=======================================================================


      subroutine dstrjp(kod,ix,iy,re,ze,the,dte,wwe,wwe2)

      implicit none

c      Use (dustdim)
      Use (dustinp)
      Use (dusttrj)
      Use (dustcntrl)
      Use (dustcom)
      Use (dustcgm)
      Use (dustout)
      Use (dustcur)
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
     ,            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      integer kod,ix,iy
      real*8 re,ze,the,dte,wwe,wwe2

      integer l,m,j
      integer idstpin
      integer k,kd,i,red,dtred
      integer nnn,nitert
      real*8 std_erf
      real*8 popt,tk1ev
      real*8 ro,cp,p,dp,pr1,pr2,pz1,pz2
      real*8 s,g
      real*8 vtn,un,un2,un2exp,unerf,etan,etaa,etaf
      real*8 unia(nimpm),un2ia(nimpm),un2expia(nimpm),
     ,       unerfia(nimpm),etiaf(nimpm)
      real*8 dr,dz,dth,dvz,dvr,dvt,dt
      real*8 fz,fr,ft
      real*8 zd,arad,invmsd,
     ,       qdh,qdc,dtds,dmsd,msdin,dtdx
      real*8 xgam,dgamm,egam,gamt
      real*8 eis,eas,gep,gee,gez
      real*8 eiz(nimpzm),eia(nimpm)
      real*8 gsub,esub,gcsp,gres,gpsp,epsp
      real*8 rnp,rep,qgamp,qgame,qgama,qgamt
      real*8 rna,rea
      real*8 qgamz,qgamia
      real*8 rnia(nimpm),reia(nimpm)
      real*8 rnz(nimpzm),rez(nimpzm)
      real*8 efinrm
      real*8 tetoti,invr
      real*8 tetotz(nimpzm)
      real*8 radc
      real*8 tmp,ekin,vabs,vpol,vtrs
      real*8 derr,tdtmp,msdtmp,rdtmp
      real*8 aui2,invaui,auz2,xisqr
      real*8 enume1,enume2
      real*8 c1,c2,c3,c4
      real*8 Tmmt,Tmpt,deltat,dttmp,dtdl,xltmp,lttmp
      real*8 ycspfc,yresfc
      real*8 psiyresfc,psiycsptf,psi_Hmelt_ev
      external psiyresfc,psiycsptf,psi_Hmelt_ev
      real*8 rs,ms,tds,qs,es,vrs,vzs,vts,vs,ros,cps,vps,drdts,dtdts,
     *       xliquds,lghts
      common /afterstep/ rs,ms,tds,qs,es,vrs,vzs,vts,vs,ros,cps,vps,
     *                   drdts,dtdts,xliquds,lghts

      real*8 eptd1,eptd2,s1
      integer itd, iradius,iradtd

      data popt   /3.d-1/        ! minimal poloidal step to search transitions in mesh
      data tk1ev  /11604.5d0/    ! temperature equivalent to 1 eV
      data radc   /679.d0/       ! threshold  for small body radiation in rd(mkm)*td(K)
      data derr   /5.d-2/        ! threshold for relative change of md & td to reduce dt
      data dtred  /100/          ! reduction factor for time step 
      data deltat /1.d0/         ! +- variation of dust temperature during melting, K

      Tmmt=Tmm0-deltat
      Tmpt=Tmm0+deltat

      tetoti=te/ti

      if (jimp .gt. 0) then
       do 40, m=1,nimpz
        tetotz(m)=te/tiz(m)
40     continue
      endif

      kd  = 0
      k   = 0
      r   = rb
      z   = zb
      th  = thb
      vr  = vrb
      vz  = vzb
      vt  = vtb
      rd  = rdb
      td  = tdb
      msd = msdb
      wt  = wtb
      t   = tb
      sp  = spb
      st  = stb
      fpd = fpb
      fmd = fmb
      qpd = qpb
      qmd = qmb
      fid = fidb
      xliqud = xliqudb
      lght = lghtb
      
      msdin=msd

      dt  = dte
      p   = 0.d0
      kaa = 0
      kd  = 0
      kod = 0

      do 1 i=1,naa

             kaa=kaa+1
         ra(kaa)=r
         za(kaa)=z
        tha(kaa)=th
        vra(kaa)=vr
        vza(kaa)=vz
        vta(kaa)=vt
        rda(kaa)=rd
        tda(kaa)=td
       msda(kaa)=msd
        wta(kaa)=wt
        tta(kaa)=t
        spa(kaa)=sp
        sta(kaa)=st
        fpa(kaa)=fpd
        fma(kaa)=fmd
        qpa(kaa)=qpd
        qma(kaa)=qmd
       fida(kaa)=fid
       xliquda(kaa)=xliqud
       lghta(kaa)=lght

       vpol=vr*vr+vz*vz
       vabs=vpol+vt*vt
       ekin=.5d0*msd*vabs
       vpol=dsqrt(vpol)
       vabs=dsqrt(vabs)

c    --initializing of interpolation on Td mesh
       call psi_tdmesh_int(td)
       call psi_tdmesh_vall
ccc       call psi_tdmesh_rocp
c    -- heat capacity, J/ g grad K
       ro=rotd
       cp=ksicp*cptd

       call psi_tmesh_int(te)
c------ion friction force and equilibruin charge------------------------
       call dstdrag

c    - equilibrium dust charge Qd/e, Zd=[rd/rBohr]*[T(eV)/2Ry]
       zd=-6.9446d2*rd*lambda*te
c       fidd=-lambda
       fidd=zd

c    -- ion projectiles
        g=2.d0*td/tk1ev
        efinrm=lambda*tetoti
        aui2=aui*aui
        invaui=1.d0/aui
        if (lambda .ge. 0.d0) then
         if (aui .lt. 5.d-3) then
          gep=.5d-20*np*dsqrt(2.d0*tierg*invmpnrm)*invpisqrt
     *       *ti*(2.d0+efinrm+(2.d0+efinrm*onethird)*aui2)
         else
          enume1=invpisqrt*(5.d0+2.d0*(aui2+efinrm))*dexp(-aui2)
          enume2=1.5d0+2.d0*aui2*(3.d0+aui2)+efinrm*(1.d0+aui2+aui2)
          enume2=enume2*std_erf(aui)*invaui
          gep=.125d-20*np*dsqrt(2.d0*tierg*invmpnrm)*ti
          gep=gep*(enume1+enume2)
         endif
        else
         if (aui .lt. 5.d-3) then
          gep=.5d-20*np*dsqrt(2.d0*tierg*invmpnrm)*invpisqrt
     *       *ti*dexp(efinrm)*(2.d0-efinrm
     *       +(3.d0+efinrm*(efinrm-2.5d0))*twothirds*aui2)
         else
          xisqr=dsqrt(-efinrm)
          c3=2.d0*aui2
          c1=5.d0+c3
          c2=(3.d0+c3)*xisqr*invaui
          c3=c1-c2
          c4=c1+c2
          c1=aui+xisqr
          c2=aui-xisqr
          enume1=c3*dexp(-c1*c1)+c4*dexp(-c2*c2)
          enume1=enume1*invpisqrt
          enume2=1.5d0+2.d0*aui2*(3.d0+aui2)+efinrm*(1.d0+2.d0*aui2)
          enume2=enume2*(std_erf(c1)+std_erf(c2))*invaui
          gep=6.25d-22*np*dsqrt(2.d0*tierg*invmpnrm)*ti
          gep=gep*(enume1+enume2)
         endif
        endif
        gep=ksip*gep
        eis=gep/gamp+lambda*te ! averaged energy of ion striking dust surface

c    -- atom projectiles
        vtn=2.d0*taerg*invmpnrm
        un=(uza-vz)*(uza-vz)
        un=un+(ura-vr)*(ura-vr)
        un=un+(utha-vt)*(utha-vt)
        un2=un/vtn
        un=dsqrt(un2)
        un2exp=dexp(-un2)
        unerf=std_erf(un)
        if (un .lt. 5.d-3) then
         etaa=2.d0+un2+un2
         etaf=1.d0+onethird*un2
        else
         etaa=((un2+3.d0)*un2+.75d0)/un*unerf*pisqrt
         etaa=.5d0*(etaa+(un2+2.5d0)*un2exp)
         etaf=.5d0*(pisqrt*(un+.5d0/un)*unerf+un2exp)
        endif
      
        eas=etaa*ta/etaf ! everaged energy of atom striking dust surface
c    -- flux of neutral atoms
        gama=0.5d-20*na*dsqrt(vtn*invpi)*etaf !atom flux 10^20/cm2/s when [na]=cm^-3

c    -- plasma particle flux on dust surface, 10^20/cm2/s
        gamt=gamp+gama

c    -- impurity particle flux on dust surface, 10^20/cm2/s
        if (jimp .gt. 0) then
         do 20, l=1,nimp
          vtn=2.d0*tiaerg(l)/mianrm(l)
          unia(l)=(uzia(l)-vz)*(uzia(l)-vz)
          unia(l)=unia(l)+(uria(l)-vr)*(uria(l)-vr)
          unia(l)=unia(l)+(uthia(l)-vt)*(uthia(l)-vt)
          un2ia(l)=unia(l)/vtn
          unia(l)=dsqrt(un2ia(l))
          un2expia(l)=dexp(-un2ia(l))
          unerfia(l)=std_erf(unia(l))
          if (unia(l) .lt. 5.d-3) then
           etaa=2.d0+un2ia(l)+un2ia(l)
           etiaf(l)=1.d0+onethird*un2ia(l)
          else
           etaa=((un2ia(l)+3.d0)*un2ia(l)+.75d0)
     *         /unia(l)*unerfia(l)*pisqrt
           etaa=.5d0*(etaa+(un2ia(l)+2.5d0)*un2expia(l))
           etiaf(l)=.5d0*(pisqrt*(unia(l)
     *             +.5d0/unia(l))*unerfia(l)+un2expia(l))
          endif
          eia(l)=etaa*tia(l)/etiaf(l) ! everaged energy of atom striking dust surface
	  
          gamia(l)=.5d-20*nia(l)*dsqrt(vtn*invpi)*etiaf(l) !atom flux 10^20/cm2/s when [na]=cm^-3
	  
          gamt=gamt+gamia(l)
20       continue

         do 30, m=1,nimpz
          efinrm=zimp(m)*lambda*tetotz(m)
          auz2=auz(m)*auz(m)
          if (lambda .ge. 0.d0) then
           if (auz(m) .lt. 5.d-3) then
            gez=.5d-20*niz(m)*dsqrt(2.d0*tizerg(m)/mznrm(m))
     *         *invpisqrt*tiz(m)*(2.d0+efinrm
     *         +(2.d0+efinrm*onethird)*auz2)
           else
            enume1=invpisqrt*(5.d0+2.d0*(auz2+efinrm))*dexp(-auz2)
            enume2=3.d0+4.d0*auz2*(3.d0+auz2)
     *            +2.d0*efinrm*(1.d0+2.d0*auz2)
            enume2=enume2*std_erf(auz(m))/(2.d0*auz(m))
            gez=.125d-20*niz(m)*dsqrt(2.d0*tizerg(m)/mznrm(m))*tiz(m)
            gez=gez*(enume1+enume2)
           endif
          else
           if (auz(m) .lt. 5.d-3) then
            gez=.5d-20*niz(m)*dsqrt(2.d0*tizerg(m)/mznrm(m))
     *         *invpisqrt*tiz(m)*dexp(efinrm)*(2.d0-efinrm
     *         +(6.d0+efinrm*(2.d0*efinrm-5.d0))*onethird*auz2)
           else
            xisqr=dsqrt(-efinrm)
            c3=2.d0*auz2
            c1=5.d0+c3
            c2=(3.d0+c3)*xisqr/auz(m)
            c3=c1-c2
            c4=c1+c2
            c1=auz(m)+xisqr
            c2=auz(m)-xisqr
            enume1=c3*dexp(-c1*c1)+c4*dexp(-c2*c2)
            enume1=enume1*invpisqrt
            enume2=3.d0+4.d0*auz2*(3.d0+auz2)
     *            +2.d0*efinrm*(1.d0+2.d0*auz2)
            enume2=enume2*(std_erf(c1)+std_erf(c2))/(2.d0*auz(m))
            gez=6.25d-22*niz(m)*dsqrt(2.d0*tizerg(m)/mznrm(m))*tiz(m)
            gez=gez*(enume1+enume2)
           endif
          endif
  
          gez=ksip*gez
          eiz(m)=gez/gamz(m)+lambda*te*zimp(m)
          gamt=gamt+gamz(m)
30       continue
        endif

c    -- redusction coefficient for RES and chemical sruttering yeilds
        yresfc=psiyresfc(gamt)
        ycspfc=psiycsptf(gamt,td)

c    -- particle and energy fluxes out off dust (10^20 atoms/cm2/s and eV)
        xgam  = 0.d0
        dgamm = 0.d0
        egam  = 0.d0

c--- plasma particle and energy reflection coefficients and
c    plasma particle and energy fluxes out off dust
c    (10^20 atoms/cm2/s and eV)

c    -- ion projectiles
        call psi_tmesh_int(.5d0*eis)
        call psi_tmesh_yt
        call psi_tmesh_ytt

        rnp=frn(ip)*ksirn
        rep=fre(ip)

c    -- associated with chemical sputtering
        gcsp=fycsp(ip)
        gcsp=gcsp*gamp*ycspfc
        xgam=xgam+gcsp
        egam=egam+gcsp*g

c    -- associated with RES
        gres=fyres(ip)
        gres=gres*gamp*yresfc
        xgam=xgam+gres
        egam=egam+gres*g

c    -- associated with physical sputtering
        gpsp=fypsp(ip)
        gpsp=gpsp*gamp
        epsp=fepsp(ip)
        xgam=xgam+gpsp
        egam=egam+gpsp*epsp
          
c    -- atom projectiles
        call psi_tmesh_int(.5d0*eas)
        call psi_tmesh_yt
        call psi_tmesh_ytt

        rna=frn(ip)*ksirn
        rea=fre(ip)

c    -- associated with chemical sputtering
        gcsp=fycsp(ip)
        gcsp=gcsp*gama*ycspfc
        xgam=xgam+gcsp
        egam=egam+gcsp*g

c    -- associated with RES
        gres=fyres(ip)
        gres=gres*gama*yresfc
        xgam=xgam+gres
        egam=egam+gres*g

c    -- associated with physical sputtering due to atoms
        gpsp=fypsp(ip)
        gpsp=gpsp*gama
        epsp=fepsp(ip)
        xgam=xgam+gpsp
        egam=egam+gpsp*epsp
        
c--- impurity particle energy reflection coefficients and
c--- energy fluxes out off dust
c---  (10^20 atoms/cm2/s and eV)
        if (jimp .gt. 0) then
c    -- ion projectiles
         do 50, m=1,nimpz
          call psi_tmesh_int(.5d0*eiz(m))
          call psi_tmesh_yt
          call psi_tmesh_ytt

          rnz(m)=frn(iiz(m))*ksirn
          rez(m)=fre(iiz(m))

c    -- associated with chemical sputtering
          gcsp=fycsp(iiz(m))
          gcsp=gcsp*gamz(m)*ycspfc
          xgam=xgam+gcsp
          egam=egam+gcsp*g

c    -- associated with RES
          gres=fyres(iiz(m))
          gres=gres*gamz(m)*yresfc
          xgam=xgam+gres
          egam=egam+gres*g

c    -- associated with physical sputtering
          gpsp=fypsp(iiz(m))
          gpsp=gpsp*gamz(m)
          epsp=fepsp(iiz(m))
          xgam=xgam+gpsp
          egam=egam+gpsp*epsp
          
c    -- deposition on dust particle
          dgamm=dgamm+gamz(m)*rnz(m)*mz(m)
50       continue
c--- end of impurity ions
 
c    -- impurity atom projectiles
         do 60, l=1,nimp
          call psi_tmesh_int(.5d0*eia(l))
          call psi_tmesh_yt
          call psi_tmesh_ytt

          rnia(l)=frn(iia(l))*ksirn
          reia(l)=fre(iia(l))

c    --- associated with chemical sputtering
          gcsp=fycsp(iia(l))
          gcsp=gcsp*gamia(l)*ycspfc
          xgam=xgam+gcsp
          egam=egam+gcsp*g

c    --- associated with RES
          gres=fyres(iia(l))
          gres=gres*gamia(l)*yresfc
          xgam=xgam+gres
          egam=egam+gres*g

c    --- associated with physical sputtering due to atoms
          gpsp=fypsp(iia(l))
          gpsp=gpsp*gamia(l)
          epsp=fepsp(iia(l))
          xgam=xgam+gpsp
          egam=egam+gpsp*epsp

c    --- deposition on dust particle
          dgamm=dgamm+gamia(l)*rnia(l)*mia(l)
60       continue
c    -- end of impurity atom projectiles
        endif 

c    -- associated with thermionic emission of electrons
        egam=egam+gamthe*(g+cwfunc)

c    -  associated with secondary e-e emission
        egam=egam+3.d0*cwfunc*gamsee

c    -- associated with plasma kinetic & potential energies
        if (lambda .ge. 0d0) then
         gee=.5d-20*ne*dsqrt(8.d0*teerg/(pi0*menrm))
     *      *te*(1+.5d0*lambda)*dexp(-lambda)
        else
         gee=.5d-20*ne*dsqrt(8.d0*teerg/(pi0*menrm))
     *      *te*(1-.5d0*lambda)
        endif
        gee=ksip*gee
        qgame=gee-lambda*te*game
        qgamp=(eis*(1.d0-rep)+ipot)*gamp
        qgama=eas*gama*(1.d0-rea)

        if (jimp .gt. 0) then
         qgamia=0.d0
         do 70, l=1,nimp
          qgamia=qgamia+eia(l)*gamia(l)*(1.d0-reia(l))
70       continue
         qgamz=0.d0
         do 80, m=1,nimpz
          qgamz=qgamz+(eiz(m)*(1.d0-rez(m))+zpot(m))*gamz(m)
80       continue
        endif

c  --- total heat flux onto dust particle surface
        qgamt=qgamp+qgama+qgame
        if (jimp .gt. 0) qgamt=qgamt+qgamia+qgamz  ! from impurities

c  --- plasma heat fluxes onto dust particle surface, erg/cm2/s
        qdh=qgamt*1.602d8*ksiq

c  --- associated with thermal sublimation
        gsub=gstd
        esub=hstd
ccc        xgam=xgam+gsub
ccc        egam=egam+gsub*esub

c---heat fluxes outfrom dust particle surface, erg/cm2/s
c   associated with desintegration and evaporation
        qdc=(egam+gsub*esub)*1.602176462d8
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   black-body radiation (SB=5.6704 10^-8 W/m2/K4)

        
        eptd1=psi_ep(1,imat) ! constant value at lowest temperature (default)
        eptd2=eptd1
        g=td-tw0

        if (g .gt. 0.d0) then
         s=td*td
         s=s*s
         s=s-tw04
         if ((jtheps .eq. 0) .or. (jtheps .eq. 1)) then
          eptd2=eptd1
          if ((jtheps .eq. 1) .and. (rd*td .le. radc*eptd1)) then
           s1=td*td*td*td*td-tw05
           s1=8.356d-8*rd*s1
           eptd2=s1/(5.6704d-5*s)
          endif
         elseif (jtheps .eq. 2) then
          eptd2=eptd
         elseif (jtheps .eq. 3) then
          call psi_rdmesh_int(rd)
          call psi_tdemesh_int(td)
          call psi_mesh_theps
          if (theps .le. 0) then
           print *, '***Error: wrong Mie emissivity',theps
c           theps=eptd1
c           stop
           return
          endif
          eptd2=theps
         else
          print *, '***Error: wrong emissivity selector',jtheps
c          stop
          return
         endif
        else
         s=0.d0
         eptd2=0.d0
        endif
        s=5.6704d-5*eptd2*s
c----
        qdc=qdc+s

3232    format(1p20e15.5)

        lghtd=s*pi4*rd*rd*1.d-8


c   T Cp d(ro 4/3 pi r^3)/dt
        s=1.d3*mol1*(dgamm-(xgam+gsub)*md)
     *   *(td*cp+xliqud*xmeltene)
c        s=1.d3*mol1*(dgamm-(xgam+gsub)*md)*xliqud*xmeltene
        qdc=qdc+s
	
c--- dust particle temperature and molten fraction
        dtdx=0.d0
        if((Tmmt    .lt. td  ) .and.
     *     (td      .lt. Tmpt) .and.
     *     (canmelt .gt. 0   )) then
         s=1.d3*rd*ro*xmeltene
         dtdx=3.d0*dt*(qdh-qdc)/s
         xliqudd=xliqud+dtdx
         tdd=(xliqudd-0.5d0)*2.d0*deltat+Tmm0
         dtds=tdd-td
         xliqudd=min(max(0.d0,xliqudd),1.0d0)
        else
         s=1.d3*rd*ro*cp
         dtds=3.0d0*dt*(qdh-qdc)/s
         tdd=td+dtds
         if(tdd     .le. Tmmt) xliqudd=0.d0
         if(Tmpt    .le. tdd ) xliqudd=1.d0
         if(canmelt .eq. 0   ) xliqudd=0.d0
        endif
        
c--- dust mass calculation
        dmsd=dt*(dgamm-(xgam+gsub)*md)*mol1*pi4*rd*rd*1.d-12
        msdd=msd+dmsd
	
        dtdl=xliqudd-xliqud

c--- reduction of time step for fast sublimation -----------------------
        if ((dabs(dtds) .gt. dabs(derr*td))  .or.
     *      (dabs(dmsd) .gt. dabs(derr*msd)) .or.
     *      (dabs(dtdl) .gt. dabs(derr)))
     *  then
         tdtmp=td
         msdtmp=msd
         rdtmp=rd
         xltmp=xliqud
         lttmp=lght
         tmp=wwe2*dt/dtred
         do 10, red=1,dtred
          if (jstt .gt. 0) then
           nsts (ix,iy)=nsts (ix,iy)+tmp
           rsts (ix,iy)=rsts (ix,iy)+rdtmp *tmp
           msts (ix,iy)=msts (ix,iy)+msdtmp*tmp
           tdsts(ix,iy)=tdsts(ix,iy)+tdtmp *tmp
           qsts (ix,iy)=qsts (ix,iy)+zd    *tmp
           ests (ix,iy)=ests (ix,iy)+ekin  *tmp
           vrsts(ix,iy)=vrsts(ix,iy)+vr    *tmp
           vzsts(ix,iy)=vzsts(ix,iy)+vz    *tmp
           vtsts(ix,iy)=vtsts(ix,iy)+vt    *tmp
           vsts (ix,iy)=vsts (ix,iy)+vabs  *tmp
           xlsts(ix,iy)=xlsts(ix,iy)+xltmp *tmp
           ltsts(ix,iy)=ltsts(ix,iy)+lttmp *tmp
c --- gather stat for distribution functions
           j=idint((dlog(rdtmp)/aln10-minlogrd)/dlogrd)
           if (j .le. 0) then
            rdfnc (ix,iy,0)=rdfnc (ix,iy,0)+tmp
            rdfnct(      0)=rdfnct(      0)+tmp
           elseif (j .lt. frdim) then
            rdfnc (ix,iy,j)=rdfnc (ix,iy,j)+tmp
            rdfnct(      j)=rdfnct(      j)+tmp
           else
            rdfnc (ix,iy,frdim)=rdfnc (ix,iy,frdim)+tmp
            rdfnct(      frdim)=rdfnct(      frdim)+tmp
           endif
           j=idint(vabs/dvd)
           if ((j .lt. fvdim) .and. (j .ge. 0)) then
            vdfnc (ix,iy,j)=vdfnc (ix,iy,j)+tmp
            vdfnct(      j)=vdfnct(      j)+tmp
           else
            vdfnc (ix,iy,fvdim)=vdfnc (ix,iy,fvdim)+tmp
            vdfnct(      fvdim)=vdfnct(      fvdim)+tmp
           endif
           vtrs=vt+0.5d0*maxvd
           if (vtrs.ge.0.d0) then
            j=idint(vtrs/dvd)
            if ((j .lt. fvdim) .and. (j .ge. 0)) then
             vdfnctr (ix,iy,j)=vdfnctr (ix,iy,j)+tmp
             vdfnctrt(      j)=vdfnctrt(      j)+tmp
            else
             vdfnctr (ix,iy,fvdim)=vdfnctr (ix,iy,fvdim)+tmp
             vdfnctrt(      fvdim)=vdfnctrt(      fvdim)+tmp
            endif
           else
            vdfnctr (ix,iy,0)=vdfnctr (ix,iy,0)+tmp
            vdfnctrt(      0)=vdfnctrt(      0)+tmp
           endif
           j=idint(vpol/dvd)
           if ((j .lt. fvdim) .and. (j .ge. 0)) then
            vdfncpl (ix,iy,j)=vdfncpl (ix,iy,j)+tmp
            vdfncplt(      j)=vdfncplt(      j)+tmp
           else
            vdfncpl (ix,iy,fvdim)=vdfncpl (ix,iy,fvdim)+tmp
            vdfncplt(      fvdim)=vdfncplt(      fvdim)+tmp
           endif
          endif
	  
          call psi_tdmesh_int(tdtmp)
          call psi_tdmesh_vall
cc          call psi_tdmesh_rocp
          ro=rotd
          cp=ksicp*cptd
          gsub=gstd
          esub=hstd
          qdc=(egam+gsub*esub)*1.602176462d8
cc====

          eptd1=psi_ep(1,imat) ! constant value at lowest temperature
          eptd2=eptd1
          g=tdtmp-tw0

          if (g .gt. 0.d0) then
           s=tdtmp*tdtmp
           s=s*s
           s=s-tw04
           if ((jtheps .eq. 0) .or. (jtheps .eq. 1)) then
            eptd2=eptd1
            if ((jtheps .eq. 1) .and.
     *          (rdtmp*tdtmp .le. radc*eptd1)) then
             s1=tdtmp*tdtmp*tdtmp*tdtmp*tdtmp-tw05
             s1=8.356d-8*rdtmp*s1
             eptd2=s1/(5.6704d-5*s)
            endif
           elseif (jtheps .eq. 2) then
            eptd2=eptd
           elseif (jtheps .eq. 3) then
            call psi_rdmesh_int(rdtmp)
            call psi_tdemesh_int(tdtmp)
            call psi_mesh_theps
            if (theps .le. 0) then
             print *, '***Error: wrong Mie emissivity',theps
c             theps=eptd1
c             stop
             return
            endif
            eptd2=theps
           else
            print *, '***Warning: wrong emissivity selector',jtheps
c            stop
            return
           endif
          else
           s=0.d0
           eptd2=0.d0
          endif
          s=5.6704d-5*eptd2*s
c=

          qdc=qdc+s
          lttmp=s*pi4*rdtmp*rdtmp*1.d-8

          s=1.d3*mol1*(dgamm-(xgam+gsub)*md)
     *     *(tdtmp*cp+xltmp*xmeltene)
c          s=1.d3*mol1*(dgamm-(xgam+gsub)*md)*xltmp*xmeltene
          qdc=qdc+s

          if ((Tmmt    .lt. tdtmp) .and.
     *        (tdtmp   .lt. Tmpt ) .and.
     *        (canmelt .gt. 0    )) then
           s=1.d3*rdtmp*ro*xmeltene
           dtdx=3.d0*dt*(qdh-qdc)/(s*dtred)
           xltmp=xltmp+dtdx
           dtds=(xltmp-0.5d0)*2.d0*deltat+Tmm0-tdtmp
           tdtmp=tdtmp+dtds
           xltmp=min(max(0.d0,xltmp),1.0d0)
          else
           s=1.d3*rdtmp*ro*cp
           dtds=3.0d0*dt*(qdh-qdc)/(s*dtred)
           tdtmp=tdtmp+dtds
           if(tdtmp   .le. Tmmt ) xltmp=0.d0
           if(Tmpt    .le. tdtmp) xltmp=1.d0
           if(canmelt .eq. 0    ) xltmp=0.d0
          endif

          dmsd=dt/dtred*(dgamm-(xgam+gsub)*md)*mol1
     *        *pi4*rdtmp*rdtmp*1.d-12
          msdtmp=msdtmp+dmsd
          if (tdtmp .le. 0.d0) then
c           print *, '*** Warning: loss of temperature precision!'
           tdtmp=tdtmp-dtds
           msdtmp=1.d-99
           goto 11
          endif
          if (msdtmp .le. 0.d0) then
c           print *, '*** Warning: loss of mass precision!'
           msdtmp=1.d-99
           goto 11
          endif
          rdtmp=1.d4*dexp(dlog(3.d0*msdtmp/(pi4*ro))*onethird)
10       continue
11       tdd=tdtmp
         msdd=msdtmp
         xliqudd=xltmp
         lghtd=lttmp
        else
         if (jstt .gt. 0) then
          tmp=wwe2*dt
          nsts (ix,iy) = nsts (ix,iy)+tmp
          rsts (ix,iy) = rsts (ix,iy)+rd    *tmp
          msts (ix,iy) = msts (ix,iy)+msd   *tmp
          tdsts(ix,iy) = tdsts(ix,iy)+td    *tmp
          qsts (ix,iy) = qsts (ix,iy)+zd    *tmp
          ests (ix,iy) = ests (ix,iy)+ekin  *tmp
          vrsts(ix,iy) = vrsts(ix,iy)+vr    *tmp
          vzsts(ix,iy) = vzsts(ix,iy)+vz    *tmp
          vtsts(ix,iy) = vtsts(ix,iy)+vt    *tmp
          vsts (ix,iy) = vsts (ix,iy)+vabs  *tmp
          xlsts(ix,iy) = xlsts(ix,iy)+xliqud*tmp
          ltsts(ix,iy) = ltsts(ix,iy)+lght  *tmp
c --- gather statistics for distribution functions
          j=idint((dlog(rd)/aln10-minlogrd)/dlogrd)
          if (j .le. 0) then
           rdfnc (ix,iy,0)=rdfnc (ix,iy,0)+tmp
           rdfnct(      0)=rdfnct(      0)+tmp
          elseif (j .lt. frdim) then
           rdfnc (ix,iy,j)=rdfnc (ix,iy,j)+tmp
           rdfnct(      j)=rdfnct(      j)+tmp
          else
           rdfnc (ix,iy,frdim)=rdfnc (ix,iy,frdim)+tmp
           rdfnct(      frdim)=rdfnct(      frdim)+tmp
          endif
          j=idint(vabs/dvd)
          if ((j .lt. fvdim) .and. (j .ge. 0)) then
           vdfnc (ix,iy,j)=vdfnc (ix,iy,j)+tmp
           vdfnct(      j)=vdfnct(      j)+tmp
          else
           vdfnc (ix,iy,fvdim)=vdfnc (ix,iy,fvdim)+tmp
           vdfnct(      fvdim)=vdfnct(      fvdim)+tmp
          endif
          vtrs=vt+0.5d0*maxvd
          if (vtrs.ge.0.d0) then
           j=idint(vtrs/dvd)
           if ((j .lt. fvdim) .and. (j .ge. 0)) then
            vdfnctr (ix,iy,j)=vdfnctr (ix,iy,j)+tmp
            vdfnctrt(      j)=vdfnctrt(      j)+tmp
           else
            vdfnctr (ix,iy,fvdim)=vdfnctr (ix,iy,fvdim)+tmp
            vdfnctrt(      fvdim)=vdfnctrt(      fvdim)+tmp
           endif
          else
           vdfnctr (ix,iy,0)=vdfnctr (ix,iy,0)+tmp
           vdfnctrt(      0)=vdfnctrt(      0)+tmp
          endif
          j=idint(vpol/dvd)
          if ((j .lt. fvdim) .and. (j .ge. 0)) then
           vdfncpl (ix,iy,j)=vdfncpl (ix,iy,j)+tmp
           vdfncplt(      j)=vdfncplt(      j)+tmp
          else
           vdfncpl (ix,iy,fvdim)=vdfncpl (ix,iy,fvdim)+tmp
           vdfncplt(      fvdim)=vdfncplt(      fvdim)+tmp
          endif
         endif
        endif
c-----------------------------------------------------------------------
        call psi_tdmesh_int(tdd)
        call psi_tdmesh_rocp
        ro=rotd
        cp=ksicp*cptd
        rdd=1.d4*dexp(dlog(3.d0*msdd/(pi4*ro))*onethird)
        wtt=wt*msdd/msd

        fpdd=gamt
        fmdd=xgam+gsub
        qpdd=qdh
        qmdd=qdc

c***********************************************************************
c     === DUST MOTION
c***********************************************************************
c    - dust surface area,cm2
       arad=pi4*1.d-8*rd*rd
c
c    - forces, dyne=g*cm/s^2
       fz = 0.d0
       fr = 0.d0
       ft = 0.d0
c
       ffz = ffz*ksifp
       ffr = ffr*ksifp
       fft = fft*ksifp

       fz = fz+ffz
       fr = fr+ffr
       ft = ft+fft

       if (jimp .gt. 0) then
        do 90, m=1,nimpz
         ffzz(m) = ffzz(m)*ksifp
         ffrz(m) = ffrz(m)*ksifp
         fftz(m) = fftz(m)*ksifp
         fz = fz+ffzz(m)
         fr = fr+ffrz(m)
         ft = ft+fftz(m)
90      continue
       endif

c    - friction force due to atoms (absorption)
       if (un .lt. 5.-3) then
        etan=8.d0*onethird*(1.d0+.2d0*un2)
       else
        etan=(1.d0+un2-.25d0/un2)*unerf*pisqrt
        etan=etan+(un+.5d0/un)*un2exp
        etan=etan/un
       endif

       s = .5d20*ksifa*etan*arad*mpnrm*gama/etaf
       ffz = s*(uza -vz)
       ffr = s*(ura -vr)
       fft = s*(utha-vt)
       fz = fz+ffz
       fr = fr+ffr
       ft = ft+fft

c    - friction force due to impurity atoms
       if (jimp .gt. 0) then
        do 100, l=1,nimp
         if (unia(l) .lt. 5.-3) then
          etan=8.d0*onethird*(1.d0+.2d0*un2ia(l))
         else
          etan=(1.d0+un2ia(l)-.25d0/un2ia(l))*unerfia(l)*pisqrt
          etan=etan+(unia(l)+.5d0/unia(l))*un2expia(l)
          etan=etan/unia(l)
         endif

         s = .5d20*ksifa*etan*arad*mianrm(l)*gamia(l)/etiaf(l)
         ffz = s*(uzia (l)-vz)
         ffr = s*(uria (l)-vr)
         fft = s*(uthia(l)-vt)
         fz  = fz+ffz
         fr  = fr+ffr
         ft  = ft+fft
100     continue
       endif

c    - electric force= e Zd E, [E]=V/m
       s=1.602176462d-14*zd
       ffz = s*ez
       ffr = s*er
       fft = s*eth
       fz = fz+ffz
       fr = fr+ffr
       ft = ft+fft

c    - gravity acceleration (ga=9.80665 m/s2)
       ffz = -msd*9.80665d2
       ffr = 0.d0
       fft = 0.d0
       fz = fz+ffz
       fr = fr+ffr
       ft = ft+fft

c    - velocities
       invmsd=1.d0/msd
       invr=1.d0/r
       dvz = dt*fz*invmsd
       dvr = dt*(vt*vt*invr+fr*invmsd)
       dvt = dt*(ft*invmsd-vt*vr*invr)
       vvz = vz+dvz
       vvr = vr+dvr
       vvt = vt+dvt

c---coordinates
       dz  = vvz*dt
       dr  = vvr*dt 
       dth = vvt*dt*invr

       zz   = z+dz
       rr   = r+dr
       thth = th+dth

       tt=t+dt

       s=dz*dz+dr*dr
       dp=dsqrt(s)
       spp=sp+dp
        
       s=s+2.d0*r*r*(1.d0-dcos(dth))
       stt=st+dsqrt(s)

       vps   = dsqrt(vvr*vvr+vvz*vvz)
       drdts = dabs(rdd-rd)/dt
       dtdts = dabs(tdd-td)/dt

       p=p+dp
       if(p .gt. popt) goto 3
       
       r   = rr
       z   = zz
       th  = thth
       vr  = vvr
       vz  = vvz
       vt  = vvt
       rd  = rdd
       td  = tdd
       msd = msdd
       wt  = wtt
       t   = tt
       sp  = spp
       st  = stt
       fpd = fpdd
       fmd = fmdd
       qpd = qpdd
       qmd = qmdd
       fid = fidd
       xliqud = xliqudd
       lght = lghtd

1     continue ! i=1,naa

3     k=kaa

c    - save trajectory into file
      call dstrjifr(k)

      re   = rr
      ze   = zz
      the  = thth
      dte  = tt-tta(1)
      wwe  = wtt

      rs  = rdd
      ms  = msdd
      tds = tdd
      qs  =-6.9446d2*rdd*lambda*te
      vs  = vvr*vvr+vvz*vvz+vvt*vvt
      es  = 0.5d0*msdd*vs
      vrs = vvr
      vzs = vvz
      vts = vvt
      vs  = dsqrt(vs)
      ros = ro
      cps = cp
      xliquds = xliqudd
      lghts = lghtd

      pr1=ra(1)
      pz1=za(1)
      pr2=rr
      pz2=zz
      kd=idstpin(pr1,pz1,pr2,pz2)
      kod=kd
      return

      end  ! dstrjp
c=======================================================================


      subroutine dstrjn(r0,z0,th0,avr,avz,avt,msd0,tdd0,wt0)

      implicit none

      real*8 r0,z0,th0,avr,avz,avt,msd0,tdd0,wt0,ro

c      Use (dustdim)
      Use (dustinp)
      Use (dustcgm)
      Use (dusttrj)
      Use (std_cns)
      Use (psitab)
      Use (psitab_s)
      Use (psiparloc)

      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     ,            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      rb  =r0
      zb  =z0
      thb =th0
      vrb =avr
      vzb =avz
      vtb =avt
      tdb =tdd
      msdb=msd0
      tdb =tdd0
      call psi_tdmesh_int(tdb)
      call psi_tdmesh_vall
cc      call psi_tdmesh_rocp
      ro=rotd
      rdb =1.d4*dexp(dlog(3.d0*msd0/(pi4*ro))*onethird)
      wtb =wt0
      tb  =tt
      spb =sp
      stb =st

      fpb    =fpdd
      fmb    =fmdd
      qpb    =qpdd
      qmb    =qmdd
      fidb   =fidd
      xliqudb=xliqudd
      lghtb  =lghtd

      end ! dstrjn
c=======================================================================
