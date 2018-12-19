      subroutine dsttrj

      implicit none

c      Use (dustdim)
      Use (dustinp)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)
      Use (dustout)
      Use (std_cns)
      
      integer idstsrs,idstsrst,idstcgp,idstwrz,idstinc,idstpin
      integer n,nn,m,mm,ix,iy,ixx,iyy
      integer ixb,iyb
      integer tp1,tp2,tp3
      integer k,kod,ker,kfin,isd,kcx
      integer inc,cellnum, icell
      integer flg
      real*8 wes,wess,wes2,wess2,tcol,t,tt,dt,r,z,th,rr,zz,thth
      real*8 avr,avz,avt,vr,vz,vt,v,vnorm
      real*8 vrm,vzm,vtm,amu
      real*8 vrd,vzd,vtd,vdm,vdmr,vdmz,vdmt
      real*8 vnl,vnlr,vnlz,vnlt,vinv
      real*8 vbrp,vrbrp,vzbrp,vtbrp
      real*8 msdd,tdd
      real*8 s,ss,g
      real*8 vsml,betbrp,nsph,gamma
      real*8 enr,esh,fish,albw,mlbw,rnv,mrv,tc
      real*8 ficrit,fitmp
      real*8 dt1,dt2,dt3,dlpmax,dtdmax,drdrel,minps,maxtime
      real*8 rs,ms,tds,qs,es,vrs,vzs,vts,vs,ros,cps,vps,drdts,dtdts,
     *       xliquds,lghts
      real*8 vsndb,vpblb,vyblb,vrblb,vzblb,vtblb,taublb
      real*8 vrbrel,vzbrel,vtbrel,vpbrel,vbrel,ub
      real*8 fblb,fbr,fbz,fbt
      real*8 bforce
      external bforce
      common /afterstep/ rs,ms,tds,qs,es,vrs,vzs,vts,vs,ros,cps,vps,
     *                   drdts,dtdts,xliquds,lghts

      data vsml   /1.d-90/  !small velocity
      data betbrp /1.d-00/  !critical dust potential coefficient
      data nsph   /1.d+00/  !coefficient of deviation of dust form from a sphere
      data gamma  /1.d+00/  !sound process thermodinamical constant
      data dlpmax /1.d-02/  !spatial resolution scale [cm]
      data drdrel /1.d-02/  !max relative dust radius deviation on timestep
      data dtdmax /1.d+01/  !max dust temperature deviation on timestep [K]
      data maxtime/1.d+00/  !max dust survival time [s]

      g=0.d0

      n=1         ! index of injected dust particle
      m=1         ! index of dust particle in the cell
      nn=1        ! index of dust particle in size distribution
      mm=1        ! index of the cell
      
      dstflux=0.d0

c    = loop over test particles ========================================
103   continue
      dt=stp0
      call dstscll(mm,ix,iy,tp1,tp2,tp3)
      ixb=ix
      iyb=iy

      if (jstt .gt. 0) then
       if (wes0(ix,iy) .eq. 0.d0) then
        m=1
        nn=1
        mm=mm+1
        goto 103
       endif
       wess =wes0(ix,iy)*dgdist(nn)/nnispcl(ix,iy)
       wess2=ndst(ix,iy)*dndist(nn)/nnispcl(ix,iy)
      else
       wess =1.d0
       wess2=1.d0
      endif

      dstflux=dstflux+wess2
      winj=winj+wess
      ker=0
      kfin=0
      t=0.d0
      call std_rndd5x(g)
      tcol=-tmfp*dlog(g)

c    - dust particle from source
      th=0.d0
      if (jstt .gt. 0) then
       ker=idstsrst(m,ix,iy,r,z,avr,avz,avt)
      else
       ker=idstsrs(ix,iy,r,z,avr,avz,avt) ! redefines ix, iy
      endif
      if (ker .gt. 0) then
       print *, '***Error of particle source ',ker
       kfin=-3         ! source geometry error
       goto 101
      endif

c    - parameters in the cell at the dust source
      call dstcgz(ix,iy)

c    - start the dust trajectory and initialize dust parameters
      call dstrjs(ixb,iyb,r,z,th,avr,avz,avt,wess,nn)
      v=vp0
      
c    - solving motion equations to poloidal distance pp or number of steps naa
100   wes =wess
      wes2=wess2
      call dstrjp(kod,ix,iy,rr,zz,thth,dt,wess,wess2)
      
c    - time to scattering collision
      tt=t+dt

c    - flight direction vector
      v=vs
      if(v .lt. vsml) then
       v=vsml
       s=0.d0
       call std_rndd4x(g) ! unity izotropic vector
       s=twopi*g
       avr=dcos(s)
       avz=dsin(s)
       call std_rndd4x(g)
       g=twopi*g
       s=dcos(g)
       ss=dsin(g)
       avt=s
       avr=avr*ss
       avz=avz*ss
       if(vrs .lt. 0.d0) then ! keep the same direction
        if(avr .gt. 0.d0) avr=-avr
       endif
       if(vzs .lt. 0.d0) then
        if(avz .gt. 0.d0) avz=-avz
       endif
       if(vts .lt. 0.d0) then
        if(avt .gt. 0.d0) avt=-avt
       endif
      else
       avr=vrs/v
       avz=vzs/v
       avt=vts/v
      endif

      k=0
c    - transition in the mesh ------------------------------------------
      if(kod .gt. 0) then
       isd=kod ! cell side

c      new cell indexes and transition code
       ixx=ixcrs(ix,iy,isd)
       iyy=iycrs(ix,iy,isd)
       k=kcrs(ix,iy,isd)
       
c      volume statistics
       s=wes-wess
       if (jstt .gt. 0) then
        xsts(ix,iy)=xsts(ix,iy)+1.d0
        tsts(ix,iy)=tsts(ix,iy)+dt
        wsts(ix,iy)=wsts(ix,iy)+s
       endif
       wabs=wabs+s

c      small size survival
       if(rs .lt. minrd) then
        if (jstt .gt. 0) then
         wsts(ix,iy)=wsts(ix,iy)+wess
        endif
        wabs=wabs+wess
        kfin=0 ! the particle is rejected
        goto 101
       endif

c      too long time survival
       if(t .gt. maxtime) then
        werr=werr+wess
        kfin=0 ! the particle is rejected
        goto 101
       endif

c --- parsing transision type
       if (k .eq. 3) then
c    -- pass through separatrix toward core
	if (jstt .gt. 0) then
         xsptr (ix)=xsptr (ix)+1.d0
         wsptr (ix)=wsptr (ix)+wess
         gsptr (ix)=gsptr (ix)+wess2
         rsptr (ix)=rsptr (ix)+rs     *wess2
         msptr (ix)=msptr (ix)+ms     *wess2
         tdsptr(ix)=tdsptr(ix)+tds    *wess2
         qsptr (ix)=qsptr (ix)+qs     *wess2
         esptr (ix)=esptr (ix)+es     *wess2
         vrsptr(ix)=vrsptr(ix)+vrs    *wess2
         vzsptr(ix)=vzsptr(ix)+vzs    *wess2
         vtsptr(ix)=vtsptr(ix)+vts    *wess2
         vsptr (ix)=vsptr (ix)+vs     *wess2
         xlsptr(ix)=xlsptr(ix)+xliquds*wess2
         ltsptr(ix)=ltsptr(ix)+lghts  *wess2
        endif
        wspr=wspr+wess
       endif

       if ((k .eq. 7) .or. (k .eq. 71)) then
c   --- pass from inner to outer bottom divertor
        pasiob=pasiob+wess
       endif
       if ((k .eq. 8) .or. (k .eq. 81)) then
c   --- pass from outer to inner bottom divertor
        pasoib=pasoib+wess
       endif
       if (k .eq. 72) then
c   --- pass from inner to outer top divertor
        pasiot=pasiot+wess
       endif
       if (k .eq. 82) then
c   --- pass from outer to inner top divertor
        pasoit=pasoit+wess
       endif

       if ((k .eq.  0) .or.
     *     (k .eq.  3) .or.
     *     (k .eq.  7) .or.
     *     (k .eq. 71) .or.
     *     (k .eq. 72) .or.
     *     (k .eq.  8) .or.
     *     (k .eq. 81) .or.
     *     (k .eq. 82)) then
c   --- new cell
        ker=idstcgp(ixx,iyy) ! geometry parameters of new cell
        if(ker .gt. 0) then
         kfin=-2 ! geometry error
         goto 101
        endif
        ix  =ixx
        iy  =iyy
        vr  =vrs
        vz  =vzs
        vt  =vts
        msdd=ms
        tdd =tds
c   --- new plasma parameters
        call dstcgz(ix,iy)
       endif

       if(k .eq. 2) then
c    ---entering the core
        if (jstt .gt. 0) then
         xcore (ix)=xcore (ix)+1.d0
         wcore (ix)=wcore (ix)+wess
         gcore (ix)=gcore (ix)+wess2
         rcore (ix)=rcore (ix)+rs     *wess2
         mcore (ix)=mcore (ix)+ms     *wess2
         tdcore(ix)=tdcore(ix)+tds    *wess2
         qcore (ix)=qcore (ix)+qs     *wess2
         ecore (ix)=ecore (ix)+es     *wess2
         vrcore(ix)=vrcore(ix)+vrs    *wess2
         vzcore(ix)=vzcore(ix)+vzs    *wess2
         vtcore(ix)=vtcore(ix)+vts    *wess2
         vcore (ix)=vcore (ix)+vs     *wess2
         xlcore(ix)=xlcore(ix)+xliquds*wess2
         ltcore(ix)=ltcore(ix)+lghts  *wess2
        endif
        wcor=wcor+wess
        kfin=1
        goto 101
       endif

       if ((k .eq.  1) .or.
     *     (k .eq. 11) .or.
     *     (k .eq. 12) .or.
     *     (k .eq.  4) .or.
     *     (k .eq. 41) .or.
     *     (k .eq. 42) .or.
     *     (k .eq.  5) .or.
     *     (k .eq. 51) .or.
     *     (k .eq. 52) .or.
     *     (k .eq.  6) .or.
     *     (k .eq. 61) .or.
     *     (k .eq. 62)) then	   ! PFCs
        vnorm=vrs*vnr+vzs*vnz+vts*vnt
        enr=.5d0*ms*vnorm*vnorm                        ! kinetic energy of dust particle in normal to the wall direction, erg
        fish=.5d0*dlog(2*pi0*me/mp*(1.d0+gamma*ti/te)) ! sheath potential drop, e*fi/Te
        esh=qs*fish*te*ev2erg                          ! sheath potential barrier for the dust, erg, qs=Qd/e
        if (enr .gt. esh) then	                       ! collision with a wall
         albw=albws*(1.d0-xliquds)+albwl*xliquds
         mlbw=mlbws*(1.d0-xliquds)+mlbwl*xliquds
         rnv =rnvs *(1.d0-xliquds)+rnvl *xliquds
         mrv =mrvs *(1.d0-xliquds)+mrvl *xliquds
         tc  =tcs  *(1.d0-xliquds)+tcl  *xliquds
         s=wess*albw*mlbw
         
c --- collision with wall
         if ((k .eq.  1) .or.
     *       (k .eq. 11) .or.
     *       (k .eq. 12)) then
          if ((rs .lt. minrd) .or. (s .eq. 0.d0)) then
           if (jstt .gt. 0) then
            xwall (ix)=xwall (ix)+1.d0
            wwall (ix)=wwall (ix)+wess
            gwall (ix)=gwall (ix)+wess2
            rwall (ix)=rwall (ix)+rs     *wess2
            mwall (ix)=mwall (ix)+ms     *wess2
            tdwall(ix)=tdwall(ix)+tds    *wess2
            qwall (ix)=qwall (ix)+qs     *wess2
            ewall (ix)=ewall (ix)+es     *wess2
            vrwall(ix)=vrwall(ix)+vrs    *wess2
            vzwall(ix)=vzwall(ix)+vzs    *wess2
            vtwall(ix)=vtwall(ix)+vts    *wess2
            vwall (ix)=vwall (ix)+vs     *wess2
            xlwall(ix)=xlwall(ix)+xliquds*wess2
            ltwall(ix)=ltwall(ix)+lghts  *wess2
           endif
           if (k .eq. 1) then
            if (ix .le. ixtop) then  ! warning: ixtop is defined only when k=1
             depwi=depwi+wess
             flxwi=flxwi+wess
            else
             depwo=depwo+wess
             flxwo=flxwo+wess
            endif
           elseif (k .eq. 11) then
            depwi=depwi+wess
            flxwi=flxwi+wess
           else
            depwo=depwo+wess
            flxwo=flxwo+wess
           endif
           kfin=0
           goto 101
          endif
          if (jstt .gt. 0) then
           xwall (ix)=xwall (ix)+1.d0
           wwall (ix)=wwall (ix)+wess-s
           gwall (ix)=gwall (ix)+wess2
           rwall (ix)=rwall (ix)+rs     *wess2
           mwall (ix)=mwall (ix)+ms     *wess2
           tdwall(ix)=tdwall(ix)+tds    *wess2
           qwall (ix)=qwall (ix)+qs     *wess2
           ewall (ix)=ewall (ix)+es     *wess2
           vrwall(ix)=vrwall(ix)+vrs    *wess2
           vzwall(ix)=vzwall(ix)+vzs    *wess2
           vtwall(ix)=vtwall(ix)+vts    *wess2
           vwall (ix)=vwall (ix)+vs     *wess2
           xlwall(ix)=xlwall(ix)+xliquds*wess2
           ltwall(ix)=ltwall(ix)+lghts  *wess2
          endif
          if (k .eq. 1) then
           if (ix .le. ixtop) then   ! warning: ixtop is defined only when k=1
            depwi=depwi+wess-s
            flxwi=flxwi+wess
            rflwi=rflwi+s
           else
            depwo=depwo+wess-s
            flxwo=flxwo+wess
            rflwo=rflwo+s
           endif
          elseif (k .eq. 11) then
           depwi=depwi+wess-s
           flxwi=flxwi+wess
           rflwi=rflwi+s
          else
           depwo=depwo+wess-s
           flxwo=flxwo+wess
           rflwo=rflwo+s
          endif
         endif
         
c --- collision with PF region wall
         if ((k .eq.  4) .or.
     *       (k .eq. 41) .or.
     *       (k .eq. 42)) then
          if((rs .lt. minrd) .or. (s .eq. 0.d0)) then
           if (jstt .gt. 0) then
            xpfwl (ix)=xpfwl (ix)+1.d0
            wpfwl (ix)=wpfwl (ix)+wess
            gpfwl (ix)=gpfwl (ix)+wess2
            rpfwl (ix)=rpfwl (ix)+rs     *wess2
            mpfwl (ix)=mpfwl (ix)+ms     *wess2
            tdpfwl(ix)=tdpfwl(ix)+tds    *wess2
            qpfwl (ix)=qpfwl (ix)+qs     *wess2
            epfwl (ix)=epfwl (ix)+es     *wess2
            vrpfwl(ix)=vrpfwl(ix)+vrs    *wess2
            vzpfwl(ix)=vzpfwl(ix)+vzs    *wess2
            vtpfwl(ix)=vtpfwl(ix)+vts    *wess2
            vpfwl (ix)=vpfwl (ix)+vs     *wess2
            xlpfwl(ix)=xlpfwl(ix)+xliquds*wess2
            ltpfwl(ix)=ltpfwl(ix)+lghts  *wess2
           endif
           if ((k .eq. 4) .or. (k .eq. 41)) then
            if (ix .le. ixpt1(1)) then
             deppib=deppib+wess
             flxpib=flxpib+wess
            else
             deppob=deppob+wess
             flxpob=flxpob+wess
            endif
           else
            if (ix .le. ixrb1) then
             deppit=deppit+wess
             flxpit=flxpit+wess
            else
             deppot=deppot+wess
             flxpot=flxpot+wess
            endif
           endif
           kfin=0
           goto 101
          endif
          if (jstt .gt. 0) then
           xpfwl (ix)=xpfwl (ix)+1.d0
           wpfwl (ix)=wpfwl (ix)+wess-s
           gpfwl (ix)=gpfwl (ix)+wess2
           rpfwl (ix)=rpfwl (ix)+rs     *wess2
           mpfwl (ix)=mpfwl (ix)+ms     *wess2
           tdpfwl(ix)=tdpfwl(ix)+tds    *wess2
           qpfwl (ix)=qpfwl (ix)+qs     *wess2
           epfwl (ix)=epfwl (ix)+es     *wess2
           vrpfwl(ix)=vrpfwl(ix)+vrs    *wess2
           vzpfwl(ix)=vzpfwl(ix)+vzs    *wess2
           vtpfwl(ix)=vtpfwl(ix)+vts    *wess2
           vpfwl (ix)=vpfwl (ix)+vs     *wess2
           xlpfwl(ix)=xlpfwl(ix)+xliquds*wess2
           ltpfwl(ix)=ltpfwl(ix)+lghts  *wess2
          endif
          if ((k .eq. 4) .or. (k .eq. 41)) then
           if (ix .le. ixpt1(1)) then
            deppib=deppib+wess-s
            flxpib=flxpib+wess
            rflpib=rflpib+s
           else
            deppob=deppob+wess-s
            flxpob=flxpob+wess
            rflpob=rflpob+s
           endif
          else
           if (ix .le. ixrb1) then
            deppit=deppit+wess-s
            flxpit=flxpit+wess
            rflpit=rflpit+s
           else
            deppot=deppot+wess-s
            flxpot=flxpot+wess
            rflpot=rflpot+s
           endif
          endif
         endif
          
c --- collision with inner divertor plates
         if ((k .eq.  5) .or. (k .eq. 51)) then
          if((rs .lt. minrd) .or. (s .eq. 0.d0)) then
           if (jstt .gt. 0) then
            xdivib (iy)=xdivib (iy)+1.d0
            wdivib (iy)=wdivib (iy)+wess
            gdivib (iy)=gdivib (iy)+wess2
            rdivib (iy)=rdivib (iy)+rs     *wess2
            mdivib (iy)=mdivib (iy)+ms     *wess2
            tddivib(iy)=tddivib(iy)+tds    *wess2
            qdivib (iy)=qdivib (iy)+qs     *wess2
            edivib (iy)=edivib (iy)+es     *wess2
            vrdivib(iy)=vrdivib(iy)+vrs    *wess2
            vzdivib(iy)=vzdivib(iy)+vzs    *wess2
            vtdivib(iy)=vtdivib(iy)+vts    *wess2
            vdivib (iy)=vdivib (iy)+vs     *wess2
            xldivib(iy)=xldivib(iy)+xliquds*wess2
            ltdivib(iy)=ltdivib(iy)+lghts  *wess2
           endif
           depdib=depdib+wess
           flxdib=flxdib+wess
           kfin=0
           goto 101
          endif
          if (jstt .gt. 0) then
           xdivib (iy)=xdivib (iy)+1.d0
           wdivib (iy)=wdivib (iy)+wess-s
           gdivib (iy)=gdivib (iy)+wess2
           rdivib (iy)=rdivib (iy)+rs     *wess2
           mdivib (iy)=mdivib (iy)+ms     *wess2
           tddivib(iy)=tddivib(iy)+tds    *wess2
           qdivib (iy)=qdivib (iy)+qs     *wess2
           edivib (iy)=edivib (iy)+es     *wess2
           vrdivib(iy)=vrdivib(iy)+vrs    *wess2
           vzdivib(iy)=vzdivib(iy)+vzs    *wess2
           vtdivib(iy)=vtdivib(iy)+vts    *wess2
           vdivib (iy)=vdivib (iy)+vs     *wess2
           xldivib(iy)=xldivib(iy)+xliquds*wess2
           ltdivib(iy)=ltdivib(iy)+lghts  *wess2
          endif
          depdib=depdib+wess-s
          flxdib=flxdib+wess
          rfldib=rfldib+s
         endif
         if (k .eq. 52) then
          if((rs .lt. minrd) .or. (s .eq. 0.d0)) then
           if (jstt .gt. 0) then
            xdivit (iy)=xdivit (iy)+1.d0
            wdivit (iy)=wdivit (iy)+wess
            gdivit (iy)=gdivit (iy)+wess2
            rdivit (iy)=rdivit (iy)+rs     *wess2
            mdivit (iy)=mdivit (iy)+ms     *wess2
            tddivit(iy)=tddivit(iy)+tds    *wess2
            qdivit (iy)=qdivit (iy)+qs     *wess2
            edivit (iy)=edivit (iy)+es     *wess2
            vrdivit(iy)=vrdivit(iy)+vrs    *wess2
            vzdivit(iy)=vzdivit(iy)+vzs    *wess2
            vtdivit(iy)=vtdivit(iy)+vts    *wess2
            vdivit (iy)=vdivit (iy)+vs     *wess2
            xldivit(iy)=xldivit(iy)+xliquds*wess2
            ltdivit(iy)=ltdivit(iy)+lghts  *wess2
           endif
           depdit=depdit+wess
           flxdit=flxdit+wess
           kfin=0
           goto 101
          endif
          if (jstt .gt. 0) then
           xdivit (iy)=xdivit (iy)+1.d0
           wdivit (iy)=wdivit (iy)+wess-s
           gdivit (iy)=gdivit (iy)+wess2
           rdivit (iy)=rdivit (iy)+rs     *wess2
           mdivit (iy)=mdivit (iy)+ms     *wess2
           tddivit(iy)=tddivit(iy)+tds    *wess2
           qdivit (iy)=qdivit (iy)+qs     *wess2
           edivit (iy)=edivit (iy)+es     *wess2
           vrdivit(iy)=vrdivit(iy)+vrs    *wess2
           vzdivit(iy)=vzdivit(iy)+vzs    *wess2
           vtdivit(iy)=vtdivit(iy)+vts    *wess2
           vdivit (iy)=vdivit (iy)+vs     *wess2
           xldivit(iy)=xldivit(iy)+xliquds*wess2
           ltdivit(iy)=ltdivit(iy)+lghts  *wess2
          endif
          depdit=depdit+wess-s
          flxdit=flxdit+wess
          rfldit=rfldit+s
         endif

c --- collision with outer divertor plates
         if ((k .eq. 6) .or. (k .eq. 61)) then
          if((rs .lt. minrd) .or. (s .eq. 0.d0)) then
           if (jstt .gt. 0) then
            xdivob (iy)=xdivob (iy)+1.d0
            wdivob (iy)=wdivob (iy)+wess
            gdivob (iy)=gdivob (iy)+wess2
            rdivob (iy)=rdivob (iy)+rs     *wess2
            mdivob (iy)=mdivob (iy)+ms     *wess2
            tddivob(iy)=tddivob(iy)+tds    *wess2
            qdivob (iy)=qdivob (iy)+qs     *wess2
            edivob (iy)=edivob (iy)+es     *wess2
            vrdivob(iy)=vrdivob(iy)+vrs    *wess2
            vzdivob(iy)=vzdivob(iy)+vzs    *wess2
            vtdivob(iy)=vtdivob(iy)+vts    *wess2
            vdivob (iy)=vdivob (iy)+vs     *wess2
            xldivob(iy)=xldivob(iy)+xliquds*wess2
            ltdivob(iy)=ltdivob(iy)+lghts  *wess2
           endif
           depdob=depdob+wess
           flxdob=flxdob+wess
           kfin=0
           goto 101
          endif
          if (jstt .gt. 0) then
           xdivob (iy)=xdivob (iy)+1.d0
           wdivob (iy)=wdivob (iy)+wess-s
           gdivob (iy)=gdivob (iy)+wess2
           rdivob (iy)=rdivob (iy)+rs     *wess2
           mdivob (iy)=mdivob (iy)+ms     *wess2
           tddivob(iy)=tddivob(iy)+tds    *wess2
           qdivob (iy)=qdivob (iy)+qs     *wess2
           edivob (iy)=edivob (iy)+es     *wess2
           vrdivob(iy)=vrdivob(iy)+vrs    *wess2
           vzdivob(iy)=vzdivob(iy)+vzs    *wess2
           vtdivob(iy)=vtdivob(iy)+vts    *wess2
           vdivob (iy)=vdivob (iy)+vs     *wess2
           xldivob(iy)=xldivob(iy)+xliquds*wess2
           ltdivob(iy)=ltdivob(iy)+lghts  *wess2
          endif
          depdob=depdob+wess-s
          flxdob=flxdob+wess
          rfldob=rfldob+s
         endif
         if (k .eq. 62) then
          if((rs .lt. minrd) .or. (s .eq. 0.d0)) then
           if (jstt .gt. 0) then
            xdivot (iy)=xdivot (iy)+1.d0
            wdivot (iy)=wdivot (iy)+wess
            gdivot (iy)=gdivot (iy)+wess2
            rdivot (iy)=rdivot (iy)+rs     *wess2
            mdivot (iy)=mdivot (iy)+ms     *wess2
            tddivot(iy)=tddivot(iy)+tds    *wess2
            qdivot (iy)=qdivot (iy)+qs     *wess2
            edivot (iy)=edivot (iy)+es     *wess2
            vrdivot(iy)=vrdivot(iy)+vrs    *wess2
            vzdivot(iy)=vzdivot(iy)+vzs    *wess2
            vtdivot(iy)=vtdivot(iy)+vts    *wess2
            vdivot (iy)=vdivot (iy)+vs     *wess2
            xldivot(iy)=xldivot(iy)+xliquds*wess2
            ltdivot(iy)=ltdivot(iy)+lghts  *wess2
           endif
           depdot=depdot+wess
           flxdot=flxdot+wess
           kfin=0
           goto 101
          endif
          if (jstt .gt. 0) then
           xdivot (iy)=xdivot (iy)+1.d0
           wdivot (iy)=wdivot (iy)+wess-s
           gdivot (iy)=gdivot (iy)+wess2
           rdivot (iy)=rdivot (iy)+rs     *wess2
           mdivot (iy)=mdivot (iy)+ms     *wess2
           tddivot(iy)=tddivot(iy)+tds    *wess2
           qdivot (iy)=qdivot (iy)+qs     *wess2
           edivot (iy)=edivot (iy)+es     *wess2
           vrdivot(iy)=vrdivot(iy)+vrs    *wess2
           vzdivot(iy)=vzdivot(iy)+vzs    *wess2
           vtdivot(iy)=vtdivot(iy)+vts    *wess2
           vdivot (iy)=vdivot (iy)+vs     *wess2
           xldivot(iy)=xldivot(iy)+xliquds*wess2
           ltdivot(iy)=ltdivot(iy)+lghts  *wess2
          endif
          depdot=depdot+wess-s
          flxdot=flxdot+wess
          rfldot=rfldot+s
         endif
         
         wess=s
         wess2=wess2*albw
         msdd=ms*mlbw

c    --- reflection from boundary at intersection point (rrt,zzt)
         vr=vrs
         vz=vzs
         vt=vts

1004     vinv=1.d0/dsqrt(vr*vr+vz*vz+vt*vt)
         vr=vr*vinv
         vz=vz*vinv
         vt=vt*vinv

c    --- finding reflection vector for elastic collision (restitution = 1)
         vrd=vr
         vzd=vz
         vtd=vt
         vrm=vr
         vzm=vz
         vtm=vt
         amu=vrm*vnr+vzm*vnz+vtm*vnt
         call rflmir(vrm,vzm,vtm, amu, vnr,vnz,vnt)
         call rfldif(vrd,vzd,vtd, vnr,vnz,vnt)
         vrd=vrd*(1.d0-mrv)
         vzd=vzd*(1.d0-mrv)
         vtd=vtd*(1.d0-mrv)
         vrm=vrm*mrv
         vzm=vzm*mrv
         vtm=vtm*mrv
         amu=mrv*(1.d0-mrv)-(vrm*vrd+vzm*vzd+vtm*vtd)
         vdm=1.d0/dsqrt(1.d0-2.d0*amu)
         vdmr=(vrd+vrm)*vdm
         vdmz=(vzd+vzm)*vdm
         vdmt=(vtd+vtm)*vdm
c    --- local normal vector at strike point
         vnlr=vdmr-vr
         vnlz=vdmz-vz
         vnlt=vdmt-vt
         vnl=1.d0/dsqrt(vnlr*vnlr+vnlz*vnlz+vnlt*vnlt)
         vnlr=vnlr*vnl
         vnlz=vnlz*vnl
         vnlt=vnlt*vnl
c    --- local plastic reflection ('rnv' is the coefficient of restitution)
         amu=(1.d0-rnv)*(vdmr*vnlr+vdmz*vnlz+vdmt*vnlt)
         vr=vdmr-vnlr*amu
         vz=vdmz-vnlz*amu
         vt=vdmt-vnlt*amu
         vr=vr*v
         vz=vz*v
         vt=vt*v
c    --- check outward reflection
         amu=vr*vnr+vz*vnz+vt*vnt
         if (amu .lt. 0.d0) goto 1004
         
         tdd=tds+0.5d0*tc*(vrs*vrs+vzs*vzs+vts*vts-vr*vr-vz*vz-vt*vt)
     *                    /cps
         tdd=tds

         rr=rrtin
         zz=zztin

        else
c    --- mirror reflection inside sheath (rrt,zzt)
         vr=avr
         vz=avz
         vt=avt
         amu=(vr*vnr+vz*vnz+vt*vnt)/v
         call rflmir(vr,vz,vt, amu, vnr,vnz,vnt)
         vr=vr*v
         vz=vz*v
         vt=vt*v
         rr=rrtin
         zz=zztin
         msdd=ms
         tdd =tds
        endif
       endif	! PFCs

      elseif(kod .eq. 0) then ! same cell ------------------------------
c
c    - point in a cell
       s=wes-wess
       if (jstt .gt. 0) then
        tsts(ix,iy)=tsts(ix,iy)+dt
        wsts(ix,iy)=wsts(ix,iy)+s
       endif
       wabs=wabs+s
c
c    - small size survival
       if(rs .lt. minrd) then
        if (jstt .gt. 0) then
         wsts(ix,iy)=wsts(ix,iy)+wess
        endif
        wabs=wabs+wess
        kfin=0 ! particle rejected
        goto 101
       endif

c      too long time survival
       if(t .gt. maxtime) then
        werr=werr+wess
        kfin=0 ! the particle is rejected
        goto 101
       endif

      else ! kod<0 -----------------------------------------------------
c    - error in geometry while transition thru boundary
       write(*,1001) 'geometry error=', ker,ix,iy
1001   format(a11,3I4)
       write (*,1002) r,z,rr,zz
1002   format(4d25.15)
       kfin=-1
       goto 101
	 
      endif ! kod ------------------------------------------------------

c --- particle electrostatic break up
      if (jdsr .gt. 0) then
       if (xliquds .lt. 0.5d0) then
        ficrit=0.3361d0*betbrp*dsqrt(tensf)*rs     ! critical dust potential for solids, V
       else
        ficrit=4.7527d3*betbrp*dsqrt(stens*rs)     ! critical dust potential for drops, V
       endif
       fitmp=1.4400d-3*dabs(qs)/rs                 ! dust potential, Volts
       if (fitmp .ge. ficrit) then
        vbrp=2.7816d-10*fitmp*fitmp*rs/(2.d0*gsplit*ms)  ! part speed, cm/s
        call std_rndd6x(g)
        vbrp=dsqrt(vbrp*g)
        call std_rndd6x(g)
        s=twopi*g
        vrbrp=dcos(s)
        vzbrp=dsin(s)
        call std_rndd6x(g)
        g=pi0*g
        s=dcos(g)
        ss=dsin(g)
        vrbrp=vbrp*ss*vrbrp
        vzbrp=vbrp*ss*vzbrp
        vtbrp=vbrp*s
        vr=vrs+vrbrp
        vz=vzs+vzbrp
        vt=vts+vtbrp
        msdd =gsplit*ms
        tdd  =tds
        wess =2.d0*gsplit*wess
        wess2=2.d0*wess2
       endif
      endif
c
c --- collision with plasma turbulence
      if ((tt .ge. tcol) .and. (jcol .gt. 0)) then
       flg=0
       if (iy .le. iysptrx(1)) then
        flg=1
       elseif (kmesh .eq. 'dnull') then
         if ((iy .le. iysptrx(2)) .and.
     *       (ix .gt.   ixpt2(1)) .and.
     *       (ix .le.   ixpt1(2))) flg=1
       endif
       if (flg .eq. 1) then
        vr  =vrs
        vz  =vzs
        vt  =vts
        msdd=ms
        tdd =tds
       else
        call std_rndd5x(g)
        if (g .le. blf*tmfp) 
     *  then
c ----   collision with plasma turbulence
         vsndb=dsqrt(2.d0*tblb/mpnrm)
         vpblb=machb*vsndb
         vyblb=vblb
         call dstvcp(ix,iy, vpblb,vyblb, vrblb,vzblb,vtblb)
         call std_rndd5x(g)
         if (g .gt. 0.5d0) vtblb=-vtblb
         vrbrel=vrblb-vrs
         vzbrel=vzblb-vzs
         vtbrel=vtblb-vts
         vpbrel=dsqrt(vrbrel*vrbrel+vzbrel*vzbrel)
         vbrel=dsqrt(vpbrel*vpbrel+vtbrel*vtbrel)
         taublb=2.d0*rblb/vpbrel
         ub=vbrel/vsndb

         fblb=bforce(ub,rs)

         if (vbrel .gt. 0.d0) then
          fbr=fblb*vrbrel/vbrel
          fbz=fblb*vzbrel/vbrel
          fbt=fblb*vtbrel/vbrel
         else
          fbr=0.d0
          fbz=0.d0
          fbt=0.d0
         endif
 
         vr=vrs+taublb*fbr/ms
         vz=vzs+taublb*fbz/ms
         vt=vts+taublb*fbt/ms

         msdd=ms
         tdd =tds
 
         tt=0.d0
         call std_rndd5x(g)
         tcol=-tmfp*dlog(g)
        else
c     -- null collision
         vr  =vrs
         vz  =vzs
         vt  =vts
         msdd=ms
         tdd =tds
        endif
       endif
      else
c   -- no change
       vr  =vrs
       vz  =vzs
       vt  =vts
       msdd=ms
       tdd =tds
      endif
c
c --- check mesh
      inc=idstwrz(rr,zz)
      if (inc .ne. 1) then
*       print *, '***Warning: cell is dropped'
       cellnum=idstinc(rr,zz)
       iyy=cellnum/nxm
       ixx=cellnum-iyy*nxm
       flg=0
       if (kmesh .eq. 'dnull') then
        if ((ixx .gt. ixrb1) .and. (ixx .le. ixlb2)) flg=1
       endif
       if ((cellnum .le. 0) .or. (flg .eq. 1)) then
        ker=idstcgp(ix,iy)
        if(ker .gt. 0) then
         kfin=-2 ! geometry error
         goto 101
        endif
        kcx=idstpin(rrtout,zztout,rr,zz)
        if (kcx .le. 0) then
         print *, '***Error: particle is out of mesh!',rr,zz,r,z
c         stop
         return
        else
         rr=rrtin
         zz=zztin
        endif
       else
        ix=ixx
        iy=iyy
        ker=idstcgp(ix,iy)
        if(ker .gt. 0) then
         kfin=-2 ! geometry error
         goto 101
        endif
        call dstcgz(ix,iy)
       endif
      endif

      r=rr
      z=zz
      th=thth
      t=tt

      call dstrjn(r,z,th,vr,vz,vt,msdd,tdd,wess)

c --- define new time step
      minps=min(dlpmax,msize(ix,iy))
      dt1=minps/vps
      dt2=drdrel*rs/drdts
      dt3=dtdmax/dtdts
      dt=min(dt1,min(dt2,dt3))
      if (dt .lt. stp0) dt=stp0
      
c    - continue motion
      goto 100

101   if (kfin .lt. 0) then
       werr=werr+wess
       write(*,1003) '***Error of mesh',-kfin
1003   format(a16,i6)
      endif

      n=n+1
      m=m+1
      nn=nn+1
      if (nn .gt. frdim) then
       nn=1
      endif
      if (m .gt. nispcl(ixb,iyb)) then
       m=1
       nn=1
       mm=mm+1
      endif
      if(n .le. nisp) goto 103
c    = end loop over test particles ====================================

      end


      real*8 function bforce(ub,rd)
      implicit none
      
c      Use (dustdim)
      Use (dustinp)
      Use (dustcom)
      Use (std_cns)

      real*8 ub,rd      
      real*8 ub2,ub2exp,berf,coeff
      real*8 dnume1,dnume2,dnume3
      real*8 c1,c2,c3,c4
      real*8 fC,fO
      real*8 xisqr,rdnrm,g,debyee,debyes,invdbys,a,eta,b90,b902,xlng,gu
      real*8 std_erf
      external std_erf

      ub2=ub*ub
      ub2exp=dexp(-ub2)
      berf=std_erf(ub)
      rdnrm=rd*1.d-4
      coeff=2.d0*rdnrm*rdnrm*nblb*tblb

c --- coefficient of Fc
      if (lambda .ge. 0.d0) then 
       if(ub .lt. 5.d-3)then
        fC=pisqrt*onethird*(8.d0+4.d0*lambda)*ub
       else
        dnume1=ub*(ub2+.5d0+lambda)*ub2exp
        dnume2=ub2*(ub2+1.d0)-.25d0-(.5d0-ub2)*lambda
        dnume2=dnume2*pisqrt*berf
        dnume3=(dnume1+dnume2)*pisqrt
        fC=dnume3/ub2
       endif
      else
       if(ub .lt. 5.d-3)then
        fC=pisqrt*onethird*(8.d0-4.d0*lambda)*ub*dexp(lambda)
       else
        xisqr=dsqrt(-lambda)
        c3=2.d0*ub2
        c1=1.d0+c3
        c2=1.d0-c3
        c3=c1
        c4=xisqr*c2/ub
        c1=c3+c4
        c2=c3-c4
        c3=ub+xisqr
        c4=ub-xisqr
        dnume1=dexp(-c3*c3)*c1+dexp(-c4*c4)*c2
        dnume1=dnume1*pisqrt/ub
        c1=1.d0+2.d0*(lambda+ub2)
        c2=1.d0+2.d0*(lambda-ub2)
        dnume2=c1-.5d0*c2/ub2
        dnume2=pi0*dnume2*(std_erf(c3)+std_erf(c4))
        fC=.25d0*(dnume1+dnume2)
       endif
      endif

c --- coulomb logarithm
      g=1.d0/(3.d0+2.d0*ub2)
      debyee=5.873d8*dsqrt(tblb/nblb)
      invdbys=1.d0+3.d0*g
      invdbys=debyee*debyee/invdbys
      debyes=dsqrt(invdbys)
      a=rdnrm/debyes
      if (lambda .ge. 0.d0) then
       eta=(1.d0+1.40824829d0*a)*debyes
       b90=rdnrm*lambda*g
       b902=b90*b90
       xlng=(b902+eta*eta)/(b902+rdnrm*rdnrm)
       xlng=.5d0*dlog(xlng)
      else
       xlng=-dlog(a)+lambda*(0.306852819d0-lambda*9.657359d-2)
      endif
      xlng=max(0.d0,xlng)
        
c --- coefficient of Fo
      if(ub .lt. 5.d-3) then 
       gu=2.d0*onethird*ub*invpisqrt
      else 
       dnume1=.5d0*berf-ub*ub2exp*invpisqrt
       gu=dnume1/ub2
      endif
      fO=lambda*lambda*(2.d0*pi0)*xlng*gu
      
      bforce=coeff*(fC+fO)
      return
      
      end
