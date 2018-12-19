      subroutine dststn

      implicit none

c      Use (dustdim)
      Use (dustinp)
      Use (dustcom)
      Use (dustout)
      Use (std_cns)

      integer i,j,k

      winj  =0.d0
      wcor  =0.d0
      wspr  =0.d0
      wabs  =0.d0
      werr  =0.d0
      depwi =0.d0
      depwo =0.d0
      depdib=0.d0
      depdob=0.d0
      depdit=0.d0
      depdot=0.d0
      deppib=0.d0
      deppob=0.d0
      deppit=0.d0
      deppot=0.d0
      flxwi =0.d0
      flxwo =0.d0
      flxdib=0.d0
      flxdob=0.d0
      flxdit=0.d0
      flxdot=0.d0
      flxpib=0.d0
      flxpob=0.d0
      flxpit=0.d0
      flxpot=0.d0
      rflwi =0.d0
      rflwo =0.d0
      rfldib=0.d0
      rfldob=0.d0
      rfldit=0.d0
      rfldot=0.d0
      rflpib=0.d0
      rflpob=0.d0
      rflpit=0.d0
      rflpot=0.d0
      pasiob=0.d0
      pasoib=0.d0
      pasiot=0.d0
      pasoit=0.d0

      do 1  i=1,nxm
       xcore (i)=0.d0
       xsptr (i)=0.d0
       xwall (i)=0.d0
       xpfwl (i)=0.d0
       wcore (i)=0.d0
       wsptr (i)=0.d0
       wwall (i)=0.d0
       wpfwl (i)=0.d0
       gcore (i)=0.d0
       gsptr (i)=0.d0
       gwall (i)=0.d0
       gpfwl (i)=0.d0
       rcore (i)=0.d0
       rsptr (i)=0.d0
       rwall (i)=0.d0
       rpfwl (i)=0.d0
       mcore (i)=0.d0
       msptr (i)=0.d0
       mwall (i)=0.d0
       mpfwl (i)=0.d0
       tdcore(i)=0.d0
       tdsptr(i)=0.d0
       tdwall(i)=0.d0
       tdpfwl(i)=0.d0
       qcore (i)=0.d0
       qsptr (i)=0.d0
       qwall (i)=0.d0
       qpfwl (i)=0.d0
       ecore (i)=0.d0
       esptr (i)=0.d0
       ewall (i)=0.d0
       epfwl (i)=0.d0
       vrcore(i)=0.d0
       vrsptr(i)=0.d0
       vrwall(i)=0.d0
       vrpfwl(i)=0.d0
       vzcore(i)=0.d0
       vzsptr(i)=0.d0
       vzwall(i)=0.d0
       vzpfwl(i)=0.d0
       vtcore(i)=0.d0
       vtsptr(i)=0.d0
       vtwall(i)=0.d0
       vtpfwl(i)=0.d0
       vcore (i)=0.d0
       vsptr (i)=0.d0
       vwall (i)=0.d0
       vpfwl (i)=0.d0
       xlcore(i)=0.d0
       xlsptr(i)=0.d0
       xlwall(i)=0.d0
       xlpfwl(i)=0.d0
       ltcore(i)=0.d0
       ltsptr(i)=0.d0
       ltwall(i)=0.d0
       ltpfwl(i)=0.d0
       do 11 j=1,nym
        xsts (i,j)=0.d0
        tsts (i,j)=0.d0
        wsts (i,j)=0.d0
        nsts (i,j)=0.d0
        rsts (i,j)=0.d0
        msts (i,j)=0.d0
        tdsts(i,j)=0.d0
        qsts (i,j)=0.d0
        ests (i,j)=0.d0
        vrsts(i,j)=0.d0
        vzsts(i,j)=0.d0
        vtsts(i,j)=0.d0
        vsts (i,j)=0.d0
        xlsts(i,j)=0.d0
        ltsts(i,j)=0.d0
        do 12, k=0,frdim
         rdfnc(i,j,k)=0.d0
12      continue
        do 13, k=0,fvdim
         vdfnc  (i,j,k)=0.d0
         vdfnctr(i,j,k)=0.d0
         vdfncpl(i,j,k)=0.d0
13      continue
11     continue
1     continue

      do 14, k=0,frdim
       rdfnct(k)=0.d0
       rdnod (k)=dexp((minlogrd+k*dlogrd)*aln10)
14    continue
      do 15, k=0,fvdim
       vdfnct  (k)=0.d0
       vdfnctrt(k)=0.d0
       vdfncplt(k)=0.d0
       vdnod   (k)=k*dvd
       vtrdnod (k)=vdnod(k)-0.5d0*maxvd
15    continue

      do 21 j=1,nym
       xdivib (j)=0.d0
       xdivob (j)=0.d0
       wdivib (j)=0.d0
       wdivob (j)=0.d0
       gdivib (j)=0.d0
       gdivob (j)=0.d0
       rdivib (j)=0.d0
       rdivob (j)=0.d0
       mdivib (j)=0.d0
       mdivob (j)=0.d0
       tddivib(j)=0.d0
       tddivob(j)=0.d0
       qdivib (j)=0.d0
       qdivob (j)=0.d0
       edivib (j)=0.d0
       edivob (j)=0.d0
       vrdivib(j)=0.d0
       vrdivob(j)=0.d0
       vzdivib(j)=0.d0
       vzdivob(j)=0.d0
       vtdivib(j)=0.d0
       vtdivob(j)=0.d0
       vdivib (j)=0.d0
       vdivob (j)=0.d0
       xldivib(j)=0.d0
       xldivob(j)=0.d0
       ltdivib(j)=0.d0
       ltdivob(j)=0.d0
       xdivit (j)=0.d0
       xdivot (j)=0.d0
       wdivit (j)=0.d0
       wdivot (j)=0.d0
       gdivit (j)=0.d0
       gdivot (j)=0.d0
       rdivit (j)=0.d0
       rdivot (j)=0.d0
       mdivit (j)=0.d0
       mdivot (j)=0.d0
       tddivit(j)=0.d0
       tddivot(j)=0.d0
       qdivit (j)=0.d0
       qdivot (j)=0.d0
       edivit (j)=0.d0
       edivot (j)=0.d0
       vrdivit(j)=0.d0
       vrdivot(j)=0.d0
       vzdivit(j)=0.d0
       vzdivot(j)=0.d0
       vtdivit(j)=0.d0
       vtdivot(j)=0.d0
       vdivit (j)=0.d0
       vdivot (j)=0.d0
       xldivit(j)=0.d0
       xldivot(j)=0.d0
       ltdivit(j)=0.d0
       ltdivot(j)=0.d0
21    continue

      end


      subroutine dststo
c
c---this subroutine is to process 
c   and printout the output data 
c
      implicit none

      Use (dustinp)
c      Use (dustdim)
      Use (dustcom)
      Use (dustout)

      integer vlm,core,sptr,wall,pfwl
      integer divib,divob,divit,divot
      integer rdf ,vdf
      integer rdft,vdft
      integer i,j,k
      integer flg
      real*8 nrm,nstst
      data vlm,core,sptr,wall,pfwl /31,32,33,34,35/
      data divib,divob,divit,divot /36,37,38,39/
      data rdf ,vdf  /41,42/
      data rdft,vdft /43,44/

      open (vlm ,file=outdir(1:len_trim(outdir))//'vlm_stat.txt' ,err=1)
      open (core,file=outdir(1:len_trim(outdir))//'core_stat.txt',err=1)
      open (sptr,file=outdir(1:len_trim(outdir))//'sptr_stat.txt',err=1)
      open (wall,file=outdir(1:len_trim(outdir))//'wall_stat.txt',err=1)
      open (pfwl,file=outdir(1:len_trim(outdir))//'pfwl_stat.txt',err=1)
      open (divib,file=outdir(1:len_trim(outdir))//
     *                 'divib_stat.txt',err=1)
      open (divob,file=outdir(1:len_trim(outdir))//
     *                 'divob_stat.txt',err=1)
      if (kmesh .eq. 'dnull') then
       open (divit,file=outdir(1:len_trim(outdir))//
     *                  'divit_stat.txt',err=1)
       open (divot,file=outdir(1:len_trim(outdir))//
     *                  'divot_stat.txt',err=1)
      endif
      open (rdf ,file=outdir(1:len_trim(outdir))//'rd_dist.txt' ,err=1)
      open (vdf ,file=outdir(1:len_trim(outdir))//'vd_dist.txt' ,err=1)
      open (rdft,file=outdir(1:len_trim(outdir))//'rdt_dist.txt',err=1)
      open (vdft,file=outdir(1:len_trim(outdir))//'vdt_dist.txt',err=1)

      write (vlm,*)  'ix iy vol x t w n r m Td Zd Ek vr vz vt v xliq Pth
     *'
      write (core ,*) 'ix s x gc g r m Td Zd Ek vr vz vt v xliq Pth'
      write (sptr ,*) 'ix s x gc g r m Td Zd Ek vr vz vt v xliq Pth'
      write (wall ,*) 'ix s x gc g r m Td Zd Ek vr vz vt v xliq Pth'
      write (pfwl ,*) 'ix s x gc g r m Td Zd Ek vr vz vt v xliq Pth'
      write (divib,*) 'iy s x gc g r m Td Zd Ek vr vz vt v xliq Pth'
      write (divob,*) 'iy s x gc g r m Td Zd Ek vr vz vt v xliq Pth'
      if (kmesh .eq. 'dnull') then
       write (divit,*) 'iy s x gc g r m Td Zd Ek vr vz vt v xliq Pth'
       write (divot,*) 'iy s x gc g r m Td Zd Ek vr vz vt v xliq Pth'
      endif
      write (rdf  ,*) 'ix iy rd frd'
      write (vdf  ,*) 'ix iy vd fvd fvpld vtrd fvtrd'
      write (rdft ,*) 'rd frd'
      write (vdft ,*) 'vd fvd fvpld vtrd fvtrd'

      print *, 'Gd_inp',dstflux
      print *,'particle balance=',winj-
     *        (wcor+werr+wabs+depwi+depwo+
     *         depdib+depdob+depdit+depdot+
     *         deppib+deppob+deppit+deppot)
      print *, '  winj=',winj  ,'  wcor=',wcor
      print *, '  wwls=',depwi+depwo+
     *                   depdib+depdob+depdit+depdot+
     *                   deppib+deppob+deppit+deppot,
     *         '  wvlm=',wabs
      print *, '  wwli=',depwi ,'  wwlo=',depwo
      print *, ' wdvib=',depdib,' wdvob=',depdob
      print *, ' wdvit=',depdit,' wdvot=',depdot
      print *, ' wpwib=',deppib,' wpwob=',deppob
      print *, ' wpwit=',deppit,' wpwot=',deppot
      print *, '  wspr=',wspr  ,'  werr=',werr
      print *, '  fwli=',flxwi ,'  fwlo=',flxwo
      print *, ' fdvib=',flxdib,' fdvob=',flxdob
      print *, ' fdvit=',flxdit,' fdvot=',flxdot
      print *, ' fpwib=',flxpib,' fpwob=',flxpob
      print *, ' fpwit=',flxpit,' fpwot=',flxpot
      print *, '  rwli=',rflwi ,'  rwlo=',rflwo
      print *, ' rdvib=',rfldib,' rdvob=',rfldob
      print *, ' rdvit=',rfldit,' rdvot=',rfldot
      print *, ' rpwib=',rflpib,' rpwob=',rflpob
      print *, ' rpwit=',rflpit,' rpwot=',rflpot
      print *, '  piob=',pasiob,'  poib=',pasoib
      print *, '  piot=',pasiot,'  poit=',pasoit

      do 10, i=2,nx-1
       flg=0
       if ((i .gt. ixpt1(1)) .and. (i .le. ixpt2(1))) flg=1
       if (kmesh .eq. 'dnull') then
        if ((i .gt. ixpt1(2)) .and. (i .le. ixpt2(2))) flg=1
       endif
       if (flg .eq. 1) then
        if (gcore(i) .ne. 0.d0) then
         nrm=1.d0/gcore(i)
         rcore (i)=rcore (i)*nrm
         mcore (i)=mcore (i)*nrm
         tdcore(i)=tdcore(i)*nrm
         qcore (i)=qcore (i)*nrm
         ecore (i)=ecore (i)*nrm
         vrcore(i)=vrcore(i)*nrm
         vzcore(i)=vzcore(i)*nrm
         vtcore(i)=vtcore(i)*nrm
         vcore (i)=vcore (i)*nrm
         xlcore(i)=xlcore(i)*nrm
         ltcore(i)=ltcore(i)*nrm
        else
         rcore (i)=0.d0
         mcore (i)=0.d0
         tdcore(i)=0.d0
         qcore (i)=0.d0
         ecore (i)=0.d0
         vrcore(i)=0.d0
         vzcore(i)=0.d0
         vtcore(i)=0.d0
         vcore (i)=0.d0
         xlcore(i)=0.d0
         ltcore(i)=0.d0
        endif
        if (score(i) .gt. 0.d0) then
         wcore(i)=wcore(i)*dstcflx/score(i)
         gcore(i)=gcore(i)/score(i)
        else
         wcore(i)=0.d0
         gcore(i)=0.d0
        endif
        if (gsptr(i) .ne. 0.d0) then
         nrm=1.d0/gsptr(i)
         rsptr (i)=rsptr (i)*nrm
         msptr (i)=msptr (i)*nrm
         tdsptr(i)=tdsptr(i)*nrm
         qsptr (i)=qsptr (i)*nrm
         esptr (i)=esptr (i)*nrm
         vrsptr(i)=vrsptr(i)*nrm
         vzsptr(i)=vzsptr(i)*nrm
         vtsptr(i)=vtsptr(i)*nrm
         vsptr (i)=vsptr (i)*nrm
         xlsptr(i)=xlsptr(i)*nrm
         ltsptr(i)=ltsptr(i)*nrm
        else
         rsptr (i)=0.d0
         msptr (i)=0.d0
         tdsptr(i)=0.d0
         qsptr (i)=0.d0
         esptr (i)=0.d0
         vrsptr(i)=0.d0
         vzsptr(i)=0.d0
         vtsptr(i)=0.d0
         vsptr (i)=0.d0
         xlsptr(i)=0.d0
         ltsptr(i)=0.d0
        endif
        if (ssptr(i) .gt. 0.d0) then
         wsptr(i)=wsptr(i)*dstcflx/ssptr(i)
         gsptr(i)=gsptr(i)/ssptr(i)
        else
         wsptr(i)=0.d0
         gsptr(i)=0.d0
        endif
       else
        wcore (i)=0.d0
        gcore (i)=0.d0
        rcore (i)=0.d0
        mcore (i)=0.d0
        tdcore(i)=0.d0
        qcore (i)=0.d0
        ecore (i)=0.d0
        vrcore(i)=0.d0
        vzcore(i)=0.d0
        vtcore(i)=0.d0
        vcore (i)=0.d0
        xlcore(i)=0.d0
        ltcore(i)=0.d0
        wsptr (i)=0.d0
        gsptr (i)=0.d0
        rsptr (i)=0.d0
        msptr (i)=0.d0
        tdsptr(i)=0.d0
        qsptr (i)=0.d0
        esptr (i)=0.d0
        vrsptr(i)=0.d0
        vzsptr(i)=0.d0
        vtsptr(i)=0.d0
        vsptr (i)=0.d0
        xlsptr(i)=0.d0
        ltsptr(i)=0.d0
       endif

       if (gwall(i) .ne. 0.d0) then
        nrm=1.d0/gwall(i)
        rwall (i)=rwall (i)*nrm
        mwall (i)=mwall (i)*nrm
        tdwall(i)=tdwall(i)*nrm
        qwall (i)=qwall (i)*nrm
        ewall (i)=ewall (i)*nrm
        vrwall(i)=vrwall(i)*nrm
        vzwall(i)=vzwall(i)*nrm
        vtwall(i)=vtwall(i)*nrm
        vwall (i)=vwall (i)*nrm
        xlwall(i)=xlwall(i)*nrm
        ltwall(i)=ltwall(i)*nrm
       else
        rwall (i)=0.d0
        mwall (i)=0.d0
        tdwall(i)=0.d0
        qwall (i)=0.d0
        ewall (i)=0.d0
        vrwall(i)=0.d0
        vzwall(i)=0.d0
        vtwall(i)=0.d0
        vwall (i)=0.d0
        xlwall(i)=0.d0
        ltwall(i)=0.d0
       endif
       if (sbnd(i,ny) .gt. 0.d0) then
        wwall(i)=wwall(i)*dstcflx/sbnd(i,ny)
        gwall(i)=gwall(i)/sbnd(i,ny)
       else
        wwall(i)=0.d0
        gwall(i)=0.d0
       endif

       if (gpfwl(i) .ne. 0.d0) then
        nrm=1.d0/gpfwl(i)
        rpfwl (i)=rpfwl (i)*nrm
        mpfwl (i)=mpfwl (i)*nrm
        tdpfwl(i)=tdpfwl(i)*nrm
        qpfwl (i)=qpfwl (i)*nrm
        epfwl (i)=epfwl (i)*nrm
        vrpfwl(i)=vrpfwl(i)*nrm
        vzpfwl(i)=vzpfwl(i)*nrm
        vtpfwl(i)=vtpfwl(i)*nrm
        vpfwl (i)=vpfwl (i)*nrm
        xlpfwl(i)=xlpfwl(i)*nrm
        ltpfwl(i)=ltpfwl(i)*nrm
       else
        rpfwl (i)=0.d0
        mpfwl (i)=0.d0
        tdpfwl(i)=0.d0
        qpfwl (i)=0.d0
        epfwl (i)=0.d0
        vrpfwl(i)=0.d0
        vzpfwl(i)=0.d0
        vtpfwl(i)=0.d0
        vpfwl (i)=0.d0
        xlpfwl(i)=0.d0
        ltpfwl(i)=0.d0
       endif
       if (sbnd(i,1) .gt. 0.d0) then
        wpfwl(i)=wpfwl(i)*dstcflx/sbnd(i,1)
        gpfwl(i)=gpfwl(i)/sbnd(i,1)
       else
        wpfwl(i)=0.d0
        gpfwl(i)=0.d0
       endif

       write (core,100) i,score(i),xcore (i),wcore (i),gcore (i),
     *                    rcore(i),mcore (i),tdcore(i),qcore (i),
     *                    ecore(i),vrcore(i),vzcore(i),vtcore(i),
     *                    vcore(i),xlcore(i),ltcore(i)
       write (sptr,100) i,ssptr(i),xsptr (i),wsptr (i),gsptr (i),
     *                    rsptr(i),msptr (i),tdsptr(i),qsptr (i),
     *                    esptr(i),vrsptr(i),vzsptr(i),vtsptr(i),
     *                    vsptr(i),xlsptr(i),ltsptr(i)
       write (wall,100) i,sbnd(i,ny),xwall (i),wwall (i),gwall (i),
     *                               rwall (i),mwall (i),tdwall(i),
     *                               qwall (i),ewall (i),vrwall(i),
     *                               vzwall(i),vtwall(i),vwall (i),
     *                               xlwall(i),ltwall(i)
       write (pfwl,100) i,sbnd(i, 1),xpfwl (i),wpfwl (i),gpfwl (i),
     *                               rpfwl (i),mpfwl (i),tdpfwl(i),
     *                               qpfwl (i),epfwl (i),vrpfwl(i),
     *                               vzpfwl(i),vtpfwl(i),vpfwl (i),
     *                               xlpfwl(i),ltpfwl(i)
10    continue

      nstst=0.d0
      
      do 19, i=1,nx
       do 20, j=1,ny
        if (nsts(i,j) .ne. 0.d0) then
         nrm=1.d0/nsts(i,j)
         rsts (i,j)=rsts (i,j)*nrm
         msts (i,j)=msts (i,j)*nrm
         tdsts(i,j)=tdsts(i,j)*nrm
         qsts (i,j)=qsts (i,j)*nrm
         ests (i,j)=ests (i,j)*nrm
         vrsts(i,j)=vrsts(i,j)*nrm
         vzsts(i,j)=vzsts(i,j)*nrm
         vtsts(i,j)=vtsts(i,j)*nrm
         vsts (i,j)=vsts (i,j)*nrm
         xlsts(i,j)=xlsts(i,j)*nrm
         ltsts(i,j)=ltsts(i,j)*nrm
         do 21, k=0,frdim
          rdfnc(i,j,k)=rdfnc(i,j,k)*nrm
21       continue
         do 22, k=0,fvdim
          vdfnc  (i,j,k)=vdfnc  (i,j,k)*nrm
          vdfnctr(i,j,k)=vdfnctr(i,j,k)*nrm
          vdfncpl(i,j,k)=vdfncpl(i,j,k)*nrm
22       continue
        else
         rsts (i,j)=0.d0
         msts (i,j)=0.d0
         tdsts(i,j)=0.d0
         qsts (i,j)=0.d0
         ests (i,j)=0.d0
         vrsts(i,j)=0.d0
         vzsts(i,j)=0.d0
         vtsts(i,j)=0.d0
         vsts (i,j)=0.d0
         xlsts(i,j)=0.d0
         ltsts(i,j)=0.d0
         do 23, k=0,frdim
          rdfnc(i,j,k)=0.d0
23       continue
         do 24, k=0,fvdim
          vdfnc  (i,j,k)=0.d0
          vdfnctr(i,j,k)=0.d0
          vdfncpl(i,j,k)=0.d0
24       continue
        endif
        if (vol(i,j) .gt. 0.d0) then
         wsts(i,j)=wsts(i,j)*dstcflx/vol(i,j)
         nsts(i,j)=nsts(i,j)/vol(i,j)
        else
         wsts(i,j)=0.d0
         nsts(i,j)=0.d0
        endif
        write (vlm,200) i,j,vol(i,j),xsts (i,j),tsts (i,j),wsts (i,j),
     *                               nsts (i,j),rsts (i,j),msts (i,j),
     *                               tdsts(i,j),qsts (i,j),ests (i,j),
     *                               vrsts(i,j),vzsts(i,j),vtsts(i,j),
     *                               vsts (i,j),xlsts(i,j),ltsts(i,j)
        do 25, k=0,frdim
         write (rdf,400) i,j,rdnod(k),rdfnc(i,j,k)
25      continue
        do 26, k=0,fvdim
         write (vdf,400) i,j,vdnod(k),vdfnc(i,j,k),vdfncpl(i,j,k),
     *                       vtrdnod(k),vdfnctr(i,j,k)
26      continue
        nstst=nstst+nsts(i,j)
20     continue
19    continue

      do 27, k=0,frdim
       write (rdft,500) rdnod(k),rdfnct(k)/nstst
27    continue
      do 28, k=0,fvdim
       write (vdft,500) vdnod(k),vdfnct(k)/nstst,vdfncplt(k)/nstst,
     *                  vtrdnod(k),vdfnctrt(k)/nstst
28    continue
      
      do 31, j=2,ny-1
       if (gdivib(j) .ne. 0.d0) then
        nrm=1.d0/gdivib(j)
        rdivib (j)=rdivib (j)*nrm
        mdivib (j)=mdivib (j)*nrm
        tddivib(j)=tddivib(j)*nrm
        qdivib (j)=qdivib (j)*nrm
        edivib (j)=edivib (j)*nrm
        vrdivib(j)=vrdivib(j)*nrm
        vzdivib(j)=vzdivib(j)*nrm
        vtdivib(j)=vtdivib(j)*nrm
        vdivib (j)=vdivib (j)*nrm
        xldivib(j)=xldivib(j)*nrm
        ltdivib(j)=ltdivib(j)*nrm
       else
        rdivib (j)=0.d0
        mdivib (j)=0.d0
        tddivib(j)=0.d0
        qdivib (j)=0.d0
        edivib (j)=0.d0
        vrdivib(j)=0.d0
        vzdivib(j)=0.d0
        vtdivib(j)=0.d0
        vdivib (j)=0.d0
        xldivib(j)=0.d0
        ltdivib(j)=0.d0
       endif
       if (sbnd(1,j) .gt. 0.d0) then
        wdivib(j)=wdivib(j)*dstcflx/sbnd(1,j)
        gdivib(j)=gdivib(j)/sbnd(1,j)
       else
        wdivib(j)=0.d0
        gdivib(j)=0.d0
       endif

       if (gdivob(j) .ne. 0.d0) then
        nrm=1.d0/gdivob(j)
        rdivob (j)=rdivob (j)*nrm
        mdivob (j)=mdivob (j)*nrm
        tddivob(j)=tddivob(j)*nrm
        qdivob (j)=qdivob (j)*nrm
        edivob (j)=edivob (j)*nrm
        vrdivob(j)=vrdivob(j)*nrm
        vzdivob(j)=vzdivob(j)*nrm
        vtdivob(j)=vtdivob(j)*nrm
        vdivob (j)=vdivob (j)*nrm
        xldivob(j)=xldivob(j)*nrm
        ltdivob(j)=ltdivob(j)*nrm
       else
        rdivob (j)=0.d0
        mdivob (j)=0.d0
        tddivob(j)=0.d0
        qdivob (j)=0.d0
        edivob (j)=0.d0
        vrdivob(j)=0.d0
        vzdivob(j)=0.d0
        vtdivob(j)=0.d0
        vdivob (j)=0.d0
        xldivob(j)=0.d0
        ltdivob(j)=0.d0
       endif
       if (sbnd(nx,j) .gt. 0.d0) then
        wdivob(j)=wdivob(j)*dstcflx/sbnd(nx,j)
        gdivob(j)=gdivob(j)/sbnd(nx,j)
       else
        wdivob(j)=0.d0
        gdivob(j)=0.d0
       endif

       write (divib,300) j,sbnd (1,j),xdivib (j),wdivib (j),gdivib (j),
     *                    rdivib  (j),mdivib (j),tddivib(j),qdivib (j),
     *                    edivib  (j),vrdivib(j),vzdivib(j),vtdivib(j),
     *                    vdivib  (j),xldivib(j),ltdivib(j)
       write (divob,300) j,sbnd(nx,j),xdivob (j),wdivob (j),gdivob (j),
     *                    rdivob  (j),mdivob (j),tddivob(j),qdivob (j),
     *                    edivob  (j),vrdivob(j),vzdivob(j),vtdivob(j),
     *                    vdivob  (j),xldivob(j),ltdivob(j)
31    continue

      if (kmesh .eq. 'dnull') then
       do 32, j=2,ny-1
        if (gdivit(j) .ne. 0.d0) then
         nrm=1.d0/gdivit(j)
         rdivit (j)=rdivit (j)*nrm
         mdivit (j)=mdivit (j)*nrm
         tddivit(j)=tddivit(j)*nrm
         qdivit (j)=qdivit (j)*nrm
         edivit (j)=edivit (j)*nrm
         vrdivit(j)=vrdivit(j)*nrm
         vzdivit(j)=vzdivit(j)*nrm
         vtdivit(j)=vtdivit(j)*nrm
         vdivit (j)=vdivit (j)*nrm
         xldivit(j)=xldivit(j)*nrm
         ltdivit(j)=ltdivit(j)*nrm
        else
         rdivit (j)=0.d0
         mdivit (j)=0.d0
         tddivit(j)=0.d0
         qdivit (j)=0.d0
         edivit (j)=0.d0
         vrdivit(j)=0.d0
         vzdivit(j)=0.d0
         vtdivit(j)=0.d0
         vdivit (j)=0.d0
         xldivit(j)=0.d0
         ltdivit(j)=0.d0
        endif
        if (sbnd(ixrb1,j) .gt. 0.d0) then
         wdivit(j)=wdivit(j)*dstcflx/sbnd(ixrb1,j)
         gdivit(j)=gdivit(j)/sbnd(ixrb1,j)
        else
         wdivit(j)=0.d0
         gdivit(j)=0.d0
        endif

        if (gdivot(j).ne.0.d0) then
         nrm=1.d0/gdivot(j)
         rdivot (j)=rdivot (j)*nrm
         mdivot (j)=mdivot (j)*nrm
         tddivot(j)=tddivot(j)*nrm
         qdivot (j)=qdivot (j)*nrm
         edivot (j)=edivot (j)*nrm
         vrdivot(j)=vrdivot(j)*nrm
         vzdivit(j)=vzdivit(j)*nrm
         vtdivot(j)=vtdivot(j)*nrm
         vdivot (j)=vdivot (j)*nrm
         xldivot(j)=xldivot(j)*nrm
         ltdivot(j)=ltdivot(j)*nrm
        else
         rdivot (j)=0.d0
         mdivot (j)=0.d0
         tddivot(j)=0.d0
         qdivot (j)=0.d0
         edivot (j)=0.d0
         vrdivot(j)=0.d0
         vzdivot(j)=0.d0
         vtdivot(j)=0.d0
         vdivot (j)=0.d0
         xldivot(j)=0.d0
         ltdivot(j)=0.d0
        endif
        if (sbnd(ixlb2+1,j) .gt. 0.d0) then
         wdivot(j)=wdivot(j)*dstcflx/sbnd(ixlb2+1,j)
         gdivot(j)=gdivot(j)/sbnd(ixlb2+1,j)
        else
         wdivot(j)=0.d0
         gdivot(j)=0.d0
        endif

        write (divit,300) j,sbnd (ixrb1,j),
     *                     xdivit  (j),wdivit (j),gdivit (j),
     *                     rdivit  (j),mdivit (j),tddivit(j),qdivit (j),
     *                     edivit  (j),vrdivit(j),vzdivit(j),vtdivit(j),
     *                     vdivit  (j),xldivit(j),ltdivit(j)
        write (divot,300) j,sbnd(ixlb2+1,j),
     *                     xdivot  (j),wdivot (j),gdivot (j),
     *                     rdivot  (j),mdivot (j),tddivot(j),qdivot (j),
     *                     edivot  (j),vrdivot(j),vzdivot(j),vtdivot(j),
     *                     vdivot  (j),xldivot(j),ltdivot(j)
32     continue
      endif

100   format(I3,15(e25.15e3))
200   format(2(I3),16(e25.15e3))
300   format(I3,15(e25.15e3))
400   format(2(I3),5(e25.15e3))
500   format(5(e25.15e3))

      close (vlm  )
      close (core )
      close (sptr )
      close (wall )
      close (pfwl )
      close (divib)
      close (divob)
      if (kmesh .eq. 'dnull') then
       close (divit)
       close (divot)
      endif
      close (rdf )
      close (vdf )
      close (rdft)
      close (vdft)
      
      goto 2
1     print *, '***Error: cannot save stat data!'
2     return
      end
