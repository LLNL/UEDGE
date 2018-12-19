      subroutine psimatinp(ierr)
      implicit none
      
      integer ierr

      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      integer nm(Nmat_ptm),np(Nprj_ptm)
      real*8  ap(Nprj_ptm)

      character*24 a24
      character*1 aa(24)
      character*(*) outdir
      equivalence (a24,aa(1))

      character*1 a

      integer NCHmax,ican
      data NCHmax /92/, ican /11/

      integer n,i,k

      parameter (outdir='./out/')
      
      ierr=0
      open(ican,file='./dat/psitab.inp', err=1)

      read(ican,100) a

      n=0
      read(ican,100) a
100   format(A1)
      read(ican,101) n
101   format(i3)
      if(n .lt. 1) then
       ierr=2
       goto 3
      endif
      if(n .gt. Nmat_ptm) then
       ierr=2
       goto 3
      endif
      Nmat_pt=n

      read(ican,100) a

      do 66 k=1,Nmat_pt
       a24='_'
       read(ican,99) a24
99     format(A24)
       i=len(a24)
       if (i .lt. 1) then
        ierr=99
        goto 3
       endif
       Namfpsi(k)=outdir//a24
66    continue

      n=0
      read(ican,100) a
      read(ican,101) n
      if(n .lt. 1) then
      ierr=3
      go to 3
      endif
      if(n .gt. Nprj_ptm) then
      ierr=3
      go to 3
      endif
      Nprj_pt=n

      do 12 i=1,Nmat_ptm
       nm(i)=0
12    continue
      read(ican,100) a
      read(ican,102) (nm(i),i=1,Nmat_pt)
102   format(10I3)
      do 11 i=1,Nmat_pt
       if(nm(i) .lt. 1) then
       ierr=4
       go to 3
       endif
       if(nm(i) .gt. NCHmax) then
       ierr=4
       go to 3
       endif
       IZmatter(i)=nm(i)
11    continue

      do 13 i=1,Nprj_ptm
      np(i)=0
13    continue
      read(ican,100) a
      read(ican,102) (np(i),i=1,Nprj_pt)
      do 14 i=1,Nprj_pt
       if(np(i) .lt. 1) then
       ierr=5
       go to 3
       endif
       if(np(i) .gt. NCHmax) then
       ierr=5
       go to 3
       endif
       IZprjctl(i)=np(i)
14    continue

      do 15 i=1,Nprj_ptm
      ap(i)=0.0
15    continue
      read(ican,100) a
      read(ican,103) (ap(i),i=1,Nprj_pt)
103   format(10(F8.2))
      do 16 i=1,Nprj_pt
       if(ap(i) .lt. 0.5) then
       ierr=6
       go to 3
       endif
       if(ap(i) .gt. 1000.0) then
       ierr=6
       go to 3
       endif
       AZprjctl(i)=ap(i)
16    continue

3     close (ican)
      go to 2 

1     ierr=1
2     if(ierr .gt. 0) then
      print *,'error psi input=',ierr
      endif

      end


      subroutine psitabicro(ierr)
      implicit none
      
      integer ierr

      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      integer i,j,k

      ierr=0
      k=Nmat_pt*Nprj_pt
      if(k .gt. Npsi_ptm) then
       print *, 'error in input'
       print *, 'Nmat=',Nmat_pt,' Nprj=',Nprj_pt,' k=',k,
     *          'Npsi=',Nprj_ptm
       ierr=1
       goto 1
      endif

      Npsi_pt=k

      do 11 i=1,Npsi_ptm
       I_wal(i)=0
       I_prj(i)=0
11    continue

      do 12 i=1,Nmat_ptm
       do 13 j=1,Nprj_ptm
        I_pwi(j,i)=0
13     continue
12    continue

      k=0
      do 10 i=1,Nmat_pt
       do 20 j=1,Nprj_pt
        k=k+1
        I_prj(k)=j
        I_wal(k)=i
        I_pwi(j,i)=k
20     continue
10    continue

1     continue

      end
