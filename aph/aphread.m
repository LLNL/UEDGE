c-----------------------------------------------------------------------
      subroutine aphread

c ... Set up tables for hydrogenic atomic-physics processes.

      implicit none

c ... Common blocks:
      Use(Dim)
      Use(Share)        # nhdf, hdfilename
      Use(Data_input)   # istabon, aphdir, data_directory
      Use(Rtdegas)      # mpd, mpe

c ... Function:
      integer utgetcl   # defined in the Basis runtime library

c ... Local variable:
      integer MAXSTRING
      parameter (MAXSTRING=500)
      character*(MAXSTRING) aphdirx
      character*(MAXSTRING) adname
      character*(MAXSTRING) dataDir

      dataDir=data_directory

c ... Get length of string containing directory name (omitting trailing
c     blanks). Basis function basfilex expands $,~.
c     NOTE: basfilex will not exist for Python; 
c     FACETS replaces aphdirx with aphdir in calls to findFile below
      call basfilex(aphdir,aphdirx)

c...  If istabon>0, set-up tables for ionization, recombination,
c...  charge exchange and energy loss.
      if (istabon .eq. 1) then
        call findFile('rates.adpak',  aphdirx, dataDir, adname, isaphdir)
        call readrt(TRIM(adname))
      elseif (istabon .eq. 2) then
        call findFile('rates.strahl', aphdirx, dataDir, adname, isaphdir)
        call readrt(TRIM(adname))
      elseif (istabon .eq. 3) then
         mpe=48
         mpd=11
         mpr=1
         call gallot("Rtdegas",0)
         call findFile('eh.dat',   aphdirx, dataDir, adname, isaphdir)
         call readeh(TRIM(adname))
         call findFile('atmc.dat', aphdirx, dataDir, adname, isaphdir)
         call readatmc(TRIM(adname))
         call setauxvar
      elseif (istabon .eq. 4) then
         mpe=60
         mpd=15
         mpr=1
         call gallot("Rtdegas",0)
         call findFile('nwfits', aphdirx, dataDir, adname,isaphdir)
         call readnw(TRIM(adname))
         call setauxvar
      elseif ( (istabon .eq. 5) .or. (istabon .eq. 6) ) then
         mpe=60
         mpd=15
         mpr=1
         call gallot("Rtdegas",0)
         call findFile('nwfits',   aphdirx, dataDir, adname, isaphdir)
         call readnw(TRIM(adname))
         call setauxvar
         call splined
      elseif (istabon .eq. 8) then
         mpe=60
         mpd=15
         mpr=1
         call gallot("Rtdegas",0)
         call findFile('ehr1.dat', aphdirx, dataDir, adname, isaphdir)
         call readehr1(TRIM(adname))
         call setauxvar
      elseif (istabon .eq. 9 .or. istabon .eq. 10) then
         mpe=60
         mpd=15
         mpr=1
         call gallot("Rtdegas",0)
         call findFile('ehr2.dat', aphdirx, dataDir, adname, isaphdir)
         call readehr1(TRIM(adname))
         call setauxvar
      elseif (istabon .eq. 11) then
         mpe=60
         mpd=15
         mpr=1
         call gallot("Rtdegas",0)
         call findFile('thin.dat', aphdirx, dataDir, adname, isaphdir)
         call readehr2(TRIM(adname))
         call setauxvar
      elseif (istabon .eq. 12) then
         mpe=60
         mpd=15
         mpr=1
         call gallot("Rtdegas",0)
         call findFile('thickLyA.dat', aphdirx, dataDir, adname,isaphdir)
         call readehr2(TRIM(adname))
         call setauxvar
      elseif (istabon .eq. 13) then
         mpe=60
         mpd=15
         mpr=1
         call gallot("Rtdegas",0)
         call findFile('thickAllLy.dat',aphdirx,dataDir,adname,isaphdir)
         call readehr2(TRIM(adname))
         call setauxvar
      elseif (istabon .eq. 14 .or. istabon .eq. 15) then
         mpe=60
         mpd=15
         mpr=30
         call gallot("Rtdegas",0)
         call findFile('ehrtau.dat',aphdirx,dataDir,adname,isaphdir)
         call readehr1(TRIM(adname))
         call setauxvar
      elseif (istabon .eq. 16) then
         call readmc(nhdf,hdfilename)
c ...       Note that the rate parameter tables (in group Multicharge) have
c ...       the same ne-dimension and te-dimension for all species, so
c ...       the hydrogenic tables must have the same dimensions as the 
c ...       impurity tables, if both are present.
      endif

      return
      end
c ... ............................................................
c ...  Determine where the data files exist.  The search order is:
c ...   dir1, then dir2, then $PWD
c ...  If flag is 1, then just assume $PWD and return
c ... ............................................................
       subroutine findFile(basename, dir1, dir2, fullname,flag)
c
       implicit none
       character(*), intent(in) ::  basename, dir1, dir2
       integer, intent(in) ::  flag
       character(*), intent(out) :: fullname
       character*(500) :: fname
       logical fileExists

cc      if (flag/=1) then 
       if (flag .ne. 1) then
         fullname=basename
         return
       endif
       fname=TRIM(dir1) // '/'//TRIM(basename)
       INQUIRE(FILE=TRIM(fname),EXIST=fileExists)
       if (fileExists) then
         fullname=fname
       else
         fname=TRIM(dir2) // '/'//TRIM(basename)
         INQUIRE(FILE=TRIM(fname),EXIST=fileExists)
         if (fileExists) then
           fullname=fname
         else
           fname=basename
           INQUIRE(FILE=TRIM(fname),EXIST=fileExists)
           if (fileExists) then
             fullname=fname
           else
             write(*,*) "Cannot find "//TRIM(fname)//" in:"
             write(*,*) TRIM(dir1)
             write(*,*) TRIM(dir2)
             write(*,*) " or current directory"
           endif
         endif
       endif
       return
cc       end subroutine
      end
c-----------------------------------------------------------------------
      subroutine readrt (fname)
      implicit none
      character*(*) fname
Use(Rtdata)

c     local variables --
      integer ios, n, nget
      character idcod*8, idtyp*8, id1*32

c     procedures --
      external freeus, kaboom, remark, gallot, readrt1

c----------------------------------------------------------------------c
c     Read ADPAK rate data from 'un*formatted' file (B. Braams)
c----------------------------------------------------------------------c

      call freeus(nget)
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call remark('**** ADPAK data file not found; set aphdir path')
         call kaboom(0)
      endif

c     read header --
*     un*formatted read for header data
      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,'(1x,1a120)') labelht

c     read dimensions --
*     un*formatted read for integer data
      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) htnt,htnn,htns

c     allocate storage --
      call gallot("Rtdata",0)

c     read abscissae and rates --
      call readrt1 (nget)

      close (nget)

      return
      end
c-----------------------------------------------------------------------
      subroutine readrt1 (nget)
      implicit none
      integer nget
Use(Rtdata)

c     local variables --
      integer i, j, k, n
      character idcod*8, idtyp*8, id1*32

c     read abscissae --
*     un*formatted read for real data

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (htza(i),i=0,htns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (htzn(i),i=0,htns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (htt(i),i=0,htnt)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (htn(i),i=0,htnn)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (htlt(i),i=0,htnt)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (htln(i),i=0,htnn)

c     read rate coefficients --
*     un*formatted read for real data

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((htlsa(i,j,k),i=0,htnt),j=0,htnn),k=0,htns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((htlra(i,j,k),i=0,htnt),j=0,htnn),k=0,htns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((htlqa(i,j,k),i=0,htnt),j=0,htnn),k=0,htns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((htlcx(i,j,k),i=0,htnt),j=0,htnn),k=0,htns-1)

      return
      end
c-----------------------------------------------------------------------
      subroutine readnw (fname)
      implicit none
      character*(*) fname
Use(Rtdegas)

c     local variables --
      integer ios, nget, jd, jt, jr

c     procedures --
      external freeus, kaboom, remark

c----------------------------------------------------------------------c
c     Read density-dependent rate data from POST93 data file 'nwfits'
c----------------------------------------------------------------------c

      call freeus(nget)
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call remark('**** data file nwfits not found; set aphdir path')
         call kaboom(0)
      endif

      jr = 1

c     ionization rate (cm**3/sec):
         read(nget,9012)((wsveh(jt,jd,jr),jt=1,mpe),jd=1,mpd)
c     recombination rate (cm**3/sec):
         read(nget,9012)((wsveh0(jt,jd,jr),jt=1,mpe),jd=1,mpd)
c     hydrogen line radiation (erg-cm**3/sec):
         read(nget,9012)((wlemiss(jt,jd),jt=1,mpe),jd=1,mpd)

      close (nget)

c     convert to SI units:
      do jt=1,mpe
         do jd=1,mpd
            wsveh(jt,jd,jr)=max(1.e-50, wsveh(jt,jd,jr))*1.0e-06
            wsveh0(jt,jd,jr)=max(1.e-50, wsveh0(jt,jd,jr))*1.0e-06
            wlemiss(jt,jd)=max(1.e-75, wlemiss(jt,jd))*1.0e-13
         enddo
      enddo

 9012 format(10(6(1x,e12.5)/))

      return
      end
c-----------------------------------------------------------------------
      subroutine readeh (fname)
      implicit none
      character*(*) fname
Use(Rtdegas)

c     local variables --
      integer ios, nget, jd, jt, jr

c     procedures --
      external freeus, kaboom, remark

c----------------------------------------------------------------------c
c     Read density-dependent rate data from DEGAS file 'eh.dat'
c----------------------------------------------------------------------c

      call freeus(nget)
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call remark('**** DEGAS file eh.dat not found; set aphdir path')
         call kaboom(0)
      endif

      jr = 1

c     ionization rate (cm**3/sec):
         read(nget,9012)((wsveh(jt,jd,jr),jt=1,mpe),jd=1,mpd)
c     recombination rate (cm**3/sec):
         read(nget,9012)((wsveh0(jt,jd,jr),jt=1,mpe),jd=1,mpd)
c     hydrogen line radiation (erg-cm**3/sec):
         read(nget,9012)((wlemiss(jt,jd),jt=1,mpe),jd=1,mpd)
c     electron energy loss (eV) per ionization:
         read(nget,9012)((welms(jt,jd),jt=1,mpe),jd=1,mpd)
c     n=3 excited state fraction:
         read(nget,9012)((pne3(jt,jd),jt=1,mpe),jd=1,mpd)
c     n=2 excited state fraction:
         read(nget,9012)((pne2(jt,jd),jt=1,mpe),jd=1,mpd)

      close (nget)

c     convert to SI units:
      do jt=1,mpe
         do jd=1,mpd
            wsveh(jt,jd,jr)=max(1.e-50, wsveh(jt,jd,jr))*1.0e-06
            wsveh0(jt,jd,jr)=max(1.e-50, wsveh0(jt,jd,jr))*1.0e-06
            wlemiss(jt,jd)=max(1.e-75, wlemiss(jt,jd))*1.0e-13
         enddo
      enddo

 9012 format(8(6(1x,e12.5)/))

      return
      end
c-----------------------------------------------------------------------
      subroutine readehr1 (fname)
      implicit none
      character*(*) fname
Use(Rtdegas)

c     local variables --
      integer ios, nget, jd, jt, jr
      character*80 zdummy

c     procedures --
      external freeus, kaboom, remark

c----------------------------------------------------------------------c
c     Read density-dependent hydrogenic rate file, e.g., 'ehr1.dat'
c----------------------------------------------------------------------c

      call freeus(nget)
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call remark('**** hydrogenic rate file not found; set aphdir 
     . path or isaphdir = 0')
         call kaboom(0)
      endif

      do 100 jr=1,mpr

      if (mpr.gt.1) then
         read(nget,9013) zdummy  # reads extra header lines in ehrtau.dat
      endif

c     ionization rate (cm**3/sec):
      read(nget,9013) zdummy
      do 12 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(wsveh(jt,jd,jr),jt=1,mpe)
 12   continue
c     recombination rate (cm**3/sec):
      read(nget,9013) zdummy
      do 14 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(wsveh0(jt,jd,jr),jt=1,mpe)
 14   continue
c     neutral-electron radiation loss rate (erg/sec):
      read(nget,9013) zdummy
      do 16 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(welms1(jt,jd,jr),jt=1,mpe)
 16   continue
c     continuum-electron radiation loss rate (erg/sec):
      read(nget,9013) zdummy
      do 18 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(welms2(jt,jd,jr),jt=1,mpe)
 18   continue

c     neutral "n=3/n=1" ratio:
      read(nget,9013) zdummy
      do 20 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne31(jt,jd),jt=1,mpe)
 20   continue
c     continuum "n=3/n=1" ratio:
      read(nget,9013) zdummy
      do 22 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne32(jt,jd),jt=1,mpe)
 22   continue
c     neutral "n=2/n=1" ratio:
      read(nget,9013) zdummy
      do 24 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne21(jt,jd),jt=1,mpe)
 24   continue
c     continuum "n=2/n=1" ratio:
      read(nget,9013) zdummy
      do 26 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne22(jt,jd),jt=1,mpe)
 26   continue

100   continue

      close (nget)

c     convert to SI units:
      do jt=1,mpe
         do jd=1,mpd
            do jr=1,mpr
               wsveh(jt,jd,jr)=max(1.e-50, wsveh(jt,jd,jr))*1.0e-06
               wsveh0(jt,jd,jr)=max(1.e-50, wsveh0(jt,jd,jr))*1.0e-06
               welms1(jt,jd,jr)=max(1.e-50, welms1(jt,jd,jr))*1.0e-07
               welms2(jt,jd,jr)=max(1.e-50, welms2(jt,jd,jr))*1.0e-07
            enddo
         enddo
      enddo

 9012 format(10(6(1x,e12.5)/))
 9013 format(a80)

      return
      end
c-----------------------------------------------------------------------
      subroutine readehr2 (fname)
      implicit none
      character*(*) fname
Use(Rtdegas)

c     local variables --
      integer ios, nget, jd, jt
      character*80 zdummy

c     procedures --
      external freeus, kaboom, remark

c----------------------------------------------------------------------c
c     Read dens-depend hydro. + n=2-9 data, e.g., thin.dat
c----------------------------------------------------------------------c

      call freeus(nget)
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call remark('**** hydrogenic rate file not found; set aphdir
     . path')
         call kaboom(0)
      endif

c     ionization rate (cm**3/sec):
      read(nget,9013) zdummy
      do 12 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(wsveh(jt,jd,1),jt=1,mpe)
 12   continue
c     recombination rate (cm**3/sec):
      read(nget,9013) zdummy
      do 14 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(wsveh0(jt,jd,1),jt=1,mpe)
 14   continue
c     neutral-electron radiation loss rate (erg/sec):
      read(nget,9013) zdummy
      do 16 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(welms1(jt,jd,1),jt=1,mpe)
 16   continue
c     continuum-electron radiation loss rate (erg/sec):
      read(nget,9013) zdummy
      do 18 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(welms2(jt,jd,1),jt=1,mpe)
 18   continue
c     neutral "n=3/n=1" ratio:
      read(nget,9013) zdummy
      do 20 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne31(jt,jd),jt=1,mpe)
 20   continue
c     continuum "n=3/n=1" ratio:
      read(nget,9013) zdummy
      do 22 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne32(jt,jd),jt=1,mpe)
 22   continue
c     neutral "n=2/n=1" ratio:
      read(nget,9013) zdummy
      do 24 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne21(jt,jd),jt=1,mpe)
 24   continue
c     continuum "n=2/n=1" ratio:
      read(nget,9013) zdummy
      do 26 jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne22(jt,jd),jt=1,mpe)
 26   continue
c     neutral "n=4/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne41(jt,jd),jt=1,mpe)
      enddo
c     continuum "n=4/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne42(jt,jd),jt=1,mpe)
      enddo
c     neutral "n=5/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne51(jt,jd),jt=1,mpe)
      enddo
c     continuum "n=5/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne52(jt,jd),jt=1,mpe)
      enddo
c     neutral "n=6/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne61(jt,jd),jt=1,mpe)
      enddo
c     continuum "n=6/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne62(jt,jd),jt=1,mpe)
      enddo
c     neutral "n=7/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne71(jt,jd),jt=1,mpe)
      enddo
c     continuum "n=7/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne72(jt,jd),jt=1,mpe)
      enddo
c     neutral "n=8/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne81(jt,jd),jt=1,mpe)
      enddo
c     continuum "n=8/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne82(jt,jd),jt=1,mpe)
      enddo
c     neutral "n=9/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne91(jt,jd),jt=1,mpe)
      enddo
c     continuum "n=9/n=1" ratio:
      read(nget,9013) zdummy
      do jd=1,mpd
         read(nget,9013) zdummy
         read(nget,9012)(pne92(jt,jd),jt=1,mpe)
      enddo

      close (nget)

c     convert to SI units:
      do jt=1,mpe
         do jd=1,mpd
            wsveh(jt,jd,1)=max(1.e-50, wsveh(jt,jd,1))*1.0e-06
            wsveh0(jt,jd,1)=max(1.e-50, wsveh0(jt,jd,1))*1.0e-06
            welms1(jt,jd,1)=max(1.e-50, welms1(jt,jd,1))*1.0e-07
            welms2(jt,jd,1)=max(1.e-50, welms2(jt,jd,1))*1.0e-07
         enddo
      enddo

 9012 format(10(6(1x,e12.5)/))
 9013 format(a80)

      return
      end
c-----------------------------------------------------------------------
      subroutine readatmc (fname)
      implicit none
      character*(*) fname
Use(Rtdegas)

c     local variables --
      integer ios, nget, je, jt

c     procedures --
      external freeus, kaboom, remark

c----------------------------------------------------------------------c
c     Read charge exchange rate data from DEGAS file 'atmc.dat'
c----------------------------------------------------------------------c

      call freeus(nget)
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call remark('**** DEGAS file atmc.dat not found; set aphdir
     . path')
         call kaboom(0)
      endif

c     unused data:
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)

      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)

      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)
      read(nget,9015)((svdum2(jt,je),jt=1,mpe),je=1,mpe)

c     charge exchange rate (cm**3/sec):
         read(nget,9015)((svphcx(jt,je),jt=1,mpe),je=1,mpe)

      close (nget)

c     convert to SI units:
      do jt=1,mpe
         do je=1,mpe
            svphcx(jt,je)=svphcx(jt,je)*1.0e-06
         enddo
      enddo

 9015 format(/,8(6(1x,e12.5)/))

      return
      end
c-----------------------------------------------------------------------
      subroutine setauxvar
      implicit none
Use(Dim)
Use(Share)        # istabon
Use(Data_input)
Use(Rtdegas)

c-----------------------------------------------------------------------
c     Set density and temperature data for DEGAS and POST93 rate tables
c-----------------------------------------------------------------------

c     local variables --
      integer jd,je
      real ddkpt,dekpt

c     dkpt = log10 of density(/m**3) :
      ddkpt = 0.5
      dkpt(1)=16.0
      do jd=2,mpd
         dkpt(jd)=dkpt(jd-1)+ddkpt
      enddo
      rldmin=dkpt(1)
      rldmax=dkpt(mpd)
      drefmin=10.0**rldmin
      drefmax=10.0**rldmax
      deldkpt=(rldmax-rldmin)/real(mpd-1)

c     ekpt = natural log of temperature(eV) :
      if (istabon .eq. 3) then
         ekpt(1)=0.0		# old DEGAS tables start at 1 eV
      else
         ekpt(1)=-1.2*log(10.0)	# new tables start at .06 eV
      endif

      dekpt = 0.1
      do je=2,mpe
         ekpt(je)=ekpt(je-1)+dekpt*log(10.0)
      enddo
      rlemin=ekpt(1)
      rlemax=ekpt(mpe)
      erefmin=exp(rlemin)
      erefmax=exp(rlemax)
      delekpt=(rlemax-rlemin)/real(mpe-1)

      taumin = 1.e-4             # minimum tabulated tau
      taumax = 1.e+10            # maximum tabulated tau
      deltau = log10(taumax/taumin)/(real(mpr-2))

      return
      end
c-----------------------------------------------------------------------

      subroutine splined
      implicit none
Use(Rtdegas)
Use(Aphwrk)
      external gchange, splined1

c     Construct 2-dimensional B-spline representation for atomic
c     physics rates as a function of temperature and density
c     (data from POST 93 tables)

c     Allocate arrays for spline fitting --
      nxdata=mpe		# temperature
      nydata=mpd		# density
      call gchange("Aphwrk",0)

      call splined1

      return
      end

c----------------------------------------------------------------------c

      subroutine splined1
      implicit none
Use(Dim)
Use(Share)        # istabon
Use(Data_input)
Use(Rtdegas)
Use(Aphwrk)
      integer i,j
 
c     Define data arrays --
      do i=1,nxdata
        xdata(i)=ekpt(i)
      enddo
      do j=1,nydata
        ydata(j)=dkpt(j)
      enddo
      ldf=nxdata

c     Define the order of the spline fit
c     kxords=4		# cubic in x=log(temperature)
c     kyords=4		# cubic in y=log10(density)

c     Compute the coefficients --
c     first, for ionization:
      do i=1,nxdata
        do j=1,nydata
           if (istabon .eq. 5) then
             fdata(i,j)=log10( wsveh(i,j,1) )
           elseif (istabon .eq. 6) then
             fdata(i,j)=wsveh(i,j,1)
           endif
        enddo
      enddo
      iflag = 1
      call s2copy (nxdata,nydata,fdata,1,nxdata,rsacoef,1,nxdata)
      call B2INhT (xdata, nxdata, ydata, nydata, kxords, kyords,
     .            xknots, yknots, rsacoef, ldf, workh, iflag)
c     next, for recombination:
      do i=1,nxdata
        do j=1,nydata
           if (istabon .eq. 5) then
             fdata(i,j)=log10( wsveh0(i,j,1) )
           elseif (istabon .eq. 6) then
             fdata(i,j)=wsveh0(i,j,1)
           endif
        enddo
      enddo
      iflag = 1
      call s2copy (nxdata,nydata,fdata,1,nxdata,rracoef,1,nxdata)
      call B2INhT (xdata, nxdata, ydata, nydata, kxords, kyords,
     .            xknots, yknots, rracoef, ldf, workh, iflag)
c     next, for hydrogen line radiation:
      do i=1,nxdata
        do j=1,nydata
           if (istabon .eq. 5) then
             fdata(i,j)=log10( wlemiss(i,j) )
           elseif (istabon .eq. 6) then
             fdata(i,j)=wlemiss(i,j)
           endif
        enddo
      enddo
      iflag = 1
      call s2copy (nxdata,nydata,fdata,1,nxdata,rqacoef,1,nxdata)
      call B2INhT (xdata, nxdata, ydata, nydata, kxords, kyords,
     .            xknots, yknots, rqacoef, ldf, workh, iflag)
      return
      end

c ------------- Routine to get set aphdir

      subroutine uedge_setDataDirectory(passedDataDirname)
      implicit none
      Use(Data_input)   # istabon, aphdir, data_directory
      character*(*), intent(in) :: passedDataDirname
c     Allows the framework to set the data directory
         data_directory=passedDataDirname
      return
      end
