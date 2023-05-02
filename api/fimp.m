      subroutine imprates(temp,kk,nzarg,rioniz,rrecomb,rcxrecom)
      real temp
      integer kk, nzarg
      real rioniz, rrecomb, rcxrecom
c ... Given temperature "temp" and a charge state "k", which is less
c     than or equal to the highest state "nzarg", interpolate from
c     tabulated rates for a particular impurity of ionization,
c     recombination, and charge-exchange recombination.
c    Note:  no scaling of temperatures is done here, so temp and tevb
c    must be provided in the same set of units.  tevb may not be in eV.
c
      Use(Multicharge)       # tevb,rsi,rre,rrcx

      integer itemp
      real xltemn, dlogte
c
      rrecomb = 0.
      rcxrecom = 0.
c
      xltemn = log10(tevb(1))
      dlogte = log10(tevb(2)) - xltemn
c
c ... Find index itemp into table such that
c        tevb(itemp) .le. temp .lt. tevb(itemp+1)
c
      itemp = int( ( log10( temp ) - xltemn ) / dlogte + 1. )
c
c ... For temperatures below minimum table temperature, extrapolate
c     downwards from table entries 1 and 2.

      itemp = max(1, itemp)
c
c ... For temperatures above maximum table temperature, extrapolate
c     upwards from table entries ntev-1 and ntev.

      itemp = min(ntev-1, itemp)
c
       if(kk .lt. nzarg)then
          rioniz = rsi(itemp,kk) + ( temp - tevb(itemp) )
     .   * ( rsi(itemp+1,kk) - rsi(itemp,kk) )
     .   / ( tevb(itemp+1) - tevb(itemp) )
          if(kk .eq. 0) return
       else
          rioniz = 0.0
       endif
c
          rrecomb = rre(itemp,kk) + ( temp - tevb(itemp) )
     .   * ( rre(itemp+1,kk) - rre(itemp,kk) )
     .   / ( tevb(itemp+1) - tevb(itemp) )
c
          rcxrecom = rrcx(itemp,kk) + ( temp - tevb(itemp) )
     .   * ( rrcx(itemp+1,kk) - rrcx(itemp,kk) )
     .   / ( tevb(itemp+1) - tevb(itemp) )

      return
      end
c-----End of subroutine imprates----------------------------------------

      subroutine mcrates(ne,tmpe,tmpi,za,zamax,zn,rion,rrec,rcxr)
      implicit none
      real ne,tmpe,tmpi
      integer za,zamax,zn
      real rion, rrec, rcxr
c
c ... Inputs are:
c        electron density, ne;
c        electron temperature, tmpe;
c        neutral hydrogen temperature, tmpi, for charge-exchange;
c        atomic charge state, za;
c        maximum atomic charge state, zamax;
c        nuclear charge state, zn;

c ... Outputs are:
c        rate parameters (sigma*v) for ionization, recombination,
c        and charge-exchange recombination on neutral hydrogen
c
c     The tables used in this subroutine are generated with a code
c     supplied by Bas Braams.  The data file it produces is called
c     'b2frates' by default.  This gives rates that may depend on both
c     density and temperature.
c
c     Input electron density is given in [parts/m**3]
c     Input electron temperature is given in [J].
c     Input temperature for neutral hydrogen is given in [J/AMU].
c     Table densities are given in [parts/m**3]
c     Table temperatures are given in [eV].
c     Output and table rates are all given in [m**3/s].
c
      Use(Share)               # cutlo
      Use(Physical_constants2) # ev2
      Use(Multicharge)         # rtln,rtlt,rtlsa,rtlra,rtlcx,iscxfit,isrtndep
c
      integer i1e,i1i,i,ii,j1
      real nenonz,y,dlogn,fy
      real tmpenonz,tmpinonz,xte,xti,dlogt,fxte,fxti,lrion,lrrec
c
      real rcxr_zn6, rcxr_zn6b
      external rcxr_zn6, rcxr_zn6b 
c
      rion = 0.
      rrec = 0.
      rcxr = 0.
c
      tmpenonz = max(tmpe, cutlo)  # to avoid possible log(0) error
      tmpinonz = max(tmpi, cutlo)  # to avoid possible log(0) error
      xte = log(tmpenonz/ev2)
      xti = log(tmpinonz/ev2)
      dlogt = rtlt(1) - rtlt(0)
c
c ... Find index i1 in temperature table such that 
c                 rtlt(i1) .le. xt .lt. rtlt(i1+1)
c     or, equivalently,
c                  rtt(i1) .le. tmp .lt. rtt(i1+1).
c
      i1e = int( (xte-rtlt(0))/dlogt )
      i1i = int( (xti-rtlt(0))/dlogt )
c
c ... For temperatures below minimum table temperature, extrapolate
c     downwards from table entries 0 and 1.
c
      i1e = max(0, i1e)
      i1i = max(0, i1i)
c
c ... For temperatures above maximum table temperature, extrapolate
c     upwards from table entries rtnt-1 and rtnt.
c
      i1e = min(rtnt-1, i1e)
      i1i = min(rtnt-1, i1i)
c
c ... Compute coefficient for linear interpolation.
c
      fxte = (xte-rtlt(i1e))/(rtlt(i1e+1)-rtlt(i1e))
      fxti = (xti-rtlt(i1i))/(rtlt(i1i+1)-rtlt(i1i))
c
c ... Now the density dependence --
c     Default is to use the lowest density in the table -
      j1 = 0
      fy = 0.
c     Otherwise, use linear interpolation -
      if (isrtndep .ne. 0) then
         nenonz = max(ne, cutlo)
         y = log(nenonz)
         dlogn = rtln(1) - rtln(0)
         j1 = int( (y-rtln(0))/dlogn )
         j1 = max(0, j1)
         j1 = min(rtnn-1, j1)
         fy = (y-rtln(j1))/(rtln(j1+1)-rtln(j1))
c        Do not extrapolate beyond table minimum and maximum densities
         fy = max(0., fy)
         fy = min(1., fy)
      endif
c
c ... For given za and zn, find the species index, ii, in the table.
c
      ii = -1
      do i=0,rtnsd-1
         if ((zn .eq. nint(rtzn(i))) .and. (za .eq. nint(rtza(i)))) then
            ii = i
            break
         endif
      enddo
      if (ii .lt. 0) then
         write (*,*) '*** mcrates could not find za=',za,' zn=',zn
         write (*,*) '*** check mcfilenames array'
         call xerrab("")
      endif
c
c     Compute rate parameters for transitions from table species ii.
c
      if (za .lt. zamax) then
         lrion = 
     .     (  fy)*((1-fxte)*rtlsa(i1e,j1+1,ii)+fxte*rtlsa(i1e+1,j1+1,ii))
     .   + (1-fy)*((1-fxte)*rtlsa(i1e,j1  ,ii)+fxte*rtlsa(i1e+1,j1  ,ii))
         rion = exp(lrion)
         if (za .eq. 0) return
      endif
         lrrec = 
     .     (  fy)*((1-fxte)*rtlra(i1e,j1+1,ii)+fxte*rtlra(i1e+1,j1+1,ii))
     .   + (1-fy)*((1-fxte)*rtlra(i1e,j1  ,ii)+fxte*rtlra(i1e+1,j1  ,ii))
         rrec = exp(lrrec)
         rcxr = exp((1-fxti)*rtlcx(i1i,0,ii)+fxti*rtlcx(i1i+1,0,ii))
c
c     Use special analytic fit for carbon c-x on neutral hydrogen.
c
      if ( (iscxfit .gt. 0) .and.
     .     (zn .eq. 6) .and. (za .le. zamax) ) then
         if (iscxfit.ge.1. .and. iscxfit.le.2.) then
            rcxr = (2.-iscxfit)*rcxr_zn6 (tmpi, za) + 
     .             (iscxfit-1.)*rcxr_zn6b(tmpi, za)
         endif
ccc         if (iscxfit .eq. 1) rcxr = rcxr_zn6 (tmpi, za)
ccc         if (iscxfit .eq. 2) rcxr = rcxr_zn6b (tmpi, za)
      endif
c
      return
      end

c-----End of subroutine mcrates----------------------------------------

      real function rcxr_zn6 (tmp, za)
      implicit none
      real tmp
      integer za
Use(Share)               # cutlo
Use(Physical_constants2) # ev2
c
c     Charge exchange rate parameter for carbon on neutral hydrogen.
c     Input (neutral hydrogen) temperature, tmp, is in [Joules/AMU].
c     Initial carbon charge state is za.
c     Output rate parameter (sigma-v) is in [m**3/sec].
c
c     This power-law fit was derived by Tom Rognlien from Figure 8.1
c     of a PhD thesis by C.F. Maggi on "Measurement and Interpretation
c     of Spectral Emission from JET Divertor Plasmas", January 1997,
c     JET report JET-IR(96)05.
c
c     local variables --
      real x, m0(6), m1(6), m2(6)
      data m0/-16.104,-18.27,-14.48,-14.85,-14.213,-17.576/
      data m1/0.5335,2.3657,0.05715,0.5219,0.42193,1.8758/
      data m2/-0.0009571,-0.29616,0.080947,0.048885,-0.033125,-0.095951/

      x = log10(max(tmp, cutlo)/ev2)
      rcxr_zn6 = 10**( m0(za) + m1(za)*x + m2(za)*x*x )
c
      return
      end
c
c-------End of function rcxr_zn6---------------------------------------

      real function rcxr_zn6b (tmp, za)
      implicit none
      real tmp
      integer za
Use(Share)               # cutlo
Use(Physical_constants2) # ev2
c
c     Charge exchange rate parameter for carbon on neutral hydrogen.
c     Input (neutral hydrogen) temperature, tmp, is in [Joules/AMU].
c     Initial carbon charge state is za.
c     Output rate parameter (sigma-v) is in [m**3/sec].
c
c     This is a modified of the function rcxr_zn6; only za=1 case is
c     changed to use a (lower) fit guided by plots from A. Pigarov.  
c     Other za's same as for rxcr_zn6 from thesis by C.F. Maggi (fit 
c     by T. Rognlien)
c
c     local variables --
      real x, m0(6), m1(6), m2(6)
      data m0/-20.027,-18.27,-14.48,-14.85,-14.213,-17.576/
      data m1/3.6433,2.3657,0.05715,0.5219,0.42193,1.8758/
      data m2/-0.59189,-0.29616,0.080947,0.048885,-0.033125,-0.095951/

      x = log10(max(tmp, cutlo)/ev2)
      rcxr_zn6b = 10**( m0(za) + m1(za)*x + m2(za)*x*x )
c
      return
      end
c
c-------End of function rcxr_zn6b---------------------------------------

      subroutine readmc(nzdf,mcfilename)

      implicit none
      integer nzdf
      character*256 mcfilename(*)
      character*256 fname
      Use(Multicharge)
      Use(Lsode)        # iprint
      Use(Impdata)      #apidir

c ... Function:
      integer utgetcl   # defined in the Basis runtime library

c     local variables --
      integer i, ios, kstart, n, nget, rtnt_old, rtnn_old, rtnsd_old, kk,
     .        stringlen1, stringlen2, iprt_imp_file
      character idcod*8, idtyp*8, id1*32, mcfnuse*256

c     procedures --
      external freeus, gchange, xerrab

c ... Initialize iprt_imp_file on
      data iprt_imp_file/1/

c----------------------------------------------------------------------c
c     Read rate data from 'un*formatted' files.
c           (file format from b2.5 code by B. Braams)
c----------------------------------------------------------------------c

      do i=1,nzdf
         rtnt_old  = rtnt
         rtnn_old  = rtnn
         rtnsd_old = rtnsd

      mcfnuse = mcfilename(i)
      call freeus(nget)
      fname = TRIM(apidir) //'/'//TRIM(mcfnuse)
      open (nget, file=TRIM(apidir)//'/'//TRIM(mcfnuse),
     .      form='formatted', iostat=ios, status='old')
      if (ios .ne. 0) then
         write(*,*)
     .     '*** Input file mcfilename = "', mcfilename(i),
     .     '" not found.'
         call xerrab("")
      endif

c     read header --
*     un*formatted read for header data
      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      if (n .lt. 0 .and. iprt_imp_file == 1) then  
        if (iprint .ne. 0) then
            write(*,*) '***Impurity file using new 2012 format is ',mcfilename(i)
        endif
        mcfformat(i) = 1
        iprt_imp_file = 0
      elseif (iprt_imp_file == 1) then
        if (iprint .ne. 0) then
            write(*,*) '***Impurity file using pre-2012 format is ',mcfilename(i)
        endif
        mcfformat(i) = 0
        iprt_imp_file = 0
      endif
      read (nget,'(1x,1a120)') labelrt(i)

c     read dimensions --
*     un*formatted read for integer data
      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) rtnt,rtnn,rtns

c     Test for compatibility of (rtnt,rtnn) from different tables:
      if ( (i .gt. 1) .and.
     .     ((rtnt .ne. rtnt_old) .or. (rtnn .ne. rtnn_old)) ) then
         write(*,*)
     .      '*** subroutine readmc: incompatible table dimensions in ',
     .      mcfilename(i),' and ',mcfilename(i-1)
         call xerrab("")
      endif

c     allocate storage --
      rtnsd=rtnsd+rtns
      call gchange("Multicharge",0)  #caution: can set vars, but not use here

c     read abscissae and rates --
      kstart=rtnsd_old		# starting species index for this table
      do kk = kstart, kstart+rtns-1
        chgstate_format(kk) = mcfformat(i)
      enddo

      call readmc1 (nget, kstart)

      close (nget)

      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine readmc1 (nget, kstart)
      implicit none
      integer nget, kstart
Use(Multicharge)

c     local variables --
      integer i, j, k, n
      character idcod*8, idtyp*8, id1*32

c     read abscissae --
*     un*formatted read for real data

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtza(k),k=kstart,kstart+rtns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtzn(k),k=kstart,kstart+rtns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtza2(k),k=kstart,kstart+rtns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtt(i),i=0,rtnt)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtn(j),j=0,rtnn)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtlt(i),i=0,rtnt)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtln(j),j=0,rtnn)

c     read rate coefficients --
*     un*formatted read for real data

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((rtlsa(i,j,k),i=0,rtnt),j=0,rtnn),
     .                              k=kstart,kstart+rtns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((rtlra(i,j,k),i=0,rtnt),j=0,rtnn),
     .                              k=kstart,kstart+rtns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((rtlqa(i,j,k),i=0,rtnt),j=0,rtnn),
     .                              k=kstart,kstart+rtns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((rtlcx(i,j,k),i=0,rtnt),j=0,rtnn),
     .                              k=kstart,kstart+rtns-1)

      return
      end

c----------------------------------------------------------------------c

      function radmc(zmax, znuc, te, dene, denz, radz)
      implicit none
      real radmc
c ... input args
      integer zmax, znuc
      real te, dene, denz(0:zmax)
c ... output args
      real radz(0:zmax)
c
c ... Compute the radiation rates, radz(0:zmax), for all charge states
c     of an impurity with nuclear charge, znuc, and return the total 
c     electron energy loss rate, radmc, including both the radiation
c     and binding energy contributions.
c
c     The tables used in this subroutine are generated with a code
c     supplied by Bas Braams.  The data file it produces is called
c     'b2frates' by default.  This gives rates that may depend on both
c     density and temperature.
c
c ... Input temperature te is in [J].
c     Input densities are in [/m**3].
c     Output rates are in [J/m**3/s].
c     Table temperatures are in [eV].
c     Table rate parameters are in [m**3/s] and [eV*m**3/s].
c
      Use(Share)               # cutlo
      Use(Physical_constants2) # ev2
      Use(Multicharge)         # rtln,rtlt,rtnn,rtnt,rtlqa,rtlsa,rtlra
                               # isrtndep,chgstate_format,ispradextrap
c
      real ebindz
      external ebindz
c
      integer i1,k,k0,j1
      real nenonz,y,dlogn,fy
      real tenonz,xt,dlogt,fxt
      real kionz,krecz,keelz,lionz,lrecz,leelz,fac_rad,temintab
c
      tenonz = max(te, cutlo)      # to avoid possible log(0) error
      xt = log(tenonz/ev2)
      dlogt = rtlt(1) - rtlt(0)

c ... Define minimum Te in table
      temintab = 0.2*ev2
c
c ... Find index i1 in table such that rtlt(i1) .le. xt .lt. rtlt(i1+1)
c     or, equivalently, rtt(i1) .le. te .lt. rtt(i1+1).
c
      i1 = int( (xt-rtlt(0))/dlogt )
c
c ... For temperatures below minimum table temperature, extrapolate
c     downwards from table entries 0 and 1.

      i1 = max(0, i1)
c
c ... For temperatures above maximum table temperature, extrapolate
c     upwards from table entries rtnt-1 and rtnt.

      i1 = min(rtnt-1, i1)
c
c ... Compute coefficient for linear interpolation.
c
      fxt = (xt-rtlt(i1))/(rtlt(i1+1)-rtlt(i1))
c
c ... Now the density dependence --
c     Default is to use the lowest density in the table -
      j1 = 0
      fy = 0.
c     Otherwise, use linear interpolation -
      if (isrtndep .ne. 0) then
         nenonz = max(dene, cutlo)
         y = log(nenonz)
         dlogn = rtln(1) - rtln(0)
         j1 = int( (y-rtln(0))/dlogn )
         j1 = max(0, j1)
         j1 = min(rtnn-1, j1)
         fy = (y-rtln(j1))/(rtln(j1+1)-rtln(j1))
c        Do not extrapolate beyond table minimum and maximum densities
         fy = max(0., fy)
         fy = min(1., fy)
      endif
c
c     Compute rates for each charge state and total rate:

c     First, find the species index of the neutral impurity atom
      k0 = -1
      do k=0,rtnsd-1
        if ((nint(rtzn(k)) .eq. znuc) .and. (nint(rtza(k)) .eq. 0)) then
           k0 = k
           break
        endif
      enddo
      if (k0 .lt. 0) then
         write (*,*) '*** radmc could not find za=',0,' zn=',znuc
         write (*,*) '*** check mcfilenames array'
         call xerrab("")
      endif

      radmc = 0.	# total electron energy loss rate
      do k=0,zmax
         lionz =
     .  (  fy)*((1-fxt)*rtlsa(i1,j1+1,k0+k)+fxt*rtlsa(i1+1,j1+1,k0+k))
     .+ (1-fy)*((1-fxt)*rtlsa(i1,j1  ,k0+k)+fxt*rtlsa(i1+1,j1  ,k0+k))
         kionz = exp(lionz)
         lrecz =
     .  (  fy)*((1-fxt)*rtlra(i1,j1+1,k0+k)+fxt*rtlra(i1+1,j1+1,k0+k))
     .+ (1-fy)*((1-fxt)*rtlra(i1,j1  ,k0+k)+fxt*rtlra(i1+1,j1  ,k0+k))
         krecz = exp(lrecz)
         leelz =
     .  (  fy)*((1-fxt)*rtlqa(i1,j1+1,k0+k)+fxt*rtlqa(i1+1,j1+1,k0+k))
     .+ (1-fy)*((1-fxt)*rtlqa(i1,j1  ,k0+k)+fxt*rtlqa(i1+1,j1  ,k0+k))
         keelz = exp(leelz)
         fac_rad = 1. 
         if(ispradextrap==1 .and. k==0 .and. te<temintab) then  #extrap below min Te
           fac_rad = (te/(temintab))**6
         endif
         radz(k) = fac_rad*dene*denz(k)*keelz*ev2
         radmc = radmc + radz(k)
c     binding energy contributions:
         if (chgstate_format(k0+k) == 0) then  #use pre-2012 format definitions
           if (k .lt. zmax) then		# ionization
              radz(k) = radz(k) - dene*denz(k)*kionz*ebindz(k,znuc)*ev2
           endif
           if (k .gt. 0) then		# recombination
              radz(k) = radz(k) + dene*denz(k)*krecz*ebindz(k-1,znuc)*ev2
           endif
         else                                  #use new 2012 format definitions
           if (k .lt. zmax) then		# ionization
             radmc = radmc + dene*denz(k)*kionz*ebindz(k,znuc)*ev2
           endif
           if (k .gt. 0) then		# recombination
             radmc = radmc - dene*denz(k)*krecz*ebindz(k-1,znuc)*ev2
           endif
         endif
      enddo

      return
      end

c----------------------------------------------------------------------c

      function ebindz(zatomic, znuclear)
      implicit none
      real ebindz
      integer zatomic, znuclear

c     This function returns the ionization energy (eV) for an atom
c     with atomic charge, zatomic, and nuclear charge, znuclear.
c     Data is from CRC Handbook of Chemistry and Physics, except as noted.

      if (zatomic .ge. znuclear) then
         write (*,*) '*** ebindz: input error'
         write (*,*) ' zatomic=',zatomic,'   .ge.   znuclear=',znuclear
         call xerrab("")
      endif

      if (znuclear .eq. 1) then		# hydrogen
         if (zatomic .eq. 0) then
            ebindz=13.59844
         endif
      elseif (znuclear .eq. 2) then	# helium
         if (zatomic .eq. 0) then
            ebindz=24.58741
         elseif (zatomic .eq. 1) then
            ebindz=54.41778
         endif
      elseif (znuclear .eq. 3) then	# lithium
         if (zatomic .eq. 0) then
            ebindz=5.39172
         elseif (zatomic .eq. 1) then
            ebindz=75.64018
         elseif (zatomic .eq. 2) then
            ebindz=122.45429
         endif
      elseif (znuclear .eq. 4) then	# beryllium
         if (zatomic .eq. 0) then
            ebindz=9.32263
         elseif (zatomic .eq. 1) then
            ebindz=18.21116
         elseif (zatomic .eq. 2) then
            ebindz=153.89661
         elseif (zatomic .eq. 3) then
            ebindz=217.71865
         endif
      elseif (znuclear .eq. 5) then	# boron
         if (zatomic .eq. 0) then
            ebindz=8.29803
         elseif (zatomic .eq. 1) then
            ebindz=25.15484
         elseif (zatomic .eq. 2) then
            ebindz=37.93064
         elseif (zatomic .eq. 3) then
            ebindz=259.37521
         elseif (zatomic .eq. 4) then
            ebindz=340.22580
         endif
      elseif (znuclear .eq. 6) then	# carbon
         if (zatomic .eq. 0) then
            ebindz=11.26030
         elseif (zatomic .eq. 1) then
            ebindz=24.38332
         elseif (zatomic .eq. 2) then
            ebindz=47.8878
         elseif (zatomic .eq. 3) then
            ebindz=64.4939
         elseif (zatomic .eq. 4) then
            ebindz=392.087
         elseif (zatomic .eq. 5) then
            ebindz=489.99334
         endif
      elseif (znuclear .eq. 7) then	# nitrogen
         if (zatomic .eq. 0) then
            ebindz=14.53414
         elseif (zatomic .eq. 1) then
            ebindz=29.6013
         elseif (zatomic .eq. 2) then
            ebindz=47.44924
         elseif (zatomic .eq. 3) then
            ebindz=77.4735
         elseif (zatomic .eq. 4) then
            ebindz=97.8902
         elseif (zatomic .eq. 5) then
            ebindz=552.0718
         elseif (zatomic .eq. 6) then
            ebindz=667.046
         endif
      elseif (znuclear .eq. 8) then	# oxygen
         if (zatomic .eq. 0) then
            ebindz=13.61806
         elseif (zatomic .eq. 1) then
            ebindz=35.11730
         elseif (zatomic .eq. 2) then
            ebindz=54.9355
         elseif (zatomic .eq. 3) then
            ebindz=77.41353
         elseif (zatomic .eq. 4) then
            ebindz=113.8990
         elseif (zatomic .eq. 5) then
            ebindz=138.1197
         elseif (zatomic .eq. 6) then
            ebindz=739.29
         elseif (zatomic .eq. 7) then
            ebindz=871.4101
         endif
      elseif (znuclear .eq. 9) then	# fluorine
         if (zatomic .eq. 0) then
            ebindz=17.42282
         elseif (zatomic .eq. 1) then
            ebindz=34.97082
         elseif (zatomic .eq. 2) then
            ebindz=62.7084
         elseif (zatomic .eq. 3) then
            ebindz=87.1398
         elseif (zatomic .eq. 4) then
            ebindz=114.2428
         elseif (zatomic .eq. 5) then
            ebindz=157.1651
         elseif (zatomic .eq. 6) then
            ebindz=185.186
         elseif (zatomic .eq. 7) then
            ebindz=953.9112
         elseif (zatomic .eq. 8) then
            ebindz=1103.1176
         endif
      elseif (znuclear .eq. 10) then	# neon
         if (zatomic .eq. 0) then
            ebindz=21.56454
         elseif (zatomic .eq. 1) then
            ebindz=40.96328
         elseif (zatomic .eq. 2) then
            ebindz=63.45
         elseif (zatomic .eq. 3) then
            ebindz=97.12
         elseif (zatomic .eq. 4) then
            ebindz=126.21
         elseif (zatomic .eq. 5) then
            ebindz=157.93
         elseif (zatomic .eq. 6) then
            ebindz=207.2759
         elseif (zatomic .eq. 7) then
            ebindz=239.0989
         elseif (zatomic .eq. 8) then
            ebindz=1195.8286
         elseif (zatomic .eq. 9) then
            ebindz=1362.1995
         endif
      elseif (znuclear .eq. 18) then	# argon
         if (zatomic .eq. 0) then
            ebindz=15.75962
         elseif (zatomic .eq. 1) then
            ebindz=27.62967
         elseif (zatomic .eq. 2) then
            ebindz=40.74
         elseif (zatomic .eq. 3) then
            ebindz=59.81
         elseif (zatomic .eq. 4) then
            ebindz=75.02
         elseif (zatomic .eq. 5) then
            ebindz=91.009
         elseif (zatomic .eq. 6) then
            ebindz=124.323
         elseif (zatomic .eq. 7) then
            ebindz=143.460
         elseif (zatomic .eq. 8) then
            ebindz=422.45
         elseif (zatomic .eq. 9) then
            ebindz=478.69
         elseif (zatomic .eq. 10) then
            ebindz=538.96
         elseif (zatomic .eq. 11) then
            ebindz=618.26
         elseif (zatomic .eq. 12) then
            ebindz=686.10
         elseif (zatomic .eq. 13) then
            ebindz=755.74
         elseif (zatomic .eq. 14) then
            ebindz=854.77
         elseif (zatomic .eq. 15) then
            ebindz=918.03
         elseif (zatomic .eq. 16) then
            ebindz=4120.8857
         elseif (zatomic .eq. 17) then
            ebindz=4426.2296
         endif
      elseif (znuclear .eq. 36) then	# krypton from CRC; state 30 and above
         if (zatomic .eq. 0) then       # are arbitrary values-need updating
            ebindz=13.99961
         elseif (zatomic .eq. 1) then
            ebindz=24.35985
         elseif (zatomic .eq. 2) then
            ebindz=36.950
         elseif (zatomic .eq. 3) then
            ebindz=52.5
         elseif (zatomic .eq. 4) then
            ebindz=64.7
         elseif (zatomic .eq. 5) then
            ebindz=78.5
         elseif (zatomic .eq. 6) then
            ebindz=111.0
         elseif (zatomic .eq. 7) then
            ebindz=125.802
         elseif (zatomic .eq. 8) then
            ebindz=230.85
         elseif (zatomic .eq. 9) then
            ebindz=268.2
         elseif (zatomic .eq. 10) then
            ebindz=308.
         elseif (zatomic .eq. 11) then
            ebindz=350.
         elseif (zatomic .eq. 12) then
            ebindz=391.
         elseif (zatomic .eq. 13) then
            ebindz=447.
         elseif (zatomic .eq. 14) then
            ebindz=492.
         elseif (zatomic .eq. 15) then
            ebindz=541.
         elseif (zatomic .eq. 16) then
            ebindz=592.
         elseif (zatomic .eq. 17) then
            ebindz=641.
         elseif (zatomic .eq. 18) then
            ebindz=786.
         elseif (zatomic .eq. 19) then
            ebindz=833.
         elseif (zatomic .eq. 20) then
            ebindz=884.
         elseif (zatomic .eq. 21) then
            ebindz=937.
         elseif (zatomic .eq. 22) then
            ebindz=998.
         elseif (zatomic .eq. 23) then
            ebindz=1051.
         elseif (zatomic .eq. 24) then
            ebindz=1151.
         elseif (zatomic .eq. 25) then
            ebindz=1205.3
         elseif (zatomic .eq. 26) then
            ebindz=2928.
         elseif (zatomic .eq. 27) then
            ebindz=3070.
         elseif (zatomic .eq. 28) then
            ebindz=3227.
         elseif (zatomic .eq. 29) then
            ebindz=3381.
         elseif (zatomic .eq. 30) then   # guesses to fill out table
            ebindz=3500.
         elseif (zatomic .eq. 31) then
            ebindz=3700.
         elseif (zatomic .eq. 32) then
            ebindz=3900.
         elseif (zatomic .eq. 33) then
            ebindz=4100.
         elseif (zatomic .eq. 34) then
            ebindz=4300.
         elseif (zatomic .eq. 35) then
            ebindz=4500.
         endif
      elseif (znuclear .eq. 42) then	# molybdenum from ADAS
         if (zatomic .eq. 0) then
            ebindz=5.737
         elseif (zatomic .eq. 1) then
            ebindz=14.56
         elseif (zatomic .eq. 2) then
            ebindz=29.15
         elseif (zatomic .eq. 3) then
            ebindz=43.37
         elseif (zatomic .eq. 4) then
            ebindz=57.59
         elseif (zatomic .eq. 5) then
            ebindz=71.81
         elseif (zatomic .eq. 6) then
            ebindz=124.0
         elseif (zatomic .eq. 7) then
            ebindz=145.8
         elseif (zatomic .eq. 8) then
            ebindz=167.5
         elseif (zatomic .eq. 9) then
            ebindz=189.2
         elseif (zatomic .eq. 10) then
            ebindz=216.7
         elseif (zatomic .eq. 11) then
            ebindz=238.7
         elseif (zatomic .eq. 12) then
            ebindz=288.1
         elseif (zatomic .eq. 13) then
            ebindz=312.5
         elseif (zatomic .eq. 14) then
            ebindz=505.3
         elseif (zatomic .eq. 15) then
            ebindz=570.8
         elseif (zatomic .eq. 16) then
            ebindz=636.2
         elseif (zatomic .eq. 17) then
            ebindz=701.7
         elseif (zatomic .eq. 18) then
            ebindz=767.1
         elseif (zatomic .eq. 19) then
            ebindz=832.6
         elseif (zatomic .eq. 20) then
            ebindz=902.2
         elseif (zatomic .eq. 21) then
            ebindz=968.4
         elseif (zatomic .eq. 22) then
            ebindz=1034.
         elseif (zatomic .eq. 23) then
            ebindz=1101.
         elseif (zatomic .eq. 24) then
            ebindz=1305.
         elseif (zatomic .eq. 25) then
            ebindz=1368.
         elseif (zatomic .eq. 26) then
            ebindz=1431.
         elseif (zatomic .eq. 27) then
            ebindz=1494.
         elseif (zatomic .eq. 28) then
            ebindz=1591.
         elseif (zatomic .eq. 29) then
            ebindz=1655.
         elseif (zatomic .eq. 30) then
            ebindz=1805.
         elseif (zatomic .eq. 31) then
            ebindz=1869.
         elseif (zatomic .eq. 32) then
            ebindz=3990.
         elseif (zatomic .eq. 33) then
            ebindz=4191.
         elseif (zatomic .eq. 34) then
            ebindz=4392.
         elseif (zatomic .eq. 35) then
            ebindz=4593.
         elseif (zatomic .eq. 36) then
            ebindz=4902.
         elseif (zatomic .eq. 37) then
            ebindz=5110.
         elseif (zatomic .eq. 38) then
            ebindz=5407.
         elseif (zatomic .eq. 39) then
            ebindz=5585.
         elseif (zatomic .eq. 40) then
            ebindz=23120.
         elseif (zatomic .eq. 41) then
            ebindz=23890.
         endif
      elseif (znuclear .eq. 50) then	# tin
                                        # Calculated by Jim Schofield at LLNL
                                        # (via Howard Scott on 13 Sep 2000)
         if (zatomic .eq. 0) then
            ebindz=5.8
c            ebindz=7.34381             # CRC handbook value
         elseif (zatomic .eq. 1) then
            ebindz=12.9
c            ebindz=14.63225            # CRC handbook value
         elseif (zatomic .eq. 2) then
            ebindz=30.1
c            ebindz=30.50260            # CRC handbook value
         elseif (zatomic .eq. 3) then
            ebindz=40.6
c            ebindz=40.73502            # CRC handbook value
         elseif (zatomic .eq. 4) then
            ebindz=76.4
c            ebindz=72.28               # CRC handbook value
         elseif (zatomic .eq. 5) then
            ebindz=96.0
         elseif (zatomic .eq. 6) then
            ebindz=116.5
         elseif (zatomic .eq. 7) then
            ebindz=137.9
         elseif (zatomic .eq. 8) then
            ebindz=160.2
         elseif (zatomic .eq. 9) then
            ebindz=183.3
         elseif (zatomic .eq. 10) then
            ebindz=208.7
         elseif (zatomic .eq. 11) then
            ebindz=233.4
         elseif (zatomic .eq. 12) then
            ebindz=258.7
         elseif (zatomic .eq. 13) then
            ebindz=284.8
         elseif (zatomic .eq. 14) then
            ebindz=382.1
         elseif (zatomic .eq. 15) then
            ebindz=410.5
         elseif (zatomic .eq. 16) then
            ebindz=439.4
         elseif (zatomic .eq. 17) then
            ebindz=468.6
         elseif (zatomic .eq. 18) then
            ebindz=509.7
         elseif (zatomic .eq. 19) then
            ebindz=540.3
         elseif (zatomic .eq. 20) then
            ebindz=615.0
         elseif (zatomic .eq. 21) then
            ebindz=647.6
         elseif (zatomic .eq. 22) then
            ebindz=1132.7
         elseif (zatomic .eq. 23) then
            ebindz=1204.6
         elseif (zatomic .eq. 24) then
            ebindz=1278.0
         elseif (zatomic .eq. 25) then
            ebindz=1352.8
         elseif (zatomic .eq. 26) then
            ebindz=1429.2
         elseif (zatomic .eq. 27) then
            ebindz=1507.0
         elseif (zatomic .eq. 28) then
            ebindz=1597.7
         elseif (zatomic .eq. 29) then
            ebindz=1679.0
         elseif (zatomic .eq. 30) then
            ebindz=1761.6
         elseif (zatomic .eq. 31) then
            ebindz=1845.6
         elseif (zatomic .eq. 32) then
            ebindz=2082.5
         elseif (zatomic .eq. 33) then
            ebindz=2157.8
         elseif (zatomic .eq. 34) then
            ebindz=2233.8
         elseif (zatomic .eq. 35) then
            ebindz=2310.2
         elseif (zatomic .eq. 36) then
            ebindz=2446.8
         elseif (zatomic .eq. 37) then
            ebindz=2526.1
         elseif (zatomic .eq. 38) then
            ebindz=2683.87
         elseif (zatomic .eq. 39) then
            ebindz=2762.49
         elseif (zatomic .eq. 40) then
            ebindz=6415.48
         elseif (zatomic .eq. 41) then
            ebindz=6627.3
         elseif (zatomic .eq. 42) then
            ebindz=6841.0
         elseif (zatomic .eq. 43) then
            ebindz=7055.9
         elseif (zatomic .eq. 44) then
            ebindz=7536.3
         elseif (zatomic .eq. 45) then
            ebindz=7762.9
         elseif (zatomic .eq. 46) then
            ebindz=8107.20
         elseif (zatomic .eq. 47) then
            ebindz=8306.99
         elseif (zatomic .eq. 48) then
            ebindz=34256.71
         elseif (zatomic .eq. 49) then
            ebindz=35192.32
         endif
      elseif (znuclear .eq. 74) then	# tungsten from ADAS
         if (zatomic .eq. 0) then
            ebindz=7.130
         elseif (zatomic .eq. 1) then
            ebindz=15.08
         elseif (zatomic .eq. 2) then
            ebindz=25.43
         elseif (zatomic .eq. 3) then
            ebindz=39.29
         elseif (zatomic .eq. 4) then
            ebindz=53.15
         elseif (zatomic .eq. 5) then
            ebindz=67.01
         elseif (zatomic .eq. 6) then
            ebindz=119.7
         elseif (zatomic .eq. 7) then
            ebindz=140.8
         elseif (zatomic .eq. 8) then
            ebindz=162.0
         elseif (zatomic .eq. 9) then
            ebindz=183.1
         elseif (zatomic .eq. 10) then
            ebindz=204.2
         elseif (zatomic .eq. 11) then
            ebindz=240.5
         elseif (zatomic .eq. 12) then
            ebindz=263.1
         elseif (zatomic .eq. 13) then
            ebindz=294.6
         elseif (zatomic .eq. 14) then
            ebindz=339.9
         elseif (zatomic .eq. 15) then
            ebindz=369.9
         elseif (zatomic .eq. 16) then
            ebindz=395.0
         elseif (zatomic .eq. 17) then
            ebindz=435.5
         elseif (zatomic .eq. 18) then
            ebindz=480.8
         elseif (zatomic .eq. 19) then
            ebindz=526.1
         elseif (zatomic .eq. 20) then
            ebindz=571.4
         elseif (zatomic .eq. 21) then
            ebindz=616.7
         elseif (zatomic .eq. 22) then
            ebindz=664.6
         elseif (zatomic .eq. 23) then
            ebindz=710.2
         elseif (zatomic .eq. 24) then
            ebindz=755.8
         elseif (zatomic .eq. 25) then
            ebindz=801.4
         elseif (zatomic .eq. 26) then
            ebindz=846.9
         elseif (zatomic .eq. 27) then
            ebindz=892.5
         elseif (zatomic .eq. 28) then
            ebindz=1154.
         elseif (zatomic .eq. 29) then
            ebindz=1206.
         elseif (zatomic .eq. 30) then
            ebindz=1259.
         elseif (zatomic .eq. 31) then
            ebindz=1312.
         elseif (zatomic .eq. 32) then
            ebindz=1365.
         elseif (zatomic .eq. 33) then
            ebindz=1417.
         elseif (zatomic .eq. 34) then
            ebindz=1483.
         elseif (zatomic .eq. 35) then
            ebindz=1537.
         elseif (zatomic .eq. 36) then
            ebindz=1591.
         elseif (zatomic .eq. 37) then
            ebindz=1645.
         elseif (zatomic .eq. 38) then
            ebindz=1870.
         elseif (zatomic .eq. 39) then
            ebindz=1926.
         elseif (zatomic .eq. 40) then
            ebindz=1981.
         elseif (zatomic .eq. 41) then
            ebindz=2037.
         elseif (zatomic .eq. 42) then
            ebindz=2163.
         elseif (zatomic .eq. 43) then
            ebindz=2223.
         elseif (zatomic .eq. 44) then
            ebindz=2386.
         elseif (zatomic .eq. 45) then
            ebindz=2447.
         elseif (zatomic .eq. 46) then
            ebindz=3734.
         elseif (zatomic .eq. 47) then
            ebindz=3882.
         elseif (zatomic .eq. 48) then
            ebindz=4029.
         elseif (zatomic .eq. 49) then
            ebindz=4177.
         elseif (zatomic .eq. 50) then
            ebindz=4325.
         elseif (zatomic .eq. 51) then
            ebindz=4472.
         elseif (zatomic .eq. 52) then
            ebindz=4684.
         elseif (zatomic .eq. 53) then
            ebindz=4836.
         elseif (zatomic .eq. 54) then
            ebindz=4987.
         elseif (zatomic .eq. 55) then
            ebindz=5139.
         elseif (zatomic .eq. 56) then
            ebindz=5538.
         elseif (zatomic .eq. 57) then
            ebindz=5671.
         elseif (zatomic .eq. 58) then
            ebindz=5803.
         elseif (zatomic .eq. 59) then
            ebindz=5936.
         elseif (zatomic .eq. 60) then
            ebindz=6468.
         elseif (zatomic .eq. 61) then
            ebindz=6611.
         elseif (zatomic .eq. 62) then
            ebindz=6919.
         elseif (zatomic .eq. 63) then
            ebindz=7055.
         elseif (zatomic .eq. 64) then
            ebindz=14760.
         elseif (zatomic .eq. 65) then
            ebindz=15140.
         elseif (zatomic .eq. 66) then
            ebindz=15520.
         elseif (zatomic .eq. 67) then
            ebindz=15900.
         elseif (zatomic .eq. 68) then
            ebindz=17630.
         elseif (zatomic .eq. 69) then
            ebindz=18060.
         elseif (zatomic .eq. 70) then
            ebindz=18800.
         elseif (zatomic .eq. 71) then
            ebindz=19150.
         elseif (zatomic .eq. 72) then
            ebindz=77510.
         elseif (zatomic .eq. 73) then
            ebindz=78990.
         endif
      else   # data not available
         write (*,*) '*** ebindz: no binding energy data'
         write (*,*) '    for znuclear=',znuclear,', zatomic=',zatomic
         call xerrab("")
      endif

      return
      end

c----------------------------------------------------------------------c

