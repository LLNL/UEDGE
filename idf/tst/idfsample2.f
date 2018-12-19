       program couettes  
c  The second sample program in FORTRAN  

c     initialization of input IDF data file String 
       character*(*)  Fname
       integer Flen 
       parameter  (Flen = 12 , Fname = 'couetteS.dat') 

c     declaration of input parameters within a common block 
       common/COUETTE/Lz,Knu,Uw,Tw,Pr,Nz,Nv 
       real*8     Lz 
       real*4     Knu,Uw,Tw 
       real*8     Pr(3) 
       integer*4  Nz 
       integer*4  Nv(2) 
       integer*4  ier, couetinp 

c     import data with IDF package 
       ier = couetinp(Fname,Flen) 

         if(ier.eq.0)  then 
c        print the imported data sets 
       print  10,Lz,Knu 
10     format('Lz=',e12.3, ' Knu=',e12.3) 
       print  20,Uw,Tw 
20     format('Uw=',e12.3, ' Tw=',e12.3) 
       print  30,Pr(1),Pr(2),Pr(3) 
30     format ('Pr=',e12.3, ',' ,e12.3, ',', e12.3) 
       print  40,Nz,Nv(1),Nv(2) 
40     format('Nz=',I5, ' Nv=',I5,',', I5) 
         end if 

       stop  
       end

  
       integer function couetinp(namef,lenf)
       character*(*) namef
       integer*4 lenf

c      function returns zero in success

       common/COUETTE/Lz,Knu,Uw,Tw,Pr,Nz,Nv 
       real*8     Lz 
       real*4     Knu,Uw,Tw 
       real*8     Pr(3) 
       integer*4  Nz 
       integer*4  Nv(2) 

       integer*4 idfinit,idfopen,idfget

c      initialization of Align Mode 
c      for F77 compiler with '-align dcommon' option
       integer amode
       parameter (amode=7)

c      initialization of Target String and Data Set Name
       integer lSN, lTS
       parameter (lSN=7 , lTS=13)
       character*13 TS /'dfffd[3]ii[2]'/
       character*7 SN /'COUETTE'/

c      start the IDF package
          couetinp=idfinit(0)
       if(couetinp.ne.0) go to 1

c      open the data file
          couetinp=idfopen(namef,lenf)
       if(couetinp.ne.0) go to 2

c      transfer the data from file to common block
          couetinp = idfget(SN, lSN, TS, lTS, amode, Lz)

       if(couetinp) 2,3,3
3      couetinp=0

c      close the data file
2      call idfclose

c      finish the IDF
1      call idfinish

       return
       end
