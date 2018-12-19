       program couette  
c  The first sample program in FORTRAN  

c    initialization of input IDF data file String  
       character*(*)  Fname 
       integer Flen 
       parameter (Flen = 11 , Fname = 'couette.dat')

c    declaration of input parameters  
       real*8  Lz,Pr(3)  
       real*4  Knu,Uw,Tw  
       integer*4  Nz,Nv(2)  
       integer*4  ier  

c    start the IDF package  
          ier=idfinit(0)  
       if(ier.ne.0) go to 1  

c    open the data file  
          ier=idfopen(Fname,Flen)  
       if(ier.ne.0) go to 2  

c    reading the data from file  
          ier=idfd('Lz', 2, Lz)  
       if(ier.ne.0)  go to  2  
          ier=idff('Knu', 3, Knu)  
       if(ier.ne.0)  go  to  2  
          ier=idff('Uw', 2, Uw)  
       if(ier.ne.0)  go  to  2  
          ier=idff('Tw', 2,  Tw)  
       if(ier.ne.0)  go  to  2  
          ier=idfdarr('Pr', 2, Pr, 3)  
       if(ier.ne.3)  go  to  2  
          ier=idfi('Nz', 2, Nz)  
       if(ier.ne.0)  go  to  2  
          ier=idfiarr('Nv', 2, Nv, 2)  
       if(ier.ne.2)  go  to  2  

c    print the imported data sets  
       print 10,Lz,Knu  
 10    format('Lz=',E12.3, 'Knu=',E12.3)  
       print 20,Uw,Tw  
 20    format('Uw=',E12.3, '  Tw=',E12.3)  
       print 30,Pr(1),Pr(2),Pr(3)  
 30    format ('Pr=',E12.3, ',', E12.3, ',', E12.3)  
       print 40,Nz,Nv(1),Nv(2)  
 40    format('Nz=',I5, 'Nv=',I5,',', I5)  

c    close the data file  
 2     call idfclose  

c    finish the IDF  
 1     call idfinish  
       stop  
       end  
