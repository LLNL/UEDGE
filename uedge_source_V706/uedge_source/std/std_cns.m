      subroutine std_cnst

      implicit none

      Use (std_cns)

      pi0       =dacos(-1.d0)
      invpi     =1.d0/pi0
      twopi     =2.d0*pi0
      pi4       =4.d0*pi0
      pisqrt    =dsqrt(pi0)
      invpisqrt =1.d0/pisqrt
      u2pi      =pi0/180.d0
      pi2u      =180.d0/pi0

      aln10     =dlog(10.d0)
      onethird  =1.d0/3.d0
      twothirds =2.d0/3.d0

      ev2erg    =1.602176462d-12
      mol1      =1.66053873d0          ! amu*10^24 or 10^24/mole 
      mol1nrm   =mol1*1.d-24           ! 1/mole
      me        =5.4858006d-4          ! electron mass, amu
      menrm     =me*mol1nrm

      end

