      subroutine psitabtheps
      implicit none
      
      Use (const)
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      Use (psitabmxx)

      integer i,j
      integer m,inf
      real*8 eps,quadr,ucut
      real*8 emiss
      integer fout
      parameter (fout=2123)
      
      real*8 rad,temp
      integer iz
      common /dust1/ rad,temp
      common /dust2/ iz
      
      real*8 radtmp,tmptmp
      
      integer imat
      real*8 fomegap_mat
      external fomegap_mat
      external radiation

      do 100, imat=1,Nmat_pt
       iz=IZmatter(imat)

       call psi_tmesh_ini
c----
       do 10, i=1,Nrd_pt
        radtmp=dlog(rdmin_pt*(rdmax_pt/rdmin_pt)
     *                       **(dble(i-1)/dble(Nrd_pt-1)))
        rdk(i)=dexp(radtmp)
        rad=rdk(i)*1.d-6 ! micron-->m
       
        do 20, j=1,Ntde_pt
         tmptmp=dlog(Tdmine_pt*(Tdmaxe_pt/Tdmine_pt)
     *                         **(dble(j-1)/dble(Ntde_pt-1)))
         tkkk(j)=dexp(tmptmp)
         temp=tkkk(j)
         m=1000
         eps=1.d-4
         inf=4
cc         inf=2
         quadr=0.d0
         ucut=hplanck*fomegap_mat(iz)/(4.d0*pi*kboltz*temp)
         call inthp(0.d0,ucut,1.d1,radiation,m,0.d0,eps,inf,quadr)
ccc         call inthp(0.d0,1.d0,1.d1,radiation,m,0.d0,eps,inf,quadr)
         emiss=quadr/(sigma*temp**4)
         gytheps_mx(i,j)=emiss
20      continue
10     continue

       do 2, j=1,Ntde_pt
        do 1, i=1,Nrd_pt
         dtheps_mx(i,j,imat)=gytheps_mx(i,j)
1       continue
2      continue

100   continue

      end
      
