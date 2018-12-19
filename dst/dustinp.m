      subroutine dstinp
      
      implicit none

c      Use (dustdim)
      Use (dustinp)
      Use (dustcntrl)
      Use (std_cns)

      integer i
      integer Flen,ierr
      character*256 Fname
      character*(*) Fnam,Ddir
      integer idfinit,idfopen,idfi,idfs,idfd
      parameter ( Ddir='./dat/',Fnam='dustt.inp' )

      outdir='./out/'
      Fname=Ddir//Fnam
      Flen=len_trim(Fname)
c    start the IDF package
      ierr=idfinit(0)
       if (ierr .ne. 0) goto 1

c    open the data file
      ierr=idfopen(Fname,Flen)
       if (ierr .ne. 0) goto 1

      ierr=idfi('jtheps',6,jtheps)
       if (ierr .ne. 0) goto 1
      ierr=idfi('jcol',4,jcol)
       if (ierr .ne. 0) goto 1
      ierr=idfi('jthe',4,jthe)
       if (ierr .ne. 0) goto 1
      ierr=idfi('jsee',4,jsee)
       if (ierr .ne. 0) goto 1
      ierr=idfi('jdsr',4,jdsr)
       if (ierr .ne. 0) goto 1
      ierr=idfi('jrds',4,jrds)
       if (ierr .ne. 0) goto 1
      ierr=idfi('jrdf',4,jrdf)
       if (ierr .ne. 0) goto 1
      ierr=idfi('jimp',4,jimp)
       if (ierr .ne. 0) goto 1
      ierr=idfi('jsrc',4,jsrc)
       if (ierr .ne. 0) goto 1
      ierr=idfi('jstt',4,jstt)
       if (ierr .ne. 0) goto 1
      
      if (jsrc .le. 0) jimp=0

      ierr=idfs('fnamgmt',7,fnamgmt,256)
       if (ierr .ne. 0) goto 1
      ierr=idfs('fnamplz',7,fnamplz,256)
       if (ierr .ne. 0) goto 1
      ierr=idfs('fnamimp',7,fnamimp,256)
       if (ierr .ne. 0) goto 1
      fnamgmt=Ddir(1:len_trim(Ddir))//fnamgmt
      fnamplz=Ddir(1:len_trim(Ddir))//fnamplz
      fnamimp=Ddir(1:len_trim(Ddir))//fnamimp

      ierr=idfi('jtrj',4,jtrj)
       if (ierr .ne. 0) goto 1
      ierr=idfi('jrartrj',7,jrartrj)
       if (ierr .ne. 0) goto 1
      ierr=idfi('jtrjfm',6,jtrjfm)
       if (ierr .ne. 0) goto 1
      ierr=idfs('namfo',5,namfo,64)
       if (ierr .ne. 0) goto 1
      namfo=outdir(1:len_trim(outdir))//namfo

      ierr=idfi('zip',3,zip)
       if (ierr .ne. 0) goto 1
      ierr=idfi('zmat',4,zmat)
       if (ierr .ne. 0) goto 1
       
      ierr=idfi('winisp',6,winisp)
       if (ierr .ne. 0) goto 1
      ierr=idfi('wonisp',6,wonisp)
       if (ierr .ne. 0) goto 1
      ierr=idfi('dibnisp',7,dibnisp)
       if (ierr .ne. 0) goto 1
      ierr=idfi('dobnisp',7,dobnisp)
       if (ierr .ne. 0) goto 1
      ierr=idfi('pibnisp',7,pibnisp)
       if (ierr .ne. 0) goto 1
      ierr=idfi('pobnisp',7,pobnisp)
       if (ierr .ne. 0) goto 1
      ierr=idfi('ditnisp',7,ditnisp)
       if (ierr .ne. 0) goto 1
      ierr=idfi('dotnisp',7,dotnisp)
       if (ierr .ne. 0) goto 1
      ierr=idfi('pitnisp',7,pitnisp)
       if (ierr .ne. 0) goto 1
      ierr=idfi('potnisp',7,potnisp)
       if (ierr .ne. 0) goto 1
       
      ierr=idfi('cells',5,cells)
       if (ierr .ne. 0) goto 1
      ierr=idfi('ixs',3,ixs)
       if (ierr .ne. 0) goto 1
      ierr=idfi('iys',3,iys)
       if (ierr .ne. 0) goto 1
       
      ierr=idfi('inipnt',6,inipnt)
       if (ierr .ne. 0) goto 1
      ierr=idfi('inivel',6,inivel)
       if (ierr .ne. 0) goto 1
      ierr=idfd('unorm0',6,unorm0)
       if (ierr .ne. 0) goto 1
      ierr=idfd('uphi0',5,uphi0)
       if (ierr .ne. 0) goto 1
      ierr=idfd('rini',4,rini)
       if (ierr .ne. 0) goto 1
      ierr=idfd('zini',4,zini)
       if (ierr .ne. 0) goto 1
      ierr=idfd('vrini',5,vrini)
       if (ierr .ne. 0) goto 1
      ierr=idfd('vzini',5,vzini)
       if (ierr .ne. 0) goto 1
      ierr=idfd('vtini',5,vtini)
       if (ierr .ne. 0) goto 1
      ierr=idfi('irndsk',6,irndsk)
       if (ierr .ne. 0) goto 1
                   
      ierr=idfd('rblb',4,rblb)
       if (ierr .ne. 0) goto 1
      ierr=idfd('vblb',4,vblb)
       if (ierr .ne. 0) goto 1
      ierr=idfd('machb',5,machb)
       if (ierr .ne. 0) goto 1

      ierr=idfd('albws',5,albws)
       if (ierr .ne. 0) goto 1
      ierr=idfd('mlbws',5,mlbws)
       if (ierr .ne. 0) goto 1
      ierr=idfd('rnvs',4,rnvs)
       if (ierr .ne. 0) goto 1
      ierr=idfd('mrvs',4,mrvs)
       if (ierr .ne. 0) goto 1
      ierr=idfd('tcs',3,tcs)
       if (ierr .ne. 0) goto 1
      ierr=idfd('albwl',5,albwl)
       if (ierr .ne. 0) goto 1
      ierr=idfd('mlbwl',5,mlbwl)
       if (ierr .ne. 0) goto 1
      ierr=idfd('rnvl',4,rnvl)
       if (ierr .ne. 0) goto 1
      ierr=idfd('mrvl',4,mrvl)
       if (ierr .ne. 0) goto 1
      ierr=idfd('tcl',3,tcl)
       if (ierr .ne. 0) goto 1

      ierr=idfd('gsplit',6,gsplit)
       if (ierr .ne. 0) goto 1

      ierr=idfd('rd0',3,rd0)
       if (ierr .ne. 0) goto 1
      ierr=idfd('vp0',3,vp0)
       if (ierr .ne. 0) goto 1
      ierr=idfd('td0',3,td0)
       if (ierr .ne. 0) goto 1
      ierr=idfd('xliqud0',7,xliqud0)
       if (ierr .ne. 0) goto 1

      ierr=idfd('radst',5,radst)
       if (ierr .ne. 0) goto 1
      ierr=idfd('wrd',3,wrd)
       if (ierr .ne. 0) goto 1
      ierr=idfd('rpow',4,rpow)
       if (ierr .ne. 0) goto 1
      ierr=idfd('minrd',5,minrd)
       if (ierr .ne. 0) goto 1
      ierr=idfd('maxrd',5,maxrd)
       if (ierr .ne. 0) goto 1

      minlogrd=dlog(minrd)/aln10
      maxlogrd=dlog(maxrd)/aln10
      fdlogrd=(maxlogrd-minlogrd)/ffrdim
      dlogrd=(maxlogrd-minlogrd)/frdim

      ierr=idfd('vdt',3,vdt)
       if (ierr .ne. 0) goto 1
      ierr=idfd('maxvd',5,maxvd)
       if (ierr .ne. 0) goto 1
      ierr=idfd('unorm1',6,unorm1)
       if (ierr .ne. 0) goto 1
      ierr=idfd('uphi1',5,uphi1)
       if (ierr .ne. 0) goto 1
      ierr=idfd('dirvc',5,dirvc)
       if (ierr .ne. 0) goto 1
       
      fdvd = maxvd/ffvdim
      dvd  = maxvd/fvdim

      lambda = 2.5d0   ! normalized dust potential to start itterations
      
      ierr=idfd('ksip',4,ksip)
       if (ierr .ne. 0) goto 1
      ierr=idfd('ksiq',4,ksiq)
       if (ierr .ne. 0) goto 1
      ierr=idfd('ksifp',5,ksifp)
       if (ierr .ne. 0) goto 1
      ierr=idfd('ksifa',5,ksifa)
       if (ierr .ne. 0) goto 1
      ierr=idfd('ksirn',5,ksirn)
       if (ierr .ne. 0) goto 1
      ierr=idfd('ksicp',5,ksicp)
       if (ierr .ne. 0) goto 1
      ierr=idfd('ksib',4,ksib)
       if (ierr .ne. 0) goto 1

      ierr=idfd('tensf',5,tensf)
       if (ierr .ne. 0) goto 1
      ierr=idfd('stens',5,stens)
       if (ierr .ne. 0) goto 1
       
      ierr=idfd('tw0',3,tw0)
       if (ierr .ne. 0) goto 1
      
      tw04=tw0*tw0*tw0*tw0
      tw05=tw04*tw0
      
      ierr=idfd('dstsc',5,dstsc)
       if (ierr .ne. 0) goto 1
      ierr=idfd('dstrc',5,dstrc)
       if (ierr .ne. 0) goto 1
      ierr=idfd('flxic',5,flxic)
       if (ierr .ne. 0) goto 1
      
      stp0=1.d-16  ! minimal timestep

1     if (ierr .ne. 0) call idferprn

      call idfclose
      call idfinish

      if (ierr .ne. 0) return  ! stop

      end
