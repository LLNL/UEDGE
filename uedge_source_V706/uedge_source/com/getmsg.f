      integer function ijmgetmr(msg,max,iwait,nbrcv)
c     
      character*(*) msg
      integer max, iwait, nbrcv
      character*80 chi
c     integer input(10)
      character*80 input
      equivalence(input,chi)
      integer i, rdfile, fcntl
      integer f_setfl, f_getfl, o_ndelay, fildes, args, curr_fl
      data fildes /0/, f_setfl /4/, f_getfl /3/, o_ndelay /4/, args /0/
c     
c     fildes=0
c     f_setfl=4
c     o_ndelay=4
c  GET CURRENT FLAGS FOR STDIN
      curr_fl=fcntl(fildes,f_getfl,args)
c     write(*,*) 'In ijmgetmr; curr_fl = ', curr_fl 
c  TURN ON O_NODELAY
      args = o_ndelay
ccc      args = curr_fl .or. o_ndelay
ccc TDR 5/26/03  replace previous line with similar that avoids some compile
c     write (*,*) 'In ijmgetmr (fcntl); args = ', args
      i=fcntl(fildes,f_setfl,args)
c     write(*,*) 'In ijmgetmr; i = ', i
      if (i.eq.-1) then
         write (102,*) "call to fcntl failed"
         stop 'fcntl'
      endif
c  DO READ 
      input = 'oops'
      i=rdfile(fildes,input,80)
c     write(*,*) 'In ijmgetmr (rdfile); i = ', i
c     if (i.eq.-1) then
c        write(102,*) "call to read failed"
c        ijmgetmr = 2
c        stop 'read'
c     elseif (i.eq.0) then
      if (i.le.0) then
	 ijmgetmr = 1
      else
	 msg = chi(1:i-1)
	 nbrcv = i-1
	 ijmgetmr = 0
      endif
c  RESET FLAGS   
      i=fcntl(fildes,f_setfl,curr_fl)
c     write(*,*) 'In ijmgetmr (reset); i = ', i
      return
      end
