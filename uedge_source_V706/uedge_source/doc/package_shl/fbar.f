	integer function fbar(a,len)
	real*4 a(*)
	integer len,i
	do i=1,len
		a(i) = i
	enddo
	fbar = len
	return
	end
