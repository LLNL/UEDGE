ccc      FUNCTION SECOND()
ccc       SECOND = 0.
ccc      return
ccc      end

      subroutine basfilex(name, namex)
c ... Basis uses this to expand special characters ($,~); here just set equal
      character*(*) name, namex
      namex = name
      return
      end

      subroutine glbwrlog(ioun)
      integer ioun
      return
      end
cc *********************************************************** cc
	function basopen(name,spec)
c
c                              basopen
c
	character*(*) name, spec
	integer basopen
        integer fd
	
c open a file as a read or write file, or see if you could.
	if(spec .eq. 'w') then
		call freeus(fd)
		open(unit=fd,file=name,err=400,status='new')
	elseif(spec .eq. 'r') then
		call freeus(fd)
		open(unit=fd,file=name,err=300,status = "old")
	else
		write (6,*) 'unknow spec for basopen'
                stop
	endif
        basopen = fd
        write (6,*) 'OPENING UNIT ', fd
	return
300     continue
        write (6,*) ('basopen: cannot open file:')
        stop
400     continue
        write (6,*) ('basopen: cannot create file:')
        stop
	end

cc ************************************************************* cc
	subroutine basclose(us)
        integer iotable(999)
        common/iotable/iotable
        
	integer us
	iotable(us) = 0
	close(unit=us,err=100)
        write (6,*) 'CLOSING UNIT ', us
	return
100 	write(6,200) us
200 	format('basclose: error in attempting to close unit ',i6,'.')
        stop
	end

cc ************************************************************* cc
        subroutine freeus(ius)
        integer iotable(999)
        common/iotable/iotable
        integer i,ius
        do i = 10, 500
           if(iotable(i) .eq. 0) then
              ius = i
              iotable(i) = 1
              if(i==500) iotable(501:999)=0
              return
           endif
        enddo
        do i = 501, 999
           if(iotable(i) .eq. 0) then
              ius = i
              iotable(i) = 1
              if(i==999) iotable(10:500)=0
              return
           endif
        enddo
c---haven't found an available unit
cSEK        write (6, *) "You have all 999 I/O units in use"
cSEK        stop
        end
      
cc ************************************************************* cc
        subroutine glbwrtim(unit,cpu,io,sys,mem)
        integer unit
        real cpu,io,sys,mem
        
        return
        end
