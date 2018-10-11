c     ------------------------------------------------------------------
c
c   These are just fortran wrappers to the C library c_dce.c
c   The f_to_c_string converts fortran character strings
c   to C char arrays.
c
c  Created:Bill Meyer 12/11/92
c
c

	subroutine dceinit
Use(DCE_vars)
	integer status
	integer iserver(100)
	call f_to_c_string(dce_server,iserver)
	call c_dceinit(iserver,status)
	return
	end

	subroutine f_init_dtm(sserver,dtmport)
	character*(*) sserver
	integer dtmport
	integer iserver(100)

	call f_to_c_string(sserver,iserver)
	call c_init_dtm(iserver,dtmport)

    	return
	end

	subroutine f_init_dce(server,status)
	character*(*) server
	integer iserver(100)
	integer status

	call f_to_c_string(server,iserver)
        call c_dceinit(iserver,status)

    	return
	end

	subroutine f_close_dtm(dtmport,status)
	integer dtmport
	integer status

	call c_close_dtm(dtmport,status)

    	return
	end

	subroutine f_dce_timeout(i,status)
	integer i
	integer status

	call c_dce_timeout(i,status)

	return
	end


	subroutine f_senddtm(port,stitle,f,ix,iy,status)
	integer port
	character *(*) stitle
	integer ititle(100)
	real f
	integer ix,iy
	integer status

	call f_to_c_string(stitle,ititle)
	call c_senddtm(port,ititle,f,ix,iy,status)


	return
	end

	subroutine f_to_c_string(fstring,cstring)
	character*(*) fstring
	integer cstring(*)
	character*400 cd
	integer id(100)
	integer j
	equivalence (cd,id)

	cd = fstring // char(0)
	do j = 1,99
		cstring(j) = id(j)
	enddo

	return
	end


	subroutine f_xgraph(sparams,x,f,ix,status)
	character*(*) sparams
	integer iparams(100)
	real f(1)
	real x(1)
	integer ix
	integer status

	call f_to_c_string(sparams,iparams)
	call c_xgraph(iparams,x,f,ix,status)

	return
	end


	subroutine f_senddtm_rzxform(port,stitle,f,r,z,xlen,ylen,oxlen,oylen, 
     x           rmin,rmax,zmin,zmax,status)
	integer port
	character *(*) stitle
	integer ititle(100)
	real f,r,z
	integer xlen,ylen,oxlen,oylen
	real rmin,rmax,zmin,zmax
	integer status
	
	call f_to_c_string(stitle,ititle)
	call c_senddtm_rzxform(port,ititle,f,r,z,xlen,ylen,oxlen,
     x 	oylen,rmin,rmax,zmin,zmax,status)
	return
	end

	subroutine f_interpdtm_rzxform(port,stitle,f,r,z,vint,rint,zint,numint,xlen,ylen,
     x           oxlen,oylen,rmin,rmax,zmin,zmax,status)
	integer port
	character *(*) stitle
	integer ititle(100)
	real f,r,z,vint,zint,rint
	integer numint,xlen,ylen,oxlen,oylen
	real rmin,rmax,zmin,zmax
	integer status
	
	call f_to_c_string(stitle,ititle)
	call c_interpdtm_rzxform(port,ititle,f,r,z,vint,rint,zint,numint,xlen,ylen,oxlen,
     x 	oylen,rmin,rmax,zmin,zmax,status)
	return
	end

	subroutine f_interprzxform(out,f,r,z,vint,rint,zint,numint,xlen,ylen,oxlen,oylen, 
     x           rmin,rmax,zmin,zmax,status)
	real vint,rint,zint
	integer numint
	real out,f,r,z
	integer xlen,ylen,oxlen,oylen
	real rmin,rmax,zmin,zmax
	integer status
	
	call c_interprzxform(out,f,r,z,vint,rint,zint,numint,xlen,ylen,oxlen,
     x 	oylen,rmin,rmax,zmin,zmax,status)
	return
	end

	subroutine f_rzxform(out,f,r,z,xlen,ylen,oxlen,oylen, 
     x           rmin,rmax,zmin,zmax,status)
	real out,f,r,z
	integer xlen,ylen,oxlen,oylen
	real rmin,rmax,zmin,zmax
	integer status
	
	call c_rzxform(out,f,r,z,xlen,ylen,oxlen,
     x 	oylen,rmin,rmax,zmin,zmax,status)
	return
	end
