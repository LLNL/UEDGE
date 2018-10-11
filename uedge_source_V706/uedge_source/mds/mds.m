ifdef([NO_INCLUDES],[],[
 include mds.mac
])



	subroutine setmdserror(err)
	integer err
Use(MDSVars)
	mdsplus_error = err
	return
	end
	subroutine getmdserror(err)
	integer err
Use(MDSVars)
        err = mdsplus_error
	return
	end
	subroutine setmdssocket(sock)
	Integer sock
Use(MDSVars)
	mdsplus_socket = sock
	return
	end
	subroutine getmdssocket(sock)
	Integer sock
Use(MDSVars)
	sock = mdsplus_socket
	return
	end

