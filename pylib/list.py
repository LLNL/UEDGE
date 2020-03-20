# Print information on a UEDGE variable                                                                                      
# Usage:                                                                                                                     
# >>> list("itermx")                                                                                                         
#                                                                                                                            
# MVU 20-mar-2020                                                                                                            
#============================================================                                                                

def list(varname):
    resbbb=bbb.listvar(varname)
    if (resbbb != None):
        print(resbbb)
#                                                                                                                            
    rescom=com.listvar(varname)
    if (rescom != None):
        print(rescom)
#                                                                                                                            
    ressvr=svr.listvar(varname)
    if (ressvr != None):
        print(ressvr)
#                                                                                                                            
    resflx=flx.listvar(varname)
    if (resflx != None):
        print(resflx)
#                                                                                                                            
    resgrd=grd.listvar(varname)
    if (resgrd != None):
        print(resgrd)
