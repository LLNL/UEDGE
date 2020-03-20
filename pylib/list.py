# Print information on a UEDGE variable                                                                                      
# Usage:                                                                                                                     
# >>> list("itermx")                                                                                                         
#                                                                                                                            
# First coding: MVU 20-mar-2020                                                                                                            
#============================================================                                                                

def list(varname):
    resbbb=bbb.listvar(varname)
    if (resbbb != ''):
        print(resbbb)
#                                                                                                                            
    rescom=com.listvar(varname)
    if (rescom != ''):
        print(rescom)
#                                                                                                                            
    ressvr=svr.listvar(varname)
    if (ressvr != ''):
        print(ressvr)
#                                                                                                                            
    resflx=flx.listvar(varname)
    if (resflx != ''):
        print(resflx)
#                                                                                                                            
    resgrd=grd.listvar(varname)
    if (resgrd != ''):
        print(resgrd)
        
