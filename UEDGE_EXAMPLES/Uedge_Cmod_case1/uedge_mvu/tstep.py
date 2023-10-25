from uedge import bbb, com
import numpy as np
import h5py



def uestep(dtnow, depth_max=10, nrefine=3, debug=False, reset=False):
# Usage example, to integrate by deltat=1.0 s, do:
# uestep(1.0, reset=True)
#======================================================#
        
    if (debug):
        print("In uedge_step() depth, depth_max, total time:", uestep.depth, depth_max, uestep.deltat)

    if (reset):
        uestep.abort=False
        uestep.depth=0
        uestep.deltat=0.0        

    if (uestep.abort):
        print("aborting")
        return
    else:
        if (uestep.depth>depth_max):
            print("Too many recursions")
            uestep.abort=True
            pass
            #return
        else:

            print("Trying time step: ", dtnow)
            bbb.dtreal=dtnow
            bbb.exmain()

            if (bbb.iterm!=1):
                #-previous step failed, try substeps now
                uestep.depth += 1
                dtnew=dtnow/nrefine
                for i in range(1,nrefine+1):
                    uestep(dtnew, depth_max, nrefine, debug)
                #-on completion of all substeps return to previous level
                uestep.depth -= 1                
            else:
                print("Successful time step with dt=", dtnow)
                uestep.deltat=uestep.deltat+dtnow
                print("Cumulative deltat=", uestep.deltat)
                print("\n")
                pass
        
    #return ##-this one does not matter?

#-attributes for the function, used as static variables
uestep.abort=False
uestep.depth=0
uestep.deltat=0.0




def run_time_evol(dt=1e0, nt=1, save=False, filename="thist.h5"):
#-run a sequence of time steps of the same size

    ni_hist = np.zeros([com.nx+2, com.ny+2, 2, nt+1])
    up_hist = np.zeros([com.nx+2, com.ny+2, 2, nt+1])
    te_hist = np.zeros([com.nx+2, com.ny+2, nt+1])
    ti_hist = np.zeros([com.nx+2, com.ny+2, nt+1])

    ni_hist[:,:,:,0]=bbb.ni
    up_hist[:,:,:,0]=bbb.up
    te_hist[:,:,0]=bbb.te
    ti_hist[:,:,0]=bbb.ti


    for it in range(1,nt+1):
        print("it=", it)
        uestep(dt)
        ni_hist[:,:,:,it]=bbb.ni
        up_hist[:,:,:,it]=bbb.up
        te_hist[:,:,it]=bbb.te
        ti_hist[:,:,it]=bbb.ti

    
    if (save):
        print("Saving in ", filename)
        
        with h5py.File(filename, 'w') as f:
            f['ni_hist'] = ni_hist
            f['up_hist'] = up_hist
            f['te_hist'] = te_hist
            f['ti_hist'] = ti_hist



            
def restore_time_evol(filename="thist.h5"):
#-restore data from time evolution file
    
    with h5py.File(filename, 'r') as f:
        dataset = f.get('ni_hist')  # get the h5py.Dataset
        ni_hist = dataset[:]  # Copy the array into memory

        dataset = f.get('up_hist')  # get the h5py.Dataset
        up_hist = dataset[:]  # Copy the array into memory

        dataset = f.get('te_hist')  # get the h5py.Dataset
        te_hist = dataset[:]  # Copy the array into memory

        dataset = f.get('ti_hist')  # get the h5py.Dataset
        ti_hist = dataset[:]  # Copy the array into memory

    return ni_hist, up_hist, te_hist, ti_hist
