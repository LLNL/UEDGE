# Holm10 Nov 5 2019, based on rdcontdt.py
# 191121 -  Created hdf5-routines to read and save time-dependent data
#           Writes and reads dictionary with multi-dimensional arrays
#           containing all restore-parameters.

def rundt(  dtreal,nfe_tot=0,savedir='../solutions',dt_tot=0,ii1max=500,ii2max=5,ftol_dt=1e-5,itermx=7,rlx=0.9,n_stor=0,
            tstor=(1e-3,4e-2),incpset=7,dtmfnk3=1e-4):
    ''' Function advancing case time-dependently: increasing time-stepping is the default to attain SS solution
    rdrundt(dtreal,**keys)

    Variables
    dtreal                  The inital time step time

    Keyword parameters:
    nfe_tot[0]              Number of function evaluations
    savedir[savedt]         Directory where hdf5 savefile is written
    dt_tot[0]               Total time accummulated: default option resets time between runs    
    ii1max[500]             Outer loop (dt-changing) iterations
    ii2max[5]               Inner loop (steps at dt) iterations
    ftol_dt[1e-5]           Time-dependent fnrm tolerance 
    itermx[7]               Max. number of linear iterations allowed
    rlx[0.9]                Max. allowed change in variable at each iteration
    n_stor[0]               Number of linearly spaced hdf5 dumps 
    tstor_s[(1e-3,4e-2)]    Tuple with start and stop times for storing snapshots to HDF5
    incpset[7]              Iterations until Jacobian is recomputed
    dtmfnk[1e-4]            dtreal for mfnksol signchange if ismfnkauto=1 (default)
    The above defaults are based on rdinitdt.

    Additional UEDGE parameters used in the function, assuming their default values are:
    bbb.rdtphidtr[1e20]     # Ratio dtphi/dtreal
    bbb.ismfnkauto[1]       # If =1, mfnksol=3 for dtreal<dtmfnk3, otherwise=-3
    bbb.mult_dt[3.4]        # Factor expanding dtreal after each successful inner loop
    bbb.itermxrdc[7]        # Itermx used by the script
    bbb.ftol_min[1e-9]      # Value of fnrm where time advance will stop
    bbb.t_stop[100]         # Value of dt_tot (sec) where calculation will stop
    bbb.dt_max[100]         # Max. time step for dtreal
    bbb.dt_kill[1e-14]      # Min. allowed time step; rdcontdt stops if reached
    bbb.deldt_min[0.04]     # Minimum relative change allowed for model_dt > 0
    bbb.numrevjmax[2]       # Number of dt reductions before Jac recalculated
    bbb.numfwdjmax[1]       # Number of dt increases before Jac recalculated
    bbb.ismmaxuc[1]         # =1 for intern calc mmaxu; =0,set mmaxu & dont chng
    bbb.irev[-1]            # Flag to allow reduced dt advance after cutback
    bbb.initjac[0]          # If=1, calc initial Jac upon reading rdcontdt
    bbb.ipt[1]              # Index of variable; value printed at step
                            # If ipt not reset from unity, ipt=idxte(nx,iysptrx+1)
   

    Additional comments (from rdcontdt):
    This file runs a time-dependent case using dtreal.  First, a converged solution for a (usually small) dtreal is obtained:
    UEDGE must report iterm=1 at the end. Then the control parameters are adjusted. If a mistake is made, to restart this file 
    without a Jacobian evaluation, be sure to reset iterm=1 (=> last step was successful)


    '''
    from uedge import bbb,com
    from uedge.hdf5 import hdf5_save
    from numpy import sqrt,append,array,expand_dims
    from os.path import exists
    from copy import copy
 

    # Store the original values
    dt_tot_o=bbb.dt_tot
    ii1max_o=bbb.ii1max
    ii2max_o=bbb.ii2max
    ftol_dt_o=bbb.ftol_dt 
    itermx_o=bbb.itermx   
    rlx_o=bbb.rlx    
    n_stor_o=bbb.n_stor   
    tstor_s_o=bbb.tstor_s  
    tstor_e_o=bbb.tstor_e 
    incpset_o=bbb.incpset 
    dtmfnk3_o=bbb.dtmfnk3
    icntnunk_o=bbb.icntnunk
    ftol_o=bbb.ftol



    # Set inital time-step to dtreal
    bbb.dtreal=dtreal

    # Check if successful time-step exists (bbb.iterm=1)
    if (bbb.iterm == 1 and bbb.ijactot>1):
        print("Initial successful time-step exists")
        bbb.dtreal = bbb.dtreal*bbb.mult_dt #compensates dtreal divided by mult_dt below
    else:
        print("*---------------------------------------------------------*")
        print("Need to take initial step with Jacobian; trying to do here")
        print("*---------------------------------------------------------*")
        bbb.icntnunk = 0
        bbb.exmain()
        bbb.dtreal = bbb.dtreal*bbb.mult_dt #compensates dtreal divided by mult_dt below

    if (bbb.iterm != 1):
        print("*--------------------------------------------------------------*")
        print("Error: converge an initial time-step first; then retry rdcontdt")
        print("*--------------------------------------------------------------*")
        return
    
    # Set UEDGE variables to the prescribed values
    bbb.dt_tot=dt_tot
    bbb.ii1max=ii1max
    bbb.ii2max=ii2max
    bbb.ftol_dt=ftol_dt 
    bbb.itermx=itermx   
    bbb.rlx=rlx    
    bbb.n_stor=n_stor   
    bbb.tstor_s=tstor[0]  
    bbb.tstor_e=tstor[1] 
    bbb.incpset=incpset 
    bbb.dtmfnk3=dtmfnk3

    # Saved intermediates counter
    i_stor=0

    # Helper variables
    nfe_tot = max(nfe_tot,0)
    deldt_0 = bbb.deldt

    # Empty dictionary for writing
    data=dict() 
    storevar=   [   ['ni',      bbb.ni],
                    ['up',      bbb.up],
                    ['te',      bbb.te],
                    ['ti',      bbb.ti],
                    ['tg',      bbb.tg],
                    ['ng',      bbb.ng],
                    ['phi',     bbb.phi],
                    ['dt_tot',  bbb.dt_tot],
                    ['nfe',     None],
                    ['dtreal',  bbb.dtreal]     ]
    # Linearly spaced time slices for writing 
    dt_stor = (bbb.tstor_e - bbb.tstor_s)/(bbb.n_stor - 1)


    isdtsf_sav = bbb.isdtsfscal


    if(bbb.ipt==1):  # No index requested
        # Check for first variable solved: order is defined as Te,Ti,ni,ng,Tg,phi
        for eq in [bbb.idxte, bbb.idxti, bbb.idxn, bbb.idxg, bbb.idxtg, bbb.idxu]:
            # If multi-species:
            if len(eq.shape)==3:
                # Loop through all species to find first solved
                for index in range(eq.shape[2]):
                    # See if equation is solved
                    if eq[:,:,index].min()!=0:
                        ipt=eq[com.nx-1,com.iysptrx+1,index]
                        break
            # If not, see if equation is solved
            else:
                if eq.min()!=0:
                    ipt=eq[com.nx-1,com.iysptrx+1]
                    break


    bbb.irev = -1         # forces second branch of irev in ii1 loop below
    if (bbb.iterm == 1):  # successful initial run with dtreal
        bbb.dtreal = bbb.dtreal/bbb.mult_dt     # gives same dtreal after irev loop
    else:                 # unsuccessful initial run; reduce dtreal
        bbb.dtreal = bbb.dtreal/(3*bbb.mult_dt) # causes dt=dt/mult_dt after irev loop
       
    if (bbb.initjac == 0): bbb.newgeo=0
    dtreal_sav = bbb.dtreal
    bbb.itermx = bbb.itermxrdc
    bbb.dtreal = bbb.dtreal/bbb.mult_dt	#adjust for mult. to follow; mult_dt in rdinitdt
    bbb.dtphi = bbb.rdtphidtr*bbb.dtreal
    bbb.ylodt = bbb.yl
    bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
    fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq]*bbb.sfscal[0:bbb.neq])**2))
    if (bbb.initjac == 1): fnrm_old=1.e20
    print("initial fnrm ={:.4E}".format(fnrm_old))

    for ii1 in range( 1, bbb.ii1max+1):
        if (bbb.ismfnkauto==1): bbb.mfnksol = 3
        # adjust the time-step
        if (bbb.irev == 0):
            # Only used after a dt reduc. success. completes loop ii2 for fixed dt
            bbb.dtreal = min(3*bbb.dtreal,bbb.t_stop)	#first move forward after reduction
            bbb.dtphi = bbb.rdtphidtr*bbb.dtreal
            if (bbb.ismfnkauto==1 and bbb.dtreal > bbb.dtmfnk3): bbb.mfnksol = -3
            bbb.deldt =  3*bbb.deldt
        else:
            # either increase or decrease dtreal; depends on mult_dt
            bbb.dtreal = min(bbb.mult_dt*bbb.dtreal,bbb.t_stop)
            bbb.dtphi = bbb.rdtphidtr*bbb.dtreal
            if (bbb.ismfnkauto==1 and bbb.dtreal > bbb.dtmfnk3): bbb.mfnksol = -3
            bbb.deldt =  bbb.mult_dt*bbb.deldt
          
        bbb.dtreal = min(bbb.dtreal,bbb.dt_max)
        bbb.dtphi = bbb.rdtphidtr*bbb.dtreal
        if (bbb.ismfnkauto==1 and bbb.dtreal > bbb.dtmfnk3): bbb.mfnksol = -3
        bbb.deldt = min(bbb.deldt,deldt_0)
        bbb.deldt = max(bbb.deldt,bbb.deldt_min)
        nsteps_nk=1
        print('--------------------------------------------------------------------')
        print('--------------------------------------------------------------------')
        print(' ')
        print('*** Number time-step changes = {} New time-step = {:.4E}'.format(ii1, bbb.dtreal))
        print('--------------------------------------------------------------------')

        bbb.itermx = bbb.itermxrdc
        if (ii1>1  or  bbb.initjac==1):	# first time calc Jac if initjac=1
            if (bbb.irev == 1):      # decrease in bbb.dtreal
                if (bbb.numrev < bbb.numrevjmax and \
                    bbb.numrfcum < bbb.numrevjmax+bbb.numfwdjmax): #dont recom bbb.jac
                    bbb.icntnunk = 1	
                    bbb.numrfcum = bbb.numrfcum + 1
                else:                          # force bbb.jac calc, reset numrev
                    bbb.icntnunk = 0
                    bbb.numrev = -1		      # yields api.zero in next statement
                    bbb.numrfcum = 0
                bbb.numrev = bbb.numrev + 1
                bbb.numfwd = 0
            else:  # increase in bbb.dtreal
                if (bbb.numfwd < bbb.numfwdjmax and \
                    bbb.numrfcum < bbb.numrevjmax+bbb.numfwdjmax): 	#dont recomp bbb.jac
                    bbb.icntnunk = 1
                    bbb.numrfcum = bbb.numrfcum + 1
                else:
                    bbb.icntnunk = 0			#recompute jacobian for increase dt
                    bbb.numfwd = -1
                    bbb.numrfcum = 0
                bbb.numfwd = bbb.numfwd + 1
                bbb.numrev = 0			#bbb.restart counter for dt reversals
            bbb.isdtsfscal = isdtsf_sav
            bbb.ftol = max(min(bbb.ftol_dt, 0.01*fnrm_old),bbb.ftol_min)
            bbb.exmain() # take a single step at the present bbb.dtreal
            if (bbb.iterm == 1):
                bbb.dt_tot += bbb.dtreal
                nfe_tot += bbb.nfe[0,0]
                bbb.ylodt = bbb.yl
                bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq-1]*bbb.sfscal[0:bbb.neq-1])**2))
                if (bbb.dt_tot>=0.9999999*bbb.t_stop  or  fnrm_old<bbb.ftol_min):
                    print(' ')
                    print('*****************************************************')
                    print('**  SUCCESS: frnm < bbb.ftol; or dt_tot >= t_stop  **')
                    print('*****************************************************')
                    break

        bbb.icntnunk = 1
        bbb.isdtsfscal = 0
        for ii2 in range( 1, bbb.ii2max+1): #take ii2max steps at the present time-step
            if (bbb.iterm == 1):
                bbb.itermx = bbb.itermxrdc
                bbb.ftol = max(min(bbb.ftol_dt, 0.01*fnrm_old),bbb.ftol_min)
                bbb.exmain()
                if (bbb.iterm == 1):
                    bbb.ylodt = bbb.yl
                    bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                    fnrm_old = sqrt(sum((bbb.yldot[0:bbb.neq-1]*bbb.sfscal[0:bbb.neq-1])**2))
                    print("Total time = {:.4E}; Timestep = {:.4E}".format(bbb.dt_tot,bbb.dtreal))
                    print("variable index ipt = {} bbb.yl[ipt] = {:.4E}".format(ipt,bbb.yl[ipt]))
                    dtreal_sav = bbb.dtreal
                    bbb.dt_tot += bbb.dtreal
                    nfe_tot += bbb.nfe[0,0]
                    if exists(savedir):        
                        hdf5_save('{}/{}_last_ii2.hdf5'.format(savedir,bbb.label[0].decode('UTF-8')))
                    else:
                        print('Folder {} not found, saving output to cwd...'.format(savedir))
                        hdf5_save('{}_last_ii2.hdf5'.format(bbb.label[0].decode('UTF-8')))
                        
                    if (bbb.dt_tot>=0.999999999999*bbb.t_stop  or  fnrm_old<bbb.ftol_min):
                        print(' ')
                        print('*****************************************************')
                        print('**  SUCCESS: frnm < bbb.ftol; or dt_tot >= t_stop  **')
                        print('*****************************************************')
                        break
                    print(" ")
    ##       Store variables if a storage time has been crossed
                    if (bbb.dt_tot >= tstor[0]+dt_stor*i_stor and i_stor<bbb.n_stor):
                        # Check if variables are already present
                        for var in storevar:
                            if var[0]=='nfe':       var[1]=nfe_tot          # Update non-pointer variable
                            if var[0]=='dtreal':    var[1]=copy(bbb.dtreal) # Update variable: unsure why pointer does not update
                            if var[0]=='dt_tot':    var[1]=copy(bbb.dt_tot) # Update variable: unsure why pointer does not update
                            # Check if array initialized
                            if var[0] in data.keys():
                                data[var[0]]=append(data[var[0]],expand_dims(array(copy(var[1])),axis=0),axis=0)
                            else:
                                data[var[0]]=expand_dims(array(copy(var[1])),axis=0)
                                                

                        i_stor = i_stor + 1
       ##          End of storage section
          
        if (bbb.dt_tot>=bbb.t_stop  or  fnrm_old<bbb.ftol_min): break   # need for both loops
        bbb.irev = bbb.irev-1
        if (bbb.iterm != 1):	#print bad eqn, cut dtreal by 3, set irev flag
            itroub()
            
            if (bbb.dtreal < bbb.dt_kill):
                print(' ')
                print('*************************************')
                print('**  FAILURE: time-step < dt_kill   **')
                print('*************************************')
                break
            bbb.irev = 1
            print('*** Converg. fails for bbb.dtreal; reduce time-step by 3, try again')
            print('----------------------------------------------------------------- ')
            bbb.dtreal = bbb.dtreal/(3*bbb.mult_dt)
            bbb.dtphi = bbb.rdtphidtr*bbb.dtreal
            if (bbb.ismfnkauto==1 and bbb.dtreal > bbb.dtmfnk3): bbb.mfnksol = -3
            bbb.deldt =  bbb.deldt/(3*bbb.mult_dt) 
            bbb.iterm = 1


    if bbb.iterm!=1:
        print('Unconverged case dropped out of loop: try again! Terminating...')
        return

    # Save the data to HDF5
    if n_stor>0:
        if exists(savedir):        
            save_dt('{}/dt_{}.hdf5'.format(savedir,bbb.label[0].decode('UTF-8')),data)
        else:
            print('Folder {} not found, saving output to cwd...'.format(savedir))
            save_dt('dt_{}.hdf5'.format(bbb.label[0].decode('UTF-8')),data)
    
    # Restore the original values
    bbb.dt_tot=dt_tot_o
    bbb.ii1max=ii1max_o
    bbb.ii2max=ii2max_o
    bbb.ftol_dt=ftol_dt_o 
    bbb.itermx=itermx_o   
    bbb.rlx=rlx_o    
    bbb.n_stor=n_stor_o   
    bbb.tstor_s=tstor_s_o  
    bbb.tstor_e=tstor_e_o 
    bbb.incpset=incpset_o 
    bbb.dtmfnk3=dtmfnk3_o
    bbb.icntnunk=icntnunk_o
    bbb.ftol=ftol_o
    bbb.dtreal=1e20


def itroub():
    ''' Function that displays information on the problematic equation '''
    from numpy import mod,argmax
    from uedge import bbb
    # Set scaling factor
    scalfac = bbb.sfscal
    if (bbb.svrpkg[0].decode('UTF-8').strip() != "nksol"): scalfac = 1/(bbb.yl + 1.e-30)  # for time-dep calc.

    # Find the fortran index of the troublemaking equation
    itrouble=argmax(abs(bbb.yldot[:bbb.neq]))+1
    print("** Fortran index of trouble making equation is:")
    print(itrouble)

    # Print equation information
    print("** Number of equations solved per cell:")
    print("numvar = {}".format(bbb.numvar))
    print(" ")
    iv_t = mod(itrouble-1,bbb.numvar) + 1 # Use basis indexing for equation number
    print("** Troublemaker equation is:")
    # Verbose troublemaker equation
    if abs(bbb.idxte-itrouble).min()==0:
        print('Electron energy equation: iv_t={}'.format(iv_t))           
    elif abs(bbb.idxti-itrouble).min()==0:
        print('Ion energy equation: iv_t={}'.format(iv_t))   
    elif abs(bbb.idxphi-itrouble).min()==0:
        print('Potential equation: iv_t={}'.format(iv_t))   
    elif abs(bbb.idxu-itrouble).min()==0:
        for species in range(bbb.idxu.shape[2]):
            if abs(bbb.idxu[:,:,species]-itrouble).min()==0:
                print('Ion momentum equation of species {}: iv_t={}'.format(species, iv_t))   
    elif abs(bbb.idxn-itrouble).min()==0:
        for species in range(bbb.idxn.shape[2]):
            if abs(bbb.idxn[:,:,species]-itrouble).min()==0:
                print('Ion density equation of species {}: iv_t={}'.format(species, iv_t))   
    elif abs(bbb.idxg-itrouble).min()==0:
        for species in range(bbb.idxg.shape[2]):
            if abs(bbb.idxg[:,:,species]-itrouble).min()==0:
                print('Gas density equation of species {}: iv_t={}'.format(species, iv_t))   
    elif abs(bbb.idxtg-itrouble).min()==0:
        for species in range(bbb.idxtg.shape[2]):
            if abs(bbb.idxtg[:,:,species]-itrouble).min()==0:
                print('Gas temperature equation of species {}: iv_t={}'.format(species, iv_t))   
    # Display additional information about troublemaker cell
    print(" ")
    print("** Troublemaker cell (ix,iy) is:")
    print(bbb.igyl[itrouble-1,])
    print(" ")
    print("** Timestep for troublemaker equation:")
    print(bbb.dtuse[itrouble-1])
    print(" ")
    print("** yl for troublemaker equation:")
    print(bbb.yl[itrouble-1])
    print(" ")


def save_dt(file,data):
    ''' 
    Save HDF5 file containing time-evolution of restore parameters and time
    Created by holm10 based on meyer8's hdf5.py
    '''
    import h5py
    from time import ctime
    from uedge import bbb

    # Open file for writing
    try:
        hf = h5py.File(file,'w')        # Open file for writing
        hfb = hf.create_group('globals')# Create group for dt data
        hfb.attrs['date'] = ctime()
        hfb.attrs['code'] = 'UEDGE'
        hfb.attrs['ver'] = bbb.uedge_ver
    except ValueError as error:
        print("HDF5 file open failed to {}".format(file))
        print(error)
    except:
        print("HDF5 file open failed to {}".format(file))
        raise

    # Store variables from dictionary data
    for var in ['dt_tot','dtreal','nfe','ni','up','te','ti','tg','ng','phi']:
        try:
            hfb.create_dataset(var, data=data[var])
        except ValueError as error:
            print("{} HDF5 file write failed to {}".format(var,file))
            print(error)
        except:
            print("{} HDF5 file write failed to {}".format(var,file))


def restore_dt(file,ReturnDict=True):
    ''' 
    Restore HDF5 file containing time-evolution of restore parameters and time
    Created by holm10 based on meyer8's hdf5.py
    '''
    import h5py
    from numpy import array    

    data=dict() # Empty dictionary for storage
    
    # Open file for reading
    try:
        hf = h5py.File(file,'r')        # Open file for reading
    except ValueError as error:
        print("HDF5 file open failed to {}".format(file))
        print(error)
        return
    except:
        print("HDF5 file open failed to {}".format(file))
        return

    try:
        dummy = hf['globals']    # Force exception if group not found
        hfb=hf.get('globals')


    except ValueError as error:
        print("HDF5 file could not find group 'globals' in {}".format(file))
        print(error)
        return
    except:
        print("HDF5 file could not find group 'globals' in {}".format(file))
        return
    # Print information on save
    print('Restored time-dependent data for {} case written {} using version {}'.format(hfb.attrs['code'], hfb.attrs['date'],hfb.attrs['ver'][0].decode('UTF-8').strip()[7:-1].replace('_','.')))
    # Loop over all variables
    for var in ['dt_tot','dtreal','nfe','ni','up','te','ti','tg','ng','phi']:
        try:
            data[var] = array(hfb.get(var))
        except ValueError as error:
            print("Couldn't read {} from {}".format(var,file))
            print(error)
        except:
            print("Couldn't read {} from {}".format(var,file))     
 
    if ReturnDict:
        return data
    else:
        return data['dt_tot'],data['dtreal'],data['nfe'],data['ni'],data['up'],data['te'],data['ti'],data['tg'],data['ng'],data['phi']


