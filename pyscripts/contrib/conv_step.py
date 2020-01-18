# Holm10 Sep 10 2019, created from scratch


def conv_ncore_step(d, name, t_stop=100,ii1max=100,nstop=1.5e20,dtreal=1e-9):
    """ Simple function to incrementally change the density of a converged solution """
    from uedge import bbb
    from uedge.rundt import rundt
    from uedge.hdf5 import hdf5_save
    
    isbcwdt=bbb.isbcwdt

    # Check if we are increasing or decreasing!
    if d>0:
        increasing=True
    else:
        increasing=False



    # Setup dt-run
    bbb.t_stop=t_stop   # Set stop time
    bbb.ii1max=ii1max   # Set max outer loop iterations 

    while True:    
        bbb.dtreal=dtreal # Small step size to ensure convergence
        bbb.ncore[0]+=d # Increase step size
        print('===================================')
        print('Solving for ncore[0]={:.2E}'.format(bbb.ncore[0]))
        print('===================================')
        '''
        bbb.exmain() # Check convergence at small dt
        # Check that we can get started
        if bbb.iterm!=1:    # The case did not converge
            bbb.isbcwdt=1       # Relax BC:s and try again
            bbb.exmain()            # Check convergence
            if bbb.iterm!=1:        # Case will not converge
                return "Case does not converge at dtreal=1e-9 w/ isbcwdt=1. Aborting..."

        # We have an initially converging case. Start dt-rund
        bbb.dt_tot=0        # Reset time
        if bbb.isbcwdt==1:  # We have relaxed BC:s - do a small step first
            # Advance to a micro-second
            bbb.t_stop=1e-5
            rundt(bbb.dtreal)
            # Set BCs and run to SS
            bbb.isbcwdt=0
            bbb.t_stop=t_stop
            rundt(bbb.dtreal)
        else:
            rundt(bbb.dtreal)
        '''
        rundt(bbb.dtreal)
        # If run with BC relaxation, ensure convergence without relaxation
        if bbb.isbcwdt==1:
            bbb.isbcwdt=0
            rundt(1e-6)

        
        # We should now have a steady-state solution: ensure this!
        bbb.dtreal=1e20
        bbb.itermx=30
        bbb.icntnunk=0
        bbb.ftol=1e-8
        bbb.exmain()
        # Check if SS or not
        if bbb.iterm==1:   # SS solution
            # Save to solutions with appropriate name
            hdf5_save("../solutions/{}_{:.3E}_ss.hdf5".format(name,bbb.ncore[0]))
        else:
            hdf5_save("../solutions/{}_{:.2E}_failed.hdf5".format(name,bbb.ncore[0]))
            print("Ramp did not converge for variable: {:.2E}. Aborting...".format(bbb.ncore[0]))
            break
        if increasing:
            if bbb.ncore[0]>nstop:
                break
        else:
            if bbb.ncore[0]<nstop:
                break
        bbb.isbcwdt=isbcwdt






        

