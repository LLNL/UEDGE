# Holm10 Nov 5 2019, based on rdcontdt.py
# 191121 -  Created hdf5-routines to read and save time-dependent data
#           Writes and reads dictionary with multi-dimensional arrays
#           containing all restore-parameters.
# 230210 - Updated old rundt function to RunData class, modularizing
#          all functionalities. Added a comprehensive diagnostics suite
#          plotting fnrm evolution as function of exmains(), plasma 
#          time, and wall-clock time. Still need to test time-slicing
#   `      procedures
# 230327 - Removed old routine, created wrapper function rundt for 
#          object. Renamed Object to UeRun.

from matplotlib.pyplot import ion
ion()

class UeRun():
    ''' Class containing information on run '''
    def __init__(self, n_stor = False):
        from time import time
        from numpy import array        
        from uedge import bbb, com
        # TODO: Add restore/recover from timeslice
        # TODO: Add plot timeslice directly
        # NOTE: No -> Utilize direct I/O from file instead
        self.tstart = time()
        self.numvar = bbb.numvar
        self.nx = com.nx
        self.ny = com.ny
        self.ixpt1 = com.ixpt1[0]
        self.ixpt2 = com.ixpt2[0]
        self.iysptrx = com.iysptrx
        self.equationkey = array([b'te', b'ti', b'phi', b'up', b'ni', b'ng', b'tg'])

        self.classvars = ['slice_ni', 'slice_ng', 'slice_up', 'slice_te', 
            'slice_ti', 'slice_tg', 'slice_phi', 'slice_dttot', 'time', 
            'fnorm', 'nfe', 'dt_tot', 'dtreal', 'ii1', 'ii2', 'ii1fail', 
            'ii2fail', 'dtrealfail', 'itrouble', 'troubleeq', 'troubleindex',
            'ylfail', 'isteon', 'istion', 'isupon', 'isphion', 'isupgon',
            'isngon', 'istgon', 'ishymol', 'nisp', 'ngsp', 'nhsp', 'nhgsp', 
            'nzsp', 'b0', 'ncore', 'pcoree', 'pcorei']

        # Intiialize all variables to empty lists in class
        for var in self.classvars:
            self.__setattr__(var, [])

    def itroub(self):
        ''' Function that displays information on the problematic equation '''
        from numpy import mod, argmax, where
        from uedge import bbb
        from copy import deepcopy

        self.equations = [bbb.idxte, bbb.idxti, bbb.idxphi, 
            bbb.idxu, bbb.idxn, bbb.idxg, bbb.idxtg]
        equationsdescription = [ 'Electron energy', 'Ion energy', 'Potential',
            'Ion momentum', 'Ion density', 'Gas density', 'Gas temperature']

        # Find the fortran index of the troublemaking equation
        self.neq = bbb.neq
        self.itrouble.append(deepcopy(argmax(abs(bbb.yldot[:self.neq]))+1))
        print("** Fortran index of trouble making equation is:\n{}".format(self.itrouble[-1]))
        # Print equation information
        print("** Number of equations solved per cell:\n    numvar = {}\n".format(self.numvar))

        # TODO: add or subtract one??
        self.troubleeq.append([abs(x-self.itrouble[-1]).min() for x in self.equations].index(0))
#        self.troubleeq.append(ideepcopy(mod(self.itrouble[-1]-1,bbb.numvar) + 1)) # Use basis indexing for equation number

        species = ''
        if self.equations[self.troubleeq[-1]].ndim == 3:
            species = ' of species {}'.format(where(self.equations[\
                self.troubleeq[-1]]-self.itrouble[-1]==0)[-1][0])

        print('** Troublemaker equation is:\n{} equation{}: iv_t={}\n'.format(\
            equationsdescription[self.troubleeq[-1]], species, self.troubleeq[-1]+1))

        # Display additional information about troublemaker cell
        self.troubleindex.append(deepcopy(bbb.igyl[self.itrouble[-1]-1,]))
        self.dtrealfail.append(deepcopy(bbb.dtreal))
        self.ylfail.append(deepcopy(bbb.yl[self.itrouble[-1]-1]))
        print('** Troublemaker cell (ix,iy) is:\n' + \
            '{}\n'.format(self.troubleindex[-1]))
        print('** Timestep for troublemaker equation:\n' + \
            '{:.4e}\n'.format(self.dtrealfail[-1]))
        print('** yl for troublemaker equation:\n' + \
            '{:.4e}\n'.format(self.ylfail[-1]))
    
    def savesuccess(self, ii1, ii2, savedir, savename, fnrm=None):
        from time import time
        from uedge import bbb
        from copy import deepcopy

        self.time.append(time())
        if fnrm is None:
            bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
            self.fnorm.append(deepcopy((sum((bbb.yldot[:bbb.neq]*\
                bbb.sfscal[:bbb.neq])**2))**0.5))
        else:
            self.fnorm.append(fnrm)
        self.nfe.append(deepcopy(bbb.nfe))
        self.dt_tot.append(deepcopy(bbb.dt_tot))
        self.dtreal.append(deepcopy(bbb.dtreal))
        self.ii1.append(ii1)
        self.ii2.append(ii2)
        self.neq = bbb.neq

        self.save_intermediate(savedir, savename)

    def store_timeslice(self):
        from copy import deepcopy
        from uedge import bbb

        self.slice_ni.append(deepcopy(bbb.ni))
        self.slice_ng.append(deepcopy(bbb.ng))
        self.slice_up.append(deepcopy(bbb.up))
        self.slice_te.append(deepcopy(bbb.te))
        self.slice_ti.append(deepcopy(bbb.ti))
        self.slice_tg.append(deepcopy(bbb.tg))
        self.slice_phi.append(deepcopy(bbb.phi))
        self.slice_dttot.append(deepcopy(bbb.dt_tot))

    def save_intermediate(self, savedir, savename):
        from uedge.hdf5 import hdf5_save
        from uedge import bbb, com
        from h5py import File

        for var in [ 'isteon', 'istion', 'isupon', 'isphion', 'isupgon',
            'isngon', 'istgon', 'ishymol']:
            self.__setattr__(var, bbb.__getattribute__(var))
        for var in [ 'nisp', 'ngsp', 'nhsp', 'nhgsp', 'nzsp']:
            self.__setattr__(var, com.__getattribute__(var))

        hdf5_save('{}/{}_last_ii2.hdf5'.format(savedir,savename))
        try:
            hdf5_save('{}/{}_last_ii2.hdf5'.format(savedir,savename))
        except:
            print('Folder {} not found, saving output to cwd...'.format(savedir))
            hdf5_save('{}_last_ii2.hdf5'.format(savename))

        try:
            file = File('{}/{}_last_ii2.hdf5'.format(savedir, savename, 'r+'))
        except:
            file = File('{}_last_ii2.hdf5'.format(savename), 'r+')

        file.require_group('convergence')
        group = file['convergence']
        group.create_dataset('t_start', data=self.tstart)
        group.create_dataset('numvar', data=self.numvar)
        group.create_dataset('neq', data=self.neq)
        group.create_dataset('nx', data=self.nx)
        group.create_dataset('ny', data=self.ny)
        group.create_dataset('ixpt1', data=self.ixpt1)
        group.create_dataset('ixpt2', data=self.ixpt2)
        group.create_dataset('iysptrx', data=self.iysptrx)
        group.create_dataset('equationkey', data=self.equationkey)
        group.create_dataset('itermx', data=self.itermx)
        group.create_dataset('incpset', data=self.incpset)
        group.create_dataset('ii1max', data=self.ii1max)
        group.create_dataset('ii2max', data=self.ii1max)
        group.create_dataset('numrevjmax', data=self.numrevjmax)
        group.create_dataset('numfwdjmax', data=self.numfwdjmax)
        group.create_dataset('numtotjmax', data=self.numtotjmax)
        group.create_dataset('rdtphidtr', data=self.rdtphidtr)
        group.create_dataset('deldt_min', data=self.deldt_min)
        group.create_dataset('rlx', data=self.rlx)

        for var in self.classvars:
            group.create_dataset(var, data=self.__getattribute__(var))

        file.close()
    
    def convergenceanalysis(savefname, savedir='../solutions', fig=None,
        xaxis = 'exmain', logx = False, color='k', label=None,
        ylim = (None, None)):
        from h5py import File
        from matplotlib.pyplot import subplots
        from datetime import timedelta
        from matplotlib.ticker import FuncFormatter
        from numpy import cumsum, ones
        if fig is None:
            f, ax = subplots(1, 3, figsize=(15, 5))
        else:
            ax = fig.get_axes()
            if len(ax) < 3:
                print('Three subplots required for plots! Aborting...')
                return
            f = fig
        try:
            file = File('{}/{}'.format(savedir, savefname), 'r')
        except:
            print('File {}/{} not found. Aborting!'.format(savedir, 
                savefname))
            return
        data = file['convergence']
        try:
            data = file['convergence']
        except:
            print('Convergence data not found in {}/{}. Aborting!'.format(\
                savedir, savefname))
            return
        
        if xaxis == 'exmain':
            xlabel = 'Exmain calls'
            xones = ones(data['ii2'][()].shape)
            x = cumsum(xones)
        elif xaxis == 'nfe':
            xlabel = 'nfe internal iterations'
            x = cumsum(data['nfe'][()][:, 0, 0])
        elif xaxis == 'time':
            xlabel = 'Total wall-clock time [HH:MM]'
            x = [timedelta(t - data['t_start'][()]) for t in data['time'][()]]
            x = data['time'][()] - data['t_start'][()]
            
        if logx is True:
            ax[0].loglog(x, data['fnorm'][()], '-', color=color, label=label)
            ax[1].loglog(data['dt_tot'][()], data['fnorm'][()], '-', color=color, label=label)
            ax[2].loglog(x, data['dtreal'][()], '-', color=color, label=label)
        else:
            ax[0].semilogy(x, data['fnorm'][()], '-', color=color, label=label)
            ax[1].semilogy(data['dt_tot'][()], data['fnorm'][()], '-', color=color, label=label)
            ax[2].semilogy(x, data['dtreal'][()], '-', color=color, label=label)

        ax[0].set_xlabel(xlabel)
        ax[1].set_xlabel('Accumulated plasma simualtion time [s]')
        ax[2].set_xlabel(xlabel)
        ax[1].set_title('Total exmain evaluations: {}'.format\
            (len(data['dtreal'][()])))
        ax[0].set_ylabel('Initial fnorm')
        ax[1].set_ylabel('Initial fnorm')
        ax[2].set_ylabel('Time-step (dtreal) [s]')
        ax[0].set_ylim(ylim)
        ax[1].set_ylim(ylim)
        if xaxis == 'time':
            ax[0].xaxis.set_major_formatter(FuncFormatter(lambda t, pos :\
                    str(timedelta(seconds=t))[:-3]))
            ax[2].xaxis.set_major_formatter(FuncFormatter(lambda t, pos :\
                    str(timedelta(seconds=t))[:-3]))
        if label is not None:
            ax[0].legend()

        return f


    def failureanalysis(savefname, savedir='../solutions', equation=None):
        from h5py import File
        from matplotlib.pyplot import subplots
        from numpy import histogram, zeros
        from matplotlib.collections import PolyCollection
        
        f, ax = subplots(2,1, figsize=(10, 7))
        try:
            file = File('{}/{}'.format(savedir, savefname), 'r')
        except:
            print('File {}/{} not found. Aborting!'.format(savedir, 
                savefname))
        data = file['convergence']
        try:
            data = file['convergence']
        except:
            print('Convergence data not found in {}/{}. Aborting!'.format(\
                savedir, savefname))
            return

        if equation is not None:
            iequation = [x.decode('UTF-8') for x in data['equationkey']].index(equation)
        # Bin the equation errors
        counts, bins = histogram(data['troubleeq'][()], bins=7, range=(-0.5,6.5))
        ax[0].hist(bins[:-1], bins, weights=counts)
        ax[0].set_xticks(range(7))
        ax[0].set_xticklabels([x.decode('UTF-8') for x in data['equationkey'][()]])
        ax[0].grid(linestyle=':', linewidth=0.5, axis='y')


        # Visualize error locations
        nx = data['nx'][()]        
        ny = data['ny'][()]        
        ixpt1 = data['ixpt1'][()]
        ixpt2 = data['ixpt2'][()]
        iysptrx = data['iysptrx'][()]
        frequency = zeros((nx, ny))

        cells = []
        for i in range(nx):
            for j in range(ny):
                cells.append([[i+0.5, j+0.5], [i+1.5, j+0.5], 
                    [i+1.5, j+1.5], [i+0.5, j+1.5]])
        polys = PolyCollection(cells, edgecolors='k', linewidth=0.5, linestyle=':')
            


        for i in range(len(data['itrouble'])):
            coord = data['troubleindex'][()][i]
            if equation is None:
                frequency[coord[0]-1, coord[1]-1] += 1
            elif iequation == data['troubleeq'][()][i]:
                frequency[coord[0]-1, coord[1]-1] += 1

        polys.set_cmap('Reds')
        polys.set_array(frequency.reshape((nx*ny,)))

        cbar = f.colorbar(polys, ax=ax[1])
        cbar.ax.set_ylabel('N trouble'+' for {}'.format(equation)*\
            (equation is not None), va='bottom', labelpad=20)

        ax[1].set_xlabel('Poloidal index')
        ax[1].set_ylabel('Radial index')
        ax[1].add_collection(polys)
        ax[1].plot([0.5, nx+0.5],[iysptrx+0.5, iysptrx+0.5], 'k-', linewidth=2)
        ax[1].plot([ixpt1+0.5, ixpt1+0.5], [0.5, iysptrx+0.5], 'k-', linewidth=2)
        ax[1].plot([ixpt2+0.5, ixpt2+0.5], [0.5, iysptrx+0.5], 'k-', linewidth=2)

        file.close()
        return f



    def converge(self, dtreal=2e-9, ii1max=5000, ii2max=5, itermx=7, ftol=1e-5,
        dt_kill=1e-14, t_stop=100, dt_max=100, ftol_min = 1e-9, incpset=7,
        n_stor=0, storedist='lin', numrevjmax=2, numfwdjmax=1, numtotjmax=0, 
        tstor=(1e-3, 4e-2), ismfnkauto=True, dtmfnk3=5e-4, mult_dt=3.4, 
        reset=True, initjac=False, rdtphidtr=1e20, deldt_min=0.04, rlx=0.9,
        tsnapshot=None, savedir='../solutions', ii2increase=1.5):

        ''' Converges the case by increasing dt 
        dtreal : float [1e-9]
            Original time-step size
        ii1max : int [500]
            Outer loop iterations, i.e. time-step changes
        ii2max : int [5]
            Inner loop iterations, i.e. time-steps per time-step change
        dt_kill : float [1e-14]
            Time-step limit for aborting simulation
        itermx : int [7]
            Maximum iterations per time-step used internally in routine
        ftol : float [1e-5]
            Internal fnrm tolerance for time-steps
        incpset : int [7]
            
        savedir : str ['../solutions']

        numtotjmax : int [None]

        ii2increase : float [1.5]
            
        ftol_min : float [1e-9]
            Value of fnrm where time-advance will stop
        t_stop : float [100.]
            Maximum total accumulated plasma-time before stopping if 
            fnorm has not decreased below ftol_min
        dt_max : float [100.]
            Maximum allowable time-step size
        numrevjmax : int [2]
            Number of time-step reducitons before Jacobian is recomputed
        numfwdjmax : int [2]
            Number of time-step increases before Jacobian is recomputed
        n_stor : int [0]
            Number of time-slices to be stored in interval tstor
        tstor : tuple of floats [(1e-3, 4e-2)]
            Time-interval in which time-slices are stored (lower, upper)
        storedist : str ['lin']
            Distribution of time-slices in tstor. Options are 'lin' and 
            'log' for linear and logarithmic distributions, respectively
        reset : bool [True]
            Switch whether to reset the total time etc of the case
        initjac : bool [False]
            Switch to re-evaluate Jacobian on first convegnec time-step
            or not
        ismfnkauto : bool [True]
            If True, sets mfnksol=3 for time-steps smaller that dtmfnk3,
            mfnksol=-3 for larger time-step sizes
        dtmfnk3 : float [5e-4]
            Time-step size for which ismfnkauto controls mfnksol if
            ismfnkauto is True
        mult_dt : float [3.4]
            Time-step size increase factor after successful inner loop
        rdtphidtr : float [1e20]
            Ratio of potential-equation time-step to plasma equaiton time-step 
            size: dtphi/dtreal
        deldt_min : float [0.04]
            Minimum relative change allowed for model_dt>0
        rlx : float [0.9]
            Maximum change in variable at each internal linear iteration
        tsnapshot : list [None]
            If None, uses linear/logarithmic interpolation according to storedist
            in the interval tstor. Snapshot times can be defined in a list and 
            supplied. Then, the tnsaphsot list defines the time-slices
         
        '''
        from numpy import linspace, logspace, log10, append
        from copy import deepcopy
        from uedge import bbb
        from os.path import exists


        # Check if requested save-directory exists: if not, write to cwd
        if not exists(savedir):
            print('Requested save-path {} not found, writing to cwd!'.format(\
                savedir))
            savedir = '.'        


        self.orig = {}
        self.orig['itermx'] = deepcopy(bbb.itermx)
        self.orig['dtreal'] = deepcopy(bbb.dtreal)
        self.orig['icntnunk'] = deepcopy(bbb.icntnunk)
        self.orig['ftol'] = deepcopy(bbb.ftol)
        self.orig['mfnksol'] = deepcopy(bbb.mfnksol)
        self.orig['rlx'] = deepcopy(bbb.rlx)
        self.orig['deldt'] = deepcopy(bbb.deldt)
        self.orig['isdtsfscal'] = deepcopy(bbb.isdtsfscal)
        self.orig['incpset'] = deepcopy(bbb.incpset)

        if numtotjmax == 0:
            numtotjmax = numrevjmax + numfwdjmax


        self.itermx = itermx
        self.incpset = incpset
        self.ii1max = ii1max
        self.ii2max = ii2max
        self.numrevjmax = numrevjmax
        self.numfwdjmax = numfwdjmax
        self.numtotjmax = numtotjmax
        self.rdtphidtr = rdtphidtr
        self.deldt_min = deldt_min
        self.rlx = rlx

# TODO: Add variable to control reduciton factor?
# TODO: Should dtreal = min(x, t_stop) actually be t_stop or dt_max?

        def restorevalues(self):
            ''' Restores the original UEDGE values '''
            for key, value in self.orig.items():
                bbb.__setattr__(key, value)

        def message(string, separator='-', pad='', seppad = '', 
            nseparator=1):
            ''' Prints formatted message to stdout '''
            # TODO: add formatting for len>75 strings
            message = pad.strip() + ' ' + string.strip() + ' ' + pad.strip()
            for i in range(nseparator):
                print(seppad + separator*(len(message)-2*len(seppad)) + seppad)
            print(message)
            print(seppad + separator*(len(message)-2*len(seppad)) + seppad)

        def scale_timestep(scaling):
            ''' Increases/decreases time-step ''' 
            bbb.dtreal *= scaling

        def exmain_isaborted(self):
            ''' Checks if abort is requested '''
            from uedge import bbb
            bbb.exmain()
            # Abort flag set, abort case
            if bbb.exmain_aborted == 1: 
                # Reset flag
                bbb.exmain_aborted == 0 
                # Restore parameters modified by script
                restorevalues(self)
                return True

        def issuccess(self, t_stop, ftol_min):
            ''' Checks if case is converged '''
            from datetime import timedelta
            from time import time
            if (bbb.iterm == 1):
                bbb.ylodt = bbb.yl
                bbb.dt_tot += bbb.dtreal
                bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                self.fnrm_old = sum((bbb.yldot[:bbb.neq-1]*bbb.sfscal[:bbb.neq-1])**2)**0.5
                self.savesuccess(ii1, ii2, savedir, bbb.label[0].strip(\
                    ).decode('UTF-8'), self.fnrm_old)
                if (bbb.dt_tot>=t_stop  or  self.fnrm_old<ftol_min):
                    print('')
                    message('SUCCESS: ' + 'fnrm < bbb.ftol'*(self.fnrm_old<ftol_min) +
                        'dt_tot >= t_stop'*(bbb.dt_tot >= t_stop), pad='**', 
                        separator='*')
                    print('Total runtime: {}'.format(timedelta(
                        seconds=round(time()-self.tstart))))
                    restorevalues(self)
                    return True
    
        def isfail(dt_kill):
            ''' Checks whether to abandon case '''
            if (bbb.dtreal < dt_kill):
                message('FAILURE: time-step < dt_kill', pad='**', 
                separator='*')
                restorevalues(self)
                return True

        def setmfnksol(ismfnkauto, dtmfnk3):
            ''' Sets mfnksol according to setup '''
            if ismfnkauto is True:
                bbb.mfnksol = 3*(-1)**(bbb.dtreal > dtmfnk3)

        def calc_fnrm():
            ''' Calculates the initial fnrm '''
            from uedge import bbb
            bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
            return sum((bbb.yldot[:bbb.neq-1]*bbb.sfscal[:bbb.neq-1])**2)**0.5
 
        ''' TIME-SLICING SETUP '''
        if tsnapshot is None:
            if storedist == 'lin':
                # Linearly spaced time slices for writing 
                dt_stor = linspace(tstor[0], tstor[1], n_stor)
            elif storedist == 'log':
                # Logarithmically spaced time-slices
                dt_stor = logspace(log10(tstor[0]), log10(tstor[1]), n_stor)
        else:
            dt_stor = tsnapshot
        # Add end-point to avoid tripping on empty arrays
        dt_stor = append(dt_stor, 1e20)

        if reset is True:
            bbb.dt_tot = 0               

        ''' TIME-STEP INITIALIZATION '''
        bbb.rlx = rlx
        bbb.dtreal = dtreal
        bbb.ftol = ftol
        if (bbb.iterm == 1) and (bbb.ijactot > 0):
            message('Initial successful time-step exists', separator='')
        else:
            message('Need to take initial step with Jacobian; ' + \
                'trying to do here', seppad='*')
            # Ensure time-step is taken
            bbb.icntnunk = 0
            # Take timestep and see if abort requested
            if exmain_isaborted(self):
                return
            # Increase time
            # Verify time-step was successful
            if (bbb.iterm != 1):
                restorevalues(self)
                message('Error: converge an initial time-step first; then ' + \
                    're-execute command', seppad='*')
                return
        bbb.incpset = incpset
        bbb.itermx = itermx
        deldt_0 = deepcopy(bbb.deldt)
        isdtsf_sav = deepcopy(bbb.isdtsfscal)
# TODO: Replace with some more useful information?
#        if (bbb.ipt==1 and bbb.isteon==1): 	# set ipt to te(nx,iysptrx+1) if no user value
#           ipt = bbb.idxte[nx-1,com.iysptrx]  #note: ipt is local, bbb.ipt global
        bbb.dtphi = rdtphidtr*bbb.dtreal
        svrpkg=bbb.svrpkg.tostring().strip()
        bbb.ylodt = bbb.yl
        self.fnrm_old = calc_fnrm()
        if initjac is True: 
            self.fnrm_old = 1e20
        else:
            bbb.newgeo=0
        # Intialize counters
        irev = -1
        numfwd = 0
        numrev = 0
        numrfcum = 0
        # Compensate for first time-step before entering loop
        scale_timestep(1/(3*(irev == 0) + mult_dt*(irev != 0)))
        ''' OUTER LOOP - MODIFY TIME-STEP SIZE'''
        # TODO: Add logic to always go back to last successful ii2 to 
        # precondition the Jacobian, to avoid downwards cascades?
        # NOTE: experomental functionality
        successivesuccesses = 0
        for ii1 in range(ii1max):
            setmfnksol(ismfnkauto, dtmfnk3)
            # adjust the time-step
            # dtmult=3 only used after a dt reduc. success. completes loop ii2 for fixed dt
            # either increase or decrease dtreal; depends on mult_dt
            scale_timestep(3*(irev == 0) + mult_dt*(irev != 0))
            bbb.dtreal = min([bbb.dtreal, dt_max]) 
            bbb.dtphi = rdtphidtr*bbb.dtreal
            bbb.deldt =  min([bbb.deldt, deldt_0, deldt_min])
            message('Number of time-step changes = ''{} New time-step: {:.2E}\n'\
                .format((ii1+1), bbb.dtreal), pad='***', nseparator=1)

            # Enter for every loop except first, unless intijac == True
            if ii1 > -int(initjac): 
                # Main time-stepping switch: controls increase/decrease in 
                # dtreal and Jacobian preconditioning
                if (irev == 1):      # decrease in bbb.dtreal
                    if (numrev < numrevjmax and \
                        numrfcum < numtotjmax): #dont recom bbb.jac
                        bbb.icntnunk = 1	
                        numrfcum += 1
                    else:                          # force bbb.jac calc, reset numrev
                        bbb.icntnunk = 0
                        numrev = -1		      # yields api.zero in next statement
                        numrfcum = 0
                    numrev += 1
                    numfwd = 0
                else:  # increase in bbb.dtreal
                    if (numfwd < numfwdjmax and \
                        numrfcum < numtotjmax): 	#dont recomp bbb.jac
                        bbb.icntnunk = 1
                        numrfcum += 1
                    else:
                        bbb.icntnunk = 0			#recompute jacobian for increase dt
                        numfwd = -1
                        numrfcum = 0
                    numfwd += 1
                    numrev = 0			#bbb.restart counter for dt reversals
                bbb.isdtsfscal = isdtsf_sav
                # Dynamically decrease ftol as the initial ftol decreases
                bbb.ftol = max(min(ftol, 0.01*self.fnrm_old),ftol_min)
                # Take timestep and see if abort requested
                if exmain_isaborted(self):
                    return
                if issuccess(self, t_stop, ftol_min):
                    return
            bbb.icntnunk = 2
            bbb.isdtsfscal = 0
            # NOTE: experomental functionality
            bbb.ii2max = ii2max + round(ii2increase*successivesuccesses)
            # Take ii2max time-steps at current time-step size while 
            # time-steps converge: if not, drop through
            for ii2 in range(bbb.ii2max): 
                if (bbb.iterm == 1):
                    bbb.ftol = max(min(ftol, 0.01*self.fnrm_old),ftol_min)
                    # Take timestep and see if abort requested
                    message("Inner iteration #{}".format(ii2+1), nseparator=0, 
                        separator='')
                    if exmain_isaborted(self):
                        return
                    if issuccess(self, t_stop, ftol_min):
                        return
                    message("Total time = {:.4E}; Timestep = {:.4E}".format(\
                        bbb.dt_tot-bbb.dtreal,bbb.dtreal), nseparator=0, 
                        separator='')
# TODO: replace with more useful information
#                       print("variable index ipt = ",ipt, " bbb.yl[ipt] = ",bbb.yl[ipt])
                    # Store variable if threshold has been passed
                    if (bbb.dt_tot >= dt_stor[0]):
                        # Remove storing time-points smaller than current 
                        # simulation time
                        while bbb.dt_tot >= dt_stor[0]:
                            dt_stor = dt_stor[1:]
                        self.store_timeslice()
            irev -= 1
            # Output and store troublemaker info
            # NOTE: experomental functionality
            successivesuccesses += 1
            if (bbb.iterm != 1):	
                # NOTE: experomental functionality
                successivesuccesses = 0
                self.itroub()
                ''' ISFAIL '''
                if isfail(dt_kill):
                    self.save_intermediate(savedir, bbb.label[0].strip().decode('UTF-8'))
                    break
                irev = 1
                message('Converg. fails for bbb.dtreal; reduce time-step by ' + \
                    '3, try again', pad = '***', nseparator=0)
                scale_timestep(1/(3*mult_dt))
                bbb.dtphi = rdtphidtr*bbb.dtreal
                bbb.deldt *=  1/(3*mult_dt) 
                setmfnksol(ismfnkauto, dtmfnk3)
#                bbb.iterm = 1
#        bbb.iterm = -1 # Ensure subsequent repetitions work as intended


def rundt(**kwargs):
    runcase=UeRun()
    runcase.converge(**kwargs)


