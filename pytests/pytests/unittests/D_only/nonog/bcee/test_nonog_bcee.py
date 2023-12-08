# Class for testing UEDGE cases



class TestClass:
    
    def perturb_solution(self, file, perturbation=1e-2):
        from h5py import File
        from numpy.random import uniform
        bbb = file['bbb']
        for var in ['nis', 'ups', 'tes', 'tis', 'ngs', 'tgs', 'phis']:
            bbb[var][...] *= uniform(low=1-perturbation, high=1+perturbation, size = bbb[var].shape)

    def test_reference(self, epsilon=1e-5):
        from h5py import File
        from uedge import bbb, com
        # Supress output
        try:
            com.iprint = 0
        except:
            bbb.iprint = 0
        from os import chdir, getcwd
        from os.path import exists
        from sys import path
        from copy import deepcopy
        from numpy import ndarray, isclose, where, mod
        import importlib.util
        # Then, we can import the UEDGE case setup
        inputspec = importlib.util.spec_from_file_location("input", "input.py")
        inputspec.loader.exec_module(importlib.util.module_from_spec(inputspec))
        bbb.issfon=1
        bbb.ftol=1e20
        # For some reason fnrm still changes after 1 exmain, but not after 2
        # Thus, assert a first exmain here
        bbb.exmain()
        # Helper for grabbing UEDGE data
        def getdata(var):
            from Forthon import package, packageobject
            for pkg in package():
                pkgobj = packageobject(pkg)
                if var in pkgobj.varlist():
                    return deepcopy(pkgobj.__getattribute__(var))

        def setdata(var, value):
            from Forthon import package, packageobject
            for pkg in package():
                pkgobj = packageobject(pkg)
                if var in pkgobj.varlist():
                    pkgobj.__setattr__(var, value)

        def seteqs(file):
            for eq in ['isnion', 'isupon', 'isteon', 'istion', 'isngon', 'istgon', 'isphion']:
                setdata(eq, file[eq][()])

        def restore(file):
            for var in ['nis', 'ups', 'tes', 'tis', 'ngs', 'tgs', 'phis']:
                setdata(var, file['bbb'][var][()])

        def recoverstate(save, refs):
            # Setup UEDGE simulations according to 
            seteqs(refs)
            # Restore save to ensure matched fnrm starting point for comparison
            restore(save)
            # Populate all arrays
            bbb.exmain()

        def matches(save, refs, epsilon):
            from uedge import bbb
            recoverstate(save, refs)
            fnrm = bbb.get_fnrm(bbb.dtreal)
            reffnrm = refs['fnrm'][()]
            return isclose(fnrm, reffnrm, atol=0.0, rtol=epsilon), fnrm, reffnrm


        def print_itroub(refs):
            eqnlist = ['ni', 'up', 'te', 'ti', 'ng', 'tg', 'phi']
            nisp = refs['casesetup/nisp'][()]
            nusp = refs['casesetup/nusp'][()]
            ngsp = refs['casesetup/ngsp'][()]
            eqnarray = [1]
            for eqns in [nisp, nusp, 1, 1, ngsp, ngsp, 1]:
                eqnarray.append(eqns+eqnarray[-1])
            numvar = defref['numvar'][()]
            sclyl = getdata('yldot')*getdata('sfscal')
            rsclyl = defref['yldot'][()]*defref['sfscal'][()]
            iva = abs(rsclyl - sclyl) / (abs(rsclyl) + abs(sclyl) + 1e-20)
            iv = where(iva == max(iva))[0][0]
            (idx, idy) = getdata('igyl')[iv]
            loc_troub_eqn = mod(iv, numvar) + 1
            troub_eqn = next(i for i,v in enumerate(eqnarray) if v >= loc_troub_eqn)
            print('TEST FAILED')
            print('Failure diagnostics for largest RHS residual difference')
            print('   - Variables solved (numvar):'.ljust(40), numvar)
            print('   - Troublemaker equation idex:'.ljust(40), loc_troub_eqn)
            if (eqnarray[troub_eqn+1] - eqnarray[troub_eqn]) > 1:
                print('   - Troublemaker equation[species]:'.ljust(40), '{}[{}]'.format(eqnlist[troub_eqn], loc_troub_eqn-eqnarray[troub_eqn]))
            else:
                print('   - Troublemaker equation:'.ljust(40), eqnlist[troub_eqn])
            print('   - Troublemaker cell index (ix, iy):'.ljust(40), '({}, {})'.format(idx, idy))       

        with File('solution.h5', 'r') as f:
            refs = f['pytests']
            defref = refs['default']
            match, fnrm, reffnrm = matches(f, defref, epsilon)
            if match:
                # Turn output back on
                try:
                    com.iprint = 1
                except:
                    bbb.iprint = 1
                assert True
            else:
                print('Returned fnrm:'.ljust(30), fnrm)
                print('Reference fnrm:'.ljust(30), reffnrm)
                print_itroub(refs)
                print('Failed equation(s):')
                for setupkey in ['ni', 'up', 'te', 'ti', 'ng', 'tg', 'phi']:
                    if setupkey in refs:
                        match, fnrm, reffnrm = matches(f, refs[setupkey], epsilon)
                        if match:
                            continue
                        elif isinstance(refs[setupkey][f'is{setupkey}on'][()], ndarray):
                            species =   refs['casesetup/nisp'][()]*(setupkey == 'ni') +\
                                        refs['casesetup/nusp'][()]*(setupkey == 'up') +\
                                        refs['casesetup/ngsp'][()]*(setupkey in ['ng', 'tg'])
                            # Loop through all variables requested
                            passed, failed = [], []
                            for var in self.variables():
                                if isinstance(getdata(var), ndarray):
                                    close = (not False in isclose(getdata(var), refs[f'{setupkey}'][var][()], atol=0.0, rtol=epsilon))
                                    if close:
                                        passed.append(var)
                                    else:
                                        if len(getdata(var).shape) != 3:
                                            failed.append(var)
                                        else:
                                            for i in range(getdata(var).shape[-1]):
                                                if not False in isclose(getdata(var)[:,:,i], refs[f'{setupkey}'][var][:,:,i], atol=0.0, rtol=epsilon):
                                                    passed.append(f'{var}[{i}]')
                                                else:
                                                    failed.append(f'{var}[{i}]')
                                else:
                                    if isclose(getdata(var), refs[f'{setupkey}'][var][()], atol=0.0, rtol=epsilon):
                                        passed.append(var)
                                    else:
                                        failed.append(var)
                            failindex = []
                            for s in range(species):
                                match, fnrm, reffnrm = matches(f, refs[f'{setupkey}-{s}'], epsilon)
                                if not match:
                                    failindex.append(s)
                            print('   - {} for indices: {}'.format(setupkey, str(failindex)[1:-1]))
                        else:
                            print(f'   - {setupkey}')
                        print('      - Failed variables: {}'.format(str(failed)[1:-1].replace("'",'')))
                # Turn output back on
                try:
                    com.iprint = 1
                except:
                    bbb.iprint = 1
                assert False                    




    def write_data(self):
        from h5py import File
        from uedge import bbb
        from uedge.hdf5 import hdf5_restore
        from copy import deepcopy
        from numpy import ndarray
        # Helper function accessing UEDGE data using variable
        def getdata(var):
            from Forthon import package, packageobject
            for pkg in package():
                pkgobj = packageobject(pkg)
                if var in pkgobj.varlist():
                    return deepcopy(pkgobj.__getattribute__(var))

        def setdata(var, value):
            from Forthon import package, packageobject
            for pkg in package():
                pkgobj = packageobject(pkg)
                if var in pkgobj.varlist():
                    pkgobj.__setattr__(var, value)
        # Write data used to construct tests to save
        casesetup = {}
        with File('solution.h5', 'a') as f:
            self.perturb_solution(f)
            try:
                group = f.create_group('pytests')
            except:
                group = f['pytests']
            try:
                group = group.create_group('casesetup')
            except:
                group = group['casesetup']
            for var in ['nisp', 'nusp', 'ngsp', 'nzsp', 'minu', 'znuclin', 'ziin']:
                casesetup[var] = getdata(var)
                try:
                    group.create_dataset(var, data = getdata(var))
                except:
                    group[var][...] = getdata(var)
        equations = {}
        equations['default'] = {}
        eqlist = ['isnion', 'isupon', 'isteon', 'istion', 'isngon', 'istgon', 'isphion']
        for var in eqlist:
            # Store the original setup
            equations['default'][var] = getdata(var)
            # Turn off all equations in preparation for cycling
            setdata(var, 0)
        for var in eqlist:
            # Get the number of active equations in array
            species =   casesetup['nisp']*(var == 'isnion') +\
                        casesetup['nusp']*(var == 'isupon') +\
                        casesetup['ngsp']*(var in ['isngon', 'istgon'])
            # Get the number of species solved for each equation
            try:
                nspecies = sum([max(x, 0) for x in equations['default'][var][:species]])
            except:
                nspecies = equations['default'][var]
            # Avoid turning on equations the original case was not converged with
            if nspecies > 0:
                # Get a pretty label
                label = var.replace('is', '').replace('on','')
                equations[label] = {}
                for writevar in eqlist:
                    equations[label][writevar] = getdata(writevar)
                # Turn on the current equation only
                equations[label][var] = equations[label][var]**0
                # Multi-species options available
                if isinstance(getdata(var), ndarray):
                    # Loop through all species, turning them on one-by-one
                    for s in range(species):
                        # Turn on only if original input used the setting
                        if equations['default'][var][s] > 0:
                            indexlabel = f'{label}-{s}'
                            equations[indexlabel] = {}
                            for writevar in eqlist:
                                equations[indexlabel][writevar] = getdata(writevar)
                            # Turn on index-wise
                            equations[indexlabel][var][s] = 1

        for key, item in equations.items():
            hdf5_restore('solution.h5')
            for equation, value in item.items():
                setdata(equation, value)
            bbb.exmain()
            fnrm = bbb.get_fnrm(bbb.dtreal)
            with File('solution.h5', 'a') as f:
                try:
                    group = f.create_group('pytests')
                except:
                    group = f['pytests']
                try:
                    group = group.create_group(key)
                except:
                    group = group[key]
                for var in self.variables() + eqlist:
                    try:
                        group.create_dataset(var, data = getdata(var))
                    except:
                        group[var][...] = getdata(var)
                try:
                    group.create_dataset('fnrm', data=fnrm)
                except:
                    group['fnrm'][...] = fnrm
        # Reset equations to default to allow consequtive runs
        for var, value in equations['default'].items():
            setdata(var, value)
        hdf5_restore('solution.h5')
        bbb.exmain()
        # Turn on output again
            

    def variables(self):
        return [
            'yldot',
            'sfscal',
            'numvar',
            'igyl',
            'flox',
            'floy',
            'floxe',
            'floye',
            'floxi',
            'floyi',
            'floxg',
            'floyg',
            'floxge',
            'floyge',
            'conx',
            'cony', 
            'conxe',
            'conye', 
            'conxi',
            'conyi', 
            'conxg',
            'conyg', 
            'conxge',
            'conyge', 
            'visx',
            'visy',
            'hcxe',
            'hcye',
            'hcxij',
            'hcyij',
            'hcxg',
            'hcyg',
            'hcxi',
            'hcxineo',
            'hcyi',
            'kxbohm',
            'kybohm',
            'vybohm',
            'fnix',
            'fniy',
            'fngx',
            'fngy',
            'feix',
            'feiy',
            'feex',
            'feey',
            'fegx',
            'fegy',
            'fmix',
            'fmiy',
            'resco',
            'resng',
            'reseg',
            'resmo',
            'resee',
            'resei',
            'resphi'
        ]
        # TODO: Add consolidated volum. sources
    

