# Automated script for making UEDGE tests
# TODO: READ ALL VARIABLES TO BE COMPARED FROM YAML FILES?



class MakeTests():

    def execute_case(self, testpath, testname):
        from h5py import File
        from uedge import bbb, com
        from uedge.hdf5 import hdf5_restore
        from os import chdir, getcwd, listdir
        from os.path import exists, join
        from sys import path
        from copy import deepcopy
        from numpy import ndarray
        import importlib.util
        # Supress output
        # TODO: settle for one option
        try:
            com.iprint = 0
        except:
            bbb.iprint = 0
        # Switch to test dir to access all files (e.g. grid, mist.dat)
        chdir(testpath)
        # Import based on path in Process - complicated but necessary
        # Set up the UEDGE case
        inputspec = importlib.util.spec_from_file_location("input", "input.py")
        inputspec.loader.exec_module(importlib.util.module_from_spec(inputspec))
        bbb.issfon = 1
        bbb.ftol = 1e20
        # For some reason fnrm still changes after 1 exmain, but not after 2
        # Thus, assert a first exmain here
        bbb.exmain()

        # Import and write test reference data
        testspec = importlib.util.spec_from_file_location("test", testname)
        testmodule = importlib.util.module_from_spec(testspec)
        testspec.loader.exec_module(testmodule)
        testmodule.TestClass().write_data()
        try:
            com.iprint = 1
        except:
            bbb.iprint = 1

    def make_all_tests(self):
        self.make_nonog()
        self.make_carbo()
        self.make_ortho()
        self.make_slab()
        self.make_fulltests()
            
    def make_carbon(self):
        self.make_unittests('pytests/unittests', 'templates', 'nonog_carbon', 
            'D+C', 'uedge_tst_template.py'
        )

    def make_nonog(self):
        self.make_unittests('pytests/unittests', 'templates', 'nonog', 
            'D_only', 'uedge_tst_template.py'
        )

    def make_ortho(self):
        self.make_unittests('pytests/unittests', 'templates', 'ortho', 
            'D_only', 'uedge_tst_template.py'
        )

    def make_slab(self):
        self.make_unittests('pytests/unittests', 'templates', 'slab', 
            'slab_D_only', 'uedge_tst_template.py'
        )

    def make_unittests(self, outdir, commondir, griddir, inputdir, testfile):
        from os import getcwd, chdir, walk, listdir
        from os.path import abspath, join, exists
        from pathlib import Path
        from shutil import copy, copytree, rmtree
        from multiprocessing import Process, Pipe

        def concat_files(file, snippet, header=None):
            """ Adds snippet to file. Args are paths/file names """
            with open(file, 'a') as outfile:
                if header is not None:
                    outfile.write(header)
                with open(snippet, 'r') as infile:
                    for line in infile:
                        outfile.write(line)

        # Concatenate the parts to create the necessary paths
        outpath = join(outdir, inputdir, griddir)
        caseinputpath = join(commondir, inputdir, "inputs")
        datapath = join(commondir, 'grids', griddir)
        # If the test directory exists, remove it to avoid carrying over data
        try:
            rmtree(outpath)
        except:
            pass
        # Create nested dir where to put tests
        Path(outpath).mkdir(parents=True, exist_ok=True)
        # Dictionary to hold test directories
        decks = {}
        # Recursivcely go through the input directories
        for inputroot, _, files in walk(caseinputpath, topdown=False):
            # Ignore all Mac temp files
            files = [i for i in files if i != '.DS_Store']
            # If a folder contains input snippets, store them
            if len(files)>0:
                # Store the list of directories with inputs to dictionary
                decks[inputroot] = files
        # Sort by key (dir name) to allow sorting snippets by naming
        decks = dict(sorted(decks.items(), key=lambda x:x[0]))
        # Open dict with all files to be written out
        files = []
        # Create default path
        default = join(outpath, 'default')
        # Ensure the intiial path exists
        Path(default).mkdir(parents=True, exist_ok=True)
        # Create the master file
        files.append(join(default, 'input.py'))
        # Initialize the grid master
        # Add the header, same for all files in inputdir
        concat_files(join(default, 'input.py'), join(commondir, inputdir, 'header/default.py'))
        # Add grid data to file, same for all files in inputdir
        concat_files(join(default, 'input.py'), join(datapath, 'default.py'))
        for file in listdir(datapath):
            # Files to omit
            if file not in ['default.py']:
                copy(join(datapath, file), default)
        confpath = join(commondir, 'tests', 'conftest.py')
        copy(confpath, default)
        # Loop through all folders containing snippets
        for path, deck in decks.items():
            # Flag for adding default deck to any unmodified input files
            isdefault = False
            # Check whether there is a default file for this setting
            if 'default.py' in deck:
                # Pop the default out for adding after children spawned
                deck.remove('default.py')
                # Set flag for later use
                isdefault = True
            # Loop through any non-default potentially present
            for file in deck:
                # Set up helper path
                newpath = join(outpath, file.split('.')[0])
                # If the file exists, append to it directly
                if exists(newpath):
                    concat_files(join(newpath, "input.py"), join(path, file), f'\n\n# {path}/{file}\n')
                # If the file does not exist, create a new copy from master
                else:
                    # Copy the current maaster files
                    copytree(default, newpath)
                    # Append the newly created file path to the list of files
                    files.append(join(newpath, 'input.py'))
                    # Concatenate the snippet in question to the new test
                    concat_files(files[-1], join(path, file), f'\n\n# {path}/{file}\n')
            # If a default file is present, update all files without a 
            # specified input
            if isdefault:
                # Loop through all files being written
                for f in files:
                    # Check whether other files in deck are to edit existing files
                    # if they are in the deck, they have already been edited
                    # by the corresponding snippet, skip them
                    if f.split('/')[-2] not in [i.split('.')[0] for i in deck]:
                        # If not, add the default snippet to create a complete files
                        concat_files(f, join(path, 'default.py'), f'\n\n# {path}/default.py\n')
        # Make final changes to all cases
        # Create list of parallel UEDGE runs
        cases = []
        # Loop through all test cases created
        for inputfile in files:
            # Add snippet restoring save
            concat_files(inputfile, join(commondir, 'restore/default.py'), '\n\n# {}/restore/default.py\n'.format(commondir))
            # Create helper paths to tests and cases
            casepath = abspath(join(*inputfile.split('/')[:-1]))
            # Ensure the path points to a .py-file
            testpath = join(commondir, 'tests', testfile.replace('.py', '')) + '.py'
            newtestname = 'test_' + '_'.join(inputfile.split('/')[-3:-1]) + '.py'
            # Copy the test-file requested to the test directory
            copy(testpath, join(casepath, newtestname) )
            # Start a subprocess in parallel to ensure order-independent execution
            cases.append(Process(target = self.execute_case, args=(casepath, newtestname), kwargs=()))
            # Start each parallel case
            cases[-1].start()
        # Finally, terminate all parallel cases
        for case in cases:
            case.join()

    def make_fulltests(self):
        self.make_fulltest('pytests/fulltests', 'templates', 'uedge_tst_template.py')

    def make_fulltest(self, outdir, commondir, testfile):
        from os import getcwd, chdir, walk
        from os.path import abspath, join, exists
        from pathlib import Path
        from shutil import copy, copytree, rmtree
        from multiprocessing import Process, Pipe
        # Concatenate the parts to create the necessary paths
        outpath = outdir
        inputdir = join(commondir, 'fulltests')
        # If the test directory exists, remove it to avoid carrying over data
        try:
            rmtree(outpath)
        except:
            pass
        # Create nested dir where to put tests
        Path(outpath).mkdir(parents=True, exist_ok=True)
        # Dictionary to hold test directories
        fulltests = {}
        # Recursivcely go through the input directories
        for inputroot, _, files in walk(inputdir, topdown=False):
            # Ignore all Mac temp files
            files = [i for i in files if i != '.DS_Store']
            # If a folder contains input snippets, store them, omittinc pycache
            if (len(files)>0) and ('__pycache__' not in inputroot):
                # Store the list of directories with case files
                fulltests[inputroot[len(inputdir)+1:]] = files
        # Create list of parallel UEDGE runs
        cases = []
        # Loop through all test cases created
        for testpath, testfiles in fulltests.items():
            # Create helper paths
            outpath = join(outdir, testpath)
            # Ensure the path points to a .py-file
            testfilepath = join(commondir, 'tests', testfile.replace('.py', '')) + '.py'
            newtestname = 'test_' + testpath + '.py'
            # Recursively create the test directory where to store test
            Path(outpath).mkdir(parents=True, exist_ok=True)
            # Loop through all files in the original directory
            for file in testfiles:
                # Copy over all test files to the test directory
                copy(join(inputdir, testpath, file), outpath)
            # Copy over the requested test
            copy(testfilepath, join(outpath, newtestname))
            conftestpath = join(commondir, 'tests', 'conftest.py')
            copy(conftestpath, outpath)
            # Start a subprocess in parallel to ensure order-independent execution
            cases.append(Process(target = self.execute_case, args=(outpath, newtestname), kwargs=()))
            # Start each parallel case
            cases[-1].start()
        # Finally, terminate all parallel cases
        for case in cases:
            case.join()

if __name__=="__main__": 
    MakeTests().make_all_tests()
    
