# The UEDGE test module

Created by A. Holm Dec 05, 2023 based on previous pytests by W. H. Meyer.

## Executing the tests


This python test area is setup to be run with the "pytest" module. You may 
need to load this into your Anaconda Python environment with:
```
    pip install pytest pytest-xdist pytest-isolate
```
Then to run all the tests type:
```
    pytest --isolate --tb=native
```
To run the test in a specific folder:
```
    pytest --isolate --tb=native <path to folder of tests>
```
To run just a specific test you can either do this while in this directory:
```
    pytest --isolate --tb=native <path to test folder>
```
or
```
    cd <path to test folder>
    pytest --isolate --tb=native
```
Repeat for any desired test.

## Description of the standard tests

The way Forthon interfaces UEDGE with Python requires the tests to be run in separate
processes in order for the test results to be order-independent. The UEDGE setup necessitates
the tests to be performed in separate folders for UEDGE to find the appropriate files required
to set up the cases. These facts dictate the design of the nested test folders.

The standard test is implemented to detect changes to the source code, detected by 
looking for changes in the unconverged fnrm of each case. The solution file is randomly perturbed 
by +/-1% for every cell, in order to detect changes in the source code influencing fnrm to 
avoid false positives caused by numerical noise due to different compilers/machines. The test 
tolerance for a successful test is a 1e-3 of the current fnrm relative to the reference 
fnrm. A warning is issued when the difference exceeds a 1e-5 tolerance. The test uses
numpy.isclose(comparison, reference, atol=0.0, ftol=1e-5) to assess whether the test
passes or not.


## Creating/updating the test module

The UEDGE test module is self-assembled by linearly combining a large number of Python
code snippets. This means the tests can be automatically created any time changes to
the source code changing the results of the code are verified and need to be consolidated
into the codebase for future tests to pass.

### Recreating the test module

To make directory of tests based on an updated uedge version, execute
```    
    python MakeTests.py
```
<MORE DOCUMENTATION TO FOLLOW>


