

import sys
try:
     from importlib import reload,import_module
except:
     from importlib import import_module


import builtins


def uexec(mname,returns=globals()):
    if mname in sys.modules:
       _m = reload(sys.modules[mname])
    else:
       _m = import_module(mname)

    # is there an __all__?  if so respect it
    if "__all__" in _m.__dict__:
        names = _m.__dict__["__all__"]
    else:
        # otherwise we import all names that don't begin with _
        names = [x for x in _m.__dict__ if not x.startswith("_")]

    # now drag them in
    for k in names:
        #print k,getattr(_m,k) 
        returns[k] = getattr(_m,k)

