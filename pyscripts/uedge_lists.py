"""
   This module uses some of the Forthon methods to provide routines for
   listing and searching the Uedge compiled packages.

"""

import uedge
import re
packages =  [uedge.com,uedge.aph,uedge.api,uedge.bbb,uedge.flx,uedge.grd,uedge.svr,uedge.wdf,uedge.ncl]


def packagename2object(package):
    for p in packages:
       if p.name() == package: return p
    return None

def list_packages(objects=None):
     """
        Return list of package string names
        or objects if object argument set.
     """
     if objects != None: return packages
     pnames = []
     for p in packages:
       pnames.append(p.name())
     return pnames


def list_package_variables(package,attribute='',vars=None):
     """
        Return list of variable string names from package.
        package - string name of Uedge package.
        attribute='search string' can be either either the 
            group name or an attribute. Search is case
            sensitive and must be exact.
        vars=[varlist] selection limited to varlist
     """
     ret = []
     if type(package) == type(''):
         p = packagename2object(package)
         if p != None: ret.extend(p.varlist(attribute))
     else:
         ret.extend(package.varlist(attribute))
     if vars == None:
         return ret
     else:
         return list(set(ret) & set(vars))

def list_variable(var):
     """
       Print variable information of name passed as a string
       Do not include the package in the variable name.
     """
     for p in packages:
        if var in p.varlist():
            print(p.listvar(var))

def list_variables_glob(s,verbose=False,veryverbose=False,vars=None):
     """
       Print variables where variable name contains string
       Case insensitive
       verbose=True will cause variable comment to print
       veryverbose=True will cause all variable info to print
       vars=[varlist] search limited to varlist
     """
     ret = []
     for p in packages:
        for var in p.varlist():
            if s.upper() in var.upper(): 
               if verbose: print(var+' : '+p.getvardoc(var))
               if veryverbose: print(p.listvar(var))
               ret.append(var)
     if vars == None:
         return ret
     else: 
         return list(set(ret) & set(vars))

def list_variables_apropos(s,verbose=False,veryverbose=False,vars=None):
     """
       Print variables where comment contains string
       Case insensitive
       verbose=True will cause variable comment to print
       veryverbose=True will cause all variable info to print
       vars=[varlist] search limited to varlist
     """
     ret = []
     for p in packages:
        for var in p.varlist():
            if s.upper() in p.getvardoc(var).upper(): 
               if verbose: print(var+' : '+p.getvardoc(var))
               if veryverbose: print(p.listvar(var))
               ret.append(var)
     if vars == None:
         return ret
     else: 
         return list(set(ret) & set(vars))

def list_variables_regex(r,verbose=False,veryverbose=False,vars=None):
     """
       Print variables where comment matches regular expression.
       verbose=True will cause variable comment to print
       veryverbose=True will cause all variable info to print
       vars=[varlist] search limited to varlist
     """
     ret = []
     for p in packages:
        for var in p.varlist():
            if re.search(r,p.getvardoc(var)): 
               if verbose: print(var+' : '+p.getvardoc(var))
               if veryverbose: print(p.listvar(var))
               ret.append(var)
     if vars == None:
         return ret
     else:
         return list(set(ret) & set(vars))

def varlistattr(a):
     """
     Return list of variables with the given attribute. Includes the package
     prefix for use in file save functions.
     """
     ret = []
     for p in packages:
        for var in p.varlist():
           if a in p.getvarattr(var).split():
              ret.append(p.name()+'.'+var)
     return ret
