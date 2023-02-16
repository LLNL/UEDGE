from uedge import *
import uedge.uedge_lists as ul
subrules = [
   ['\(','['],
   ['\)',']'],
   [';','\n'],
   ['^!','#!'],
   ['^ *',''],
   ['^\t*',''],
   [r'\ballocate\b','bbb.allocate()'],
   [r'\bexmain\b','bbb.exmain()'],
   [r'\bexponseed\b','grd.exponseed()'],
   ]
warnrules = []
def raw_string(s):
    s = s.encode('unicode-escape').decode()
    return s

for p in ul.list_packages():
   subrules.append([r'\b'+raw_string('package '+ p)+r'\b','from uedge import '+p])
   po = ul.packagename2object(p)
   for v in ul.list_package_variables(p):
       subrules.append([r'\b'+raw_string(v)+r'\b',p+'.'+v])

       if "Dimension:" in po.listvar(v):
          d = po.listvar(v).split("Dimension:")[1].split("\n")
          if "0:" in d[0]:
             warnrules.append([r'\b'+raw_string(v)+r'\b','base 0, '+d[0]])

       

subrules.append([r'\bbbb.del\b','bbb.delpy'])
