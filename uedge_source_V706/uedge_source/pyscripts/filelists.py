# Methods to create a list of files in a directory, and a sub-list containing sspecified suffixes.

import os
import string
def filelist(path):
   return os.listdir(path)

def filesublist(path,suffix):
   sublist = []
   biglist = filelist(path)
   for filename in biglist:
       filesplit = string.split(filename,".")
       if filesplit[-1] == suffix:
           sublist.append(filename)
   return sublist


