# A Python translator that will convert files following a set
# of file-suffix-specific substitution rules
# Written by Ron Cohen, February 2007.  Last revision March 15, 2007.

# import utilities to create filtered list of files 
from filelists import *
import os
import types
import filecmp

"""
create a list of classes to enstantiate which will have conversion rules
for different language files.

Usage:

  Create a file "localrules.py" that defines globalsubrules to apply to
  all file types and any of cppsubrules, pythonsubrules,fortransubrules,
  mpplsubrules,F90subrules that you want to apply to specific file
  types.  These rules are of the form:
  globalsubrules = [("in1","out1"),("in2","out2"),
    ...] where string in1 will be converted to string out1, etc.  It is
    an ordered list which will be executed first to last.
    Another option is that an entry in globalsubrules can be the name
    of a function of one argument (the string being processed).
  Note currently, if the ordered list returns "None" (the python None,
    not the string "None")for a line,
    no further proceessing of that line is done and no line is written.
    This can be used to re-order lines, since a function entry used in
    globalsubrules can be used to return a list of lines, or "None". But
    no further processing is done of the returned list of lines, so any
    such re-ordering must be done as the last of the rules..
  Note in additon to duples, one can also include in the subrules methods
    which take as its sole argument the string the method operates on;
    use this for more complicated operations than simple replacements
  Also create a local script in the directory where you want to
  convert files. Begin that script with "from convert import *",
  making sure that convert.py and localrules.py are in your path.
  Warning: don't end the local script with the suffix .py or the convertor
  will appy your rules to `your local script (assuming you use the default
  suffix as the suffxin for class Py). 

  Further notes

  1.If desired, the default suffixes searched to designate files of
    a specific type can be edited in your local script.  to do so
    import convert, edit the suffix for the sub-class (for example,
    convert.py.suffix = "newsuffix"), then type "from convert import *".

  2.Create sub-classes that inherit from class generic, for any additional
    file suffix that will be converted.  In each subclass, set suffixin
    to a string that gives the suffix of files to be converted, e.g. ".cc"
    or ".py".  If the processed files are to have a different suffix,
    set suffixout to the desired suffix.

  3.Append processall.classlist, a list of sub-class names that will
    be processed to do the conversion, to add the names of any classes
    added in step 2.

  Then, to proces files:
  1.Execute processall(indir,outdir) to process all the specified file
    types.  indir and outdir are optional arguments to specify paths
    for the directory to be processed and the directory where processed
    files will be written.  By default indir = ".", the current directory,
    and outdir = "./converted".
  2.You can also create instances of specific classes created in step 2, and
    just process those files.  For example a = Py(indir,outdir).  See
    documentation for class generic to see available methods.
    WARNING: If you've run the script before and already created the outdir,
    the script will overwrite files in outdir with the same name
  3.The script will only process files if the source file is newer than
    the target (or the target doesn't exist).
  4.You can also use processdirs(dirlist) to run processall in a list of
    directory names specified in dirlist (names in quotes).
"""

# Global substitution dictionary, globalsubdict = {"instring1:outstring1",...}
globalsubdict = {}

# subrules of form [("in1","out1"),("in2","out2"), ...]
globalsubrules = []

# substitution rules to be appended to globalsubrules for each language
cppsubrules = []
# e.g. cppsubrules = [ ("hier::Box<NDIM>",       "Box"),
#    ("hier::IntVector<NDIM>", "IntVect"),
#    ("hier::Index<NDIM>",     "IntVect") ]

pythonsubrules = []

fortransubrules = []

F90subrules = []

F90_90subrules = []

MPPLsubrules = []

try:
  execfile("localrules.py")
except:
  try:
    execfile("../localrules.py")
  except:
    print("No file or problem with 'localrules.py'; proceeding with default rules")

def fnconvert(name,suffixout):
    # Converts a file name "name" by substituting suffixout for the
    # existing suffix, or if there is no suffix, appending suffixout
    suffixin = name.split(".")[-1]
    if (suffixin != name):
        nameout = name.rstrip(suffixin)+suffixout
        # this coding makes sure only last occurence of suffix is repalced
    else:
        nameout = name+"."+suffixout
    return nameout

from stat import ST_MTIME
from numpy import greater
def newer(file0,file1):
    # returns 1 if file0 is newer than file1, 0 otherwise
    # if file0 does not exist, raise a standard python error
    # if file1 does not exist, return 2
    time0 = os.stat(file0)[ST_MTIME]
    try:
        time1 = os.stat(file1)[ST_MTIME]
    except:
        return 2
    return greater(time0,time1)

class processdirs:
    # Process a list of directories to do conversions with processall
    # The default list of directories is ".", the current directory
    # Establishes a set of subdirectories with processed files
    # in a root directory whose default is called "converted" and
    # is parallel to ".".
    def __init__(self,indirs = ["."],outroot = "../converted",clean=None):
        # create the output root directory if it needs to be
        curpath = os.getenv("PWD")
        try:
            os.mkdir(outroot)
            print "Creating directory "+outroot
        except:
            try:
                os.chdir(outroot)
                print "Output directory already exists; proceeding"
                os.chdir(curpath)
            except:
                raise "Can't create output directory"
        for direc in indirs:
            print "Entering directory "+direc
            if direc == ".":
                curdir = curpath.split("/")[-1]
                outdir = outroot+"/"+curdir
            else:
                outdir = outroot + "/"+direc
            processall(direc,outdir,clean)

class processall:
    # Process all files in a directory with rules according to file type
    # as designated by suffix.
    # List of classes of distict file types
    classlist = ["Py","Cpp","Cpp_h","Cpp_hh","Fortran","F90","F90_90","MPPL"]
    def __init__(self,indir=".",outdir="./converted",clean=None):
        for entry in self.classlist:
            a=eval(entry+"("+`indir`+","+`outdir`+",clean="+`clean`+")")
            # Only process file types with non-empty substition rules
            if (a.subrules != []):
                print "processing for file type ",entry
                a.process()

class generic:
    """
generic class for replacing strings in a series of files.
Methods:
   process(): processes the list of files in the directory indir with
     suffix suffixin, and for each one creates a new file in directory
     outdir with the same root name appended by
     the string suffixout (which by default is suffixin).
   processfile(filename): processes a specific file
Notes:
   indir is string with path.  By default it is ".", the current directory
   outdir is by default "./converted".   If the output directory does not
      exist it will be created.  Can be overwritten after instance created
    """

    suffixin = ""
    suffixout = ""
    subrules = globalsubrules
    def __init__(self,indir=".",outdir="./converted",clean=None):
        self.indir = indir
        self.outdir = outdir
        self.clean=clean
        self.doclean = 0   # default is no removal of duplicate files
        if (self.suffixout == ""):
            self.suffixout = self.suffixin
        # create the output directory if it hasn't been created
        try:
            os.mkdir(self.outdir)
        except:
            pass
    def process(self):
        # get the list of files to process
        self.filelist = filesublist(self.indir,self.suffixin)
        print "processing directory",self.indir,"to directory",self.outdir
        # exclude convert.py, localrules.py and filelists.py from filelist
        if (self.suffixin == "py"):
            try:
                self.filelist.remove("convert.py")
            except:
                pass
            try:
                self.filelist.remove("localrules.py")
            except:
                pass
            try:
                self.filelist.remove("filelists.py")
            except:
                pass
        # Set a flag to remove unchanged files if input and output directories are different
        if self.clean and (self.outdir != self.indir):
            # if the input and output directories are different, eliminate files
            # from the output directory that duplicate those in the input directory.
            self.doclean = 1
        else:
            self.doclean = 0
        # Now process files remaining
        for file in self.filelist:
            self.processfile(file)
    def processfile(self,filename):
        # Convert an individual file
        self.infile=self.indir+"/"+filename
        if (self.suffixout == self.suffixin):
            outfilename=filename
        else:
            outfilename=fnconvert(filename,self.suffixout)
        self.outfile = self.outdir+"/"+outfilename
        if (newer(self.infile,self.outfile)):
            # Only process infile if it is newer than outfile
            print "converting file "+self.infile+" to "+self.outfile
            f=open(self.infile,"r")
            g=open(self.outfile,"w")
            lines=f.readlines()
            iline = 0   # just a diagnostic, counts lines
            for line in lines:
              iline = iline + 1
              for rule in self.subrules:
                if (line != None):
                  if type(rule)==types.TupleType or type(rule)==types.ListType:
                    line = string.replace(line,rule[0],rule[1])
                  else:
                    # if it's not a tuple or list, assume it is a method
                    # that executes something
                    # more complicated than a simple replace
                    line = rule(line)
              if (line != None):
                  if (type(line) == types.StringType):
                      g.write(line)
                  if (type(line) == types.ListType):
                      g.writelines(line)
                    
            g.close()
            f.close()
            if self.doclean == 1:
              if filecmp.cmp(self.outfile,self.infile) and self.outfile != self.infile:
                print self.outfile + " is the same as " + self.infile + ", removing"
                os.remove(self.outfile)

class Py(generic):
    suffixin="py"
    subrules = globalsubrules + pythonsubrules

class Cpp(generic):
    suffixin = "C"
    subrules = globalsubrules + cppsubrules

class Cpp_h(Cpp):
    suffixin = "h"

class Cpp_hh(Cpp):
    suffixin = "hh"

class Fortran(generic):
    suffixin = "f"
    subrules = globalsubrules + fortransubrules

class F90(generic):
    suffixin = "F"
    subrules = globalsubrules + F90subrules

class F90_90(F90):
    suffixin = "F90"

class MPPL(generic):
    suffixin = "m"
    subrules = globalsubrules + MPPLsubrules

