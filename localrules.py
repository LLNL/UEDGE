#rules for converting mppl to f90
# This is generic down to the UEDGE section
import copy
import os
def Use2use(s):
    #Return s if this is a comment
    if (not s[0].isspace()) and s[0] != "U":
        return s
    #Do a substitution if line contains "Use" and it is not
    #part of a comment
    if (s.find("Use")+1):
        if (s.find("!")==-1) or (s.find("Use")<s.find("!")):
            sout=s.replace("("," ")
            sout=sout.replace(")"," ")
            sout=sout.replace("Use","      use")
            return sout
    return s
def Allot(s):
    # converts MPPL allot to F90 allocate
    if (s.find("allot")+1):
        if (s.find("call allot")+1):
            s=s.replace("call allot","allocate")
            s=s.replace('"','')
            s=s.replace("'","")
            s=s.replace(",","(")
            s=s.replace(")","))")
            return s
        s=s.replace(", allot","")
        s=s.replace(",allot","")
    return s
def Nopdb(s):
    if os.getenv('PACT_DIR') == None:
        if (s.find("c!nopdb")+1):
            s=s.replace("c!nopdb","       ")
    return s
def Petsc(s):
    if os.getenv('PETSC_DIR') != None:
        if os.getenv('PARALLEL') == None:
           if (s.find("cunipetsc")+1):
               s=s.replace("cunipetsc","")
        if (s.find("cpetsc")+1):
            s=s.replace("cpetsc","")
    return s

saved_dec=0
in_uses=0
savedlines=[]
def MoveDecs(s):
    global saved_dec,in_uses
    # Return if this is a comment
    if (not s[0].isspace()) and (not saved_dec) and (not in_uses):
        return s
    # collect lines for declarations
    # if we find an "implicit none" statement, store it and remove the line
    sls=s.lstrip().lower()
    indfunc=sls.find("function")
    indcom=sls.find("!")
    functest= (indfunc == -1 or -1<indcom<indfunc)
    # tests to exclude "real function" but allow "real ! function as
    # part of declaration block
    if (sls[0:8]=="implicit") or (sls[0:4]=="real" and functest) \
       or (sls[0:7]=="integer" and functest) \
       or (sls[0:9]=="character") or (sls[0:9]=="parameter") or \
       (sls[0:8]=="external") or (sls[0:9] == "intrinsic") or \
       (sls[0:7]=="logical" and functest) or (sls[0:9]=="dimension") or \
       (sls[0:4] == "data"):
        savedlines.append(s)
        saved_dec=1
        in_uses=0
        return None
    # if we are in the midst of declarations, save also comments (except for
    # "Common block") and continuations and blank lines as part of
    #  what is moved.)
    if (saved_dec==1) and (sls == "" or s[0].lower() == "c" or sls[0]=="!" \
      or s[0]=="*") and (in_uses == 0) \
      and (s.find("Common block")==-1):
        savedlines.append(s)
        return None
    # Check for continuation line in midst of declarations
    if (saved_dec==1) and (len(s)>6):
        if (s[5].isspace() == 0):
            savedlines.append(s)
            return None

    if (sls[0:3] == "use"):
        in_uses=1

    if (saved_dec==1) and (sls != "") and s[0] != "c" and sls[0] != "!" and \
       (sls[0:3] != "use"):
        #This is our first executable statement.  Add it to our saved
        # declarations lines and return them now
        templines = copy.copy(savedlines)
        templines.append(s)
        #empty out savedlines
        del savedlines[0:]
        saved_dec = 0
        in_uses=0
        return templines
    return s
        
inelseif = 0
savedelselines=""
def Elseifthen(s):
    # put a "then" at the end of an elseif if it isn't already there
    # need to check to see if next line is a continue
    global inelseif,savedelselines
    # return s if this is a comment
    if (not s[0].isspace()):
        return s
    if s.find("elseif")+1:
        if s.find("then")+1:
            return s
        # set a flag that we are in an "in-else-if" block that needs work
        inelseif = 1
        # If there is no "then" we need to save it to test if the next
        # line is a continuation
        savedelselines = s
        return(None)
    if (inelseif and len(s)>6 and not s[5].isspace()):
        # This is a continue line, add it to savedelselines
        savedelselines += s 
        return(None)
    if (inelseif and (len(s)<6 or s[5].isspace())):
        # No longer in a continue, so process lines
        if savedelselines.find("then")+1:
            savedelselines +=  s
            inelseif = 0
            return savedelselines
        if savedelselines.split("\n")[-2].find("!")+1:
            # if last line in saved lines has a comment,
            #  find index of last comment sign
            last = savedelselines.rfind("!")
            savedelselines=savedelselines[0:last]+ \
                " then "+savedelselines[last:] + s
            inelseif=0
            return savedelselines
        #Otherwise the last line has no comment so insert "then" at end
        savedelselines = savedelselines[0:-1]+" then\n" + s
        inelseif=0
        return savedelselines
    return s

M2Fsubrules = [("#","!"),Use2use,
               ("c!ifdef","#ifdef"),
               ("c!else","#else"),
               ("c!endif","#endif"),
               ("(Size4)","(kind=4)::"),
               (":: function"," function"),
               (" break "," exit "),
               (" break\n"," exit\n"),
               ("while (","do while ("),
               ("endwhile","end do"),
               ("      call ruthere","c     call ruthere"),
               ("c!include ","#include "),
               Nopdb,
               Petsc,
               Allot,
               Elseifthen,
               MoveDecs
               ]

#-------------------------------------
# Special for UEDGE
wordsizectr=0
def grdproc(s):
    # process to eliminate ifelse write construction
    global wordsizectr
    if (s.find("ifelse([WORDSIZE]")+1):
        s="#if WORDSIZE == 64\n 2001 format(1p3e23.15)\n#else\n 2001 format(1p3d23.15)\n#endif\n"
#        s="#ifndef WORDSIZE\n 2001 format(1p3e23.15)\n#else\n 2001 format(1p3d23.15)\n#endif\n"
        wordsizectr=4
        return s
    elif wordsizectr > 0:
        wordsizectr -= 1
        return None
    else:
        wordsizectr = 0
        return s

M2Fsubrules.insert(6,grdproc)
M2Fsubrules.insert(3,("do i1=","do_i1: do i1="))
M2Fsubrules.insert(4,("break (2) ! exit do_i1","exit do_i1"))
M2Fsubrules.insert(5,("enddo ! do_i1","enddo do_i1"))
M2Fsubrules.insert(6,("float","real"))
M2Fsubrules.insert(6,("dfloat","real"))


