#
#
# idf_tst.make is a makefile to create test programs for IDF package
#
# usage:       make -f idf_test.make testXXX
# where testXXX stands for a name of program
#
# libidf.a library is in ../lib 
# when created by Makefile located in the home IDF directory
#
# 
#=======================================================
#       Import of Data from File (IDF) package
#     for time-saving programming of data input
#              to scientific codes.
#
#                 IDF Version 1.1   
#              Date: February 1, 2000
#
#=======================================================
#
#

# specify suffixes
A         = .a   # for object module library
O         = .o   # for object file
EXE       = .exe # for executable module

# specify C compliler
CC        = gcc
COPT      = -O2   # optimization option
CFLAGS    = -o -ansi

# specify C++ compliler
CCCC      = g++
CCOPT     = -O3   # optimization option
CCFLAGS   = -o

# specify FORTRAN compliler
FC        = g77
FOPT      = -O    # optimization option
FFLAGS    = -o

# specify FORTRAN compiler option for aggregate target alignment
FALIGN    =  
#FALIGN    =  -align dcommon  # DIGITAL alpha
#FALIGN    =  -dalign         # sun 


# specify the name (and directory) for IDF library.
IDFLIB  = ../lib/libidf$A


# specify pathway to directory in which the IDF header file
# "idfuser.h" was put.
IDFINC  = ../include/


LIBS = ${IDFLIB} -lm

TC1	= testc1${EXE}
TC2	= testc2${EXE}
TCC	= testcc${EXE}
TF1	= testf1${EXE}
TF2	= testf2${EXE}
TC      =   test${EXE}
	
all	 : testc testcpp testfor
	
# C programs for IDF test
testc	 : testc1 testc2 test
	
# C++ programs for IDF test
testcpp	 : testcc
	
# FORTRAN-77 programs for IDF test
testfor  : testf1 testf2
	
testc1	 :
	 ${CC} ${COPT} ${TARGET_ARCH} -I${IDFINC} idfsample1.c \
         ${LIBS} ${CFLAGS} ${TC1}
	
testc2	 :
	 ${CC} ${COPT} ${TARGET_ARCH} -I${IDFINC} idfsample2.c \
	 ${LIBS} ${CFLAGS} ${TC2}
	
testcc	 :
	 ${CCCC} ${CCOPT} ${TARGET_ARCH} -I${IDFINC} idfsample.cc \
         ${LIBS} ${CCFLAGS} ${TCC}
	
testf1	 :
	 ${FC} ${FOPT} ${TARGET_ARCH} idfsample1.f \
         ${LIBS} ${FFLAGS} ${TF1}
	
testf2	 :
	 ${FC} ${FOPT} ${TARGET_ARCH} ${FALIGN} idfsample2.f \
         ${LIBS} ${FFLAGS} ${TF2}
	
test	 :
	 ${CC} ${COPT} ${TARGET_ARCH} -I${IDFINC} idftst.c \
         ${LIBS} ${CFLAGS} ${TC}
	
delete	 :
	 rm ${TC} ${TC1} ${TC2} ${TCC} ${TF1} ${TF2}
	
