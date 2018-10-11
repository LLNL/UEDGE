##########################################################################
#
# Pyuedge/Makefile -- top-level makefile for building uedge as a shared-
#   object python extension. Directory structure uedge/Pyuedge is assumed.
# User must set TOPALL below to parent directory of /uedge & /pyUtils.
#
#   $Id: Makefile,v 7.0 2018/02/28 18:32:42 meyer8 Exp $
#   $Id: Makefile,v 7.0 2018/02/28 18:32:42 meyer8 Exp $
#
##########################################################################


# comment out the following line to make uedge.so wo python's distutils
USEdistutils = $(DEFINED)

# here must set path to directory containing /uedge and /pyUtils 
TOPALL = /usr/local/Uedge/Ver_4.41_mds

# TOP - the uedge directory.
TOP = $(TOPALL)/uedge
##TOP = ..
 
CONFIGP  = '$$(TOP)'/Pyuedge/config
CONFIG   = $(subst '$$(TOP)',$(TOP),$(CONFIGP))
#CONFIGP  = '$$(TOP)'/../config
##CONFIG   = $(CONFIGP)

PYUTILSP = $(TOPALL)/pyUtils
PYUTILS  = $(PYUTILSP)
#PYUTILSP = '$$(TOP)'/../pyUtils
#PYUTILS  = $(subst '$$(TOP)',$(TOP),$(PYUTILSP))

include $(CONFIG)/Makefile.tmpl

# make-variables passed to make's below are inherited by *all* subsequent make's

PASS_ALL     = ARCH_INC="$(ARCH_INC)" CPUNAME="$(CPUNAME)" DEBUG="$(DEBUG)" \
               CONFIG=$(CONFIGP) PYUTILS=$(PYUTILSP)
PASSFLAGS    = TOP=$(TOP) $(PASS_ALL)
##PASSFLAGS_X1 = TOP=../$(TOP) $(PASS_ALL)
PASSFLAGS_X1 = TOP=$(TOP) $(PASS_ALL)

Packages = com bbb aph api flx grd
Dirs     = Csor Fsor svr

FFLAGS = $(FFFLAGS)
CFLAGS = $(CCFLAGS) $(PYTHON_INCLUDES)

ifndef USEdistutils
  PYTHONSOS = $(PYTHON_ROOT)/lib/python$(PYVERS)/site-packages/Numeric/*.so \
              $(PYTHON_ROOT)/lib/python$(PYVERS)/site-packages/Numeric/*/*.so
  MKMOD.O =  MKMOD.O=DEFINED
endif

PYTHON_LINCLUDES = -I$(PYWRAPPER)

all:: dirtst $(PYTHON_INSTALL)/uedge.so

ifdef USEdistutils
# can't control python setup.py dependencies.
# rm below is too much if there is more than one platform under build
# and too little if  ...pymodule.c files are new (then rm build/*/*.o).
# Extension(<name>, makes <name>.so
$(PYTHON_INSTALL)/uedge.so: Packages Csor/$(ARCH)/handlers.o \
                            Fsor/$(ARCH)/aux_tot.o $(ARCH)/uedgepymod.o
	rm -f $(PYTHON_INSTALL)/uedge*.so build/*/uedge*.so
	env PYUTILS='$(PYUTILS)' FORTLIB='$(FORTLIB)' ARCH='$(ARCH)' \
	  Packages='$(Packages)' \
	  $(PYTHON) setup.py install --install-lib=$(PYTHON_INSTALL)

else
$(PYTHON_INSTALL)/uedge.so: Packages Csor/$(ARCH)/handlers.o \
                            Fsor/$(ARCH)/aux_tot.o $(ARCH)/uedgepymod.o
	$(LD) $(LDFLAGS) -o $(PYTHON_INSTALL)/uedge.so \
	 $(foreach pp,$(Packages),$(pp)/$(ARCH)/$(pp).sor.o) \
	 $(foreach pp,$(Packages),$(pp)/$(ARCH)/$(pp).p.o) \
	 $(foreach pp,$(Packages),$(pp)/$(ARCH)/$(pp)pymodule.o) \
	 Csor/$(ARCH)/handlers.o Fsor/$(ARCH)/aux_tot.o $(ARCH)/uedgepymod.o \
	 $(PYWRAPPER)/../$(ARCH)/pybasisC.so $(PYTHONSOS) $(FORTLIB)

endif

Packages:
	@for pp in $(Packages); \
	do \
	  cd $$pp; $(MAKE) $(PASSFLAGS_X1) PKG=$$pp $(MKMOD.O); cd ..; \
	done

Fsor/$(ARCH)/aux_tot.o: Fsor/$(ARCH)/aux_tot.f 
	cd Fsor/$(ARCH); $(FORT) $(FFLAGS) -c -w aux_tot.f 

Csor/$(ARCH)/handlers.o: Csor/$(ARCH)/handlers.c
	cd Csor/$(ARCH); $(CC) $(CFLAGS) -c handlers.c; cd ..

Fsor/$(ARCH)/aux_tot.f: svr/$(ARCH)/svr.sor.f Fsor/ex_das_isa.f $(TOP)/svr/daux1.f \
	  $(TOP)/com/dsum.f $(TOP)/com/error.f $(TOP)/com/getmsg.f $(TOP)/com/dummy_py.f
	cat $(TOP)/com/dsum.f $(TOP)/com/error.f $(TOP)/com/getmsg.f \
            $(TOP)/com/dummy_py.f  Fsor/ex_das_isa.f > Fsor/$(ARCH)/aux_tot.f
	cat $(TOP)/svr/daux1.f >> Fsor/$(ARCH)/aux_tot.f
	cat svr/$(ARCH)/svr.sor.f >> Fsor/$(ARCH)/aux_tot.f

Csor/$(ARCH)/handlers.c: $(TOP)/com/handlers.c
	cp $(TOP)/com/handlers.c Csor/$(ARCH)/handlers.c

svr/$(ARCH)/svr.sor.f: svr/$(ARCH)/svr.sor.m 
	$(MPPL_ROOT)/bin/mppl -I$(MPPL_ROOT)/include -DNO_INCLUDES=1 -t$(ARCH_PB) \
	    -r$(RealSize) -DCOMPILER=F90 $(TOP)/$(Macros) $(OTHER_DFILES) \
	    svr/$(ARCH)/svr.sor.m > svr/$(ARCH)/svr.sor.f.tmp 2>svr/$(ARCH)/svr.sor.f.err
	$(MV) svr/$(ARCH)/svr.sor.f.tmp svr/$(ARCH)/svr.sor.f

svr/$(ARCH)/svr.sor.m: $(TOP)/svr/*.m 
	cat $(TOP)/svr/*.m > $@

$(ARCH)/%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:: clean.distutils clean.pyUtils.rt
	@for pp in $(Packages); \
	do \
	  cd $$pp; $(MAKE) $(PASSFLAGS_X1) PKG=$$pp clean; cd ..; \
	done
	rm -f $(ARCH)/*.o $(ARCH)/*.so \
	  svr/$(ARCH)/svr.sor.f svr/$(ARCH)/svr.sor.m \
	  Csor/$(ARCH)/*.o Fsor/$(ARCH)/aux_tot.f Fsor/$(ARCH)/*.o

dirtst:
	@for pp in $(Dirs); \
	do \
	  cd $$pp; if [ ! -d $(ARCH) ]; then mkdir $(ARCH); fi; cd ..; \
	done

make.all: pyUtils all

pyUtils:
#	cd $(PYUTILS); $(MAKE) TOP=$(TOP)/uedge $(PASS_ALL)
	cd $(PYUTILS); $(MAKE) TOP=$(TOP) $(PASS_ALL)

clean.all:: clean
#	cd $(PYUTILS); $(MAKE) TOP=$(TOP)/uedge $(PASS_ALL) clean
	cd $(PYUTILS); $(MAKE) TOP=$(TOP) $(PASS_ALL) clean

.PHONY:: dirtst Packages pyUtils
