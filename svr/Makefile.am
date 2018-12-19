## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am,v 7.0 2018/02/28 18:59:57 meyer8 Exp $
##
## UEDGE project
##
## ############################################################

BUILDDIR	= ../build
FORTHON_FLAGS	=
FORTHON_INCS	= @FORTHON_INCS@
## PETSC_INC	= @PETSC_INC@
SO	= so

if HAVE_PETSC
PETSC_F90	= petsc-uedge.F90
endif

SOURCES = \
	daspk.m \
	nksol.m \
	svrut1.m \
	svrut2.m \
	svrut3.m \
	svrut4.m \
	uoa.m \
	vodpk.m

DISTSOURCES = \
	svr.v

DISTSOURCES += petsc-uedge.F90

EXTRA_DIST = \
	$(SOURCES) \
	$(DISTSOURCES)

INCLUDES = \
	-I$(srcdir) \
	-I$(top_srcdir)/grd \
	-I$(top_srcdir)/flx \
	-I$(top_srcdir)/bbb \
	-I$(top_srcdir)/api \
	-I$(top_srcdir)/aph \
	-I$(top_srcdir)/com \
	-I$(top_srcdir)

.m.F:
	$(top_srcdir)/pyscripts/convert1.py -i $(srcdir) -o \
	  $$PWD `basename $<`

all-local: dist-local $(BUILDDIR)/svrpy.$(SO)

$(BUILDDIR)/svrpy.$(SO): \
  	daspk.F nksol.F svrut1.F svrut2.F svrut3.F svrut4.F \
  	vodpk.F uoa.F \
  	$(PETSC_F90) \
  	svr.v
	echo Makefile CC = $(CC); \
	source $(top_srcdir)/debugenv.sh; \
	$(FORTHON_CMD) $(INCLUDES) --builddir $(BUILDDIR) \
	  --noimplicitnone \
	  $(FORTHON_FLAGS) -f daspk.F svr nksol.F svrut1.F \
	  svrut2.F svrut3.F svrut4.F vodpk.F uoa.F \
	  $(PETSC_F90)
	if test -f svrpy.$(SO) -a $(BUILDDIR) != .; then \
          mv svrpy.$(SO) $(BUILDDIR) ; \
        fi

clean-local:
	rm -f *.F forthonf2c.h *.mod
	for i in $(DISTSOURCES); do \
	  if test -L $$i; then \
	    rm -f $$i; \
	  fi \
	done
	-cd $(BUILDDIR) && rm -f daspk.o nksol.o svrut1.o svrut2.o \
	  svrut3.o svrut4.o uoa.o vodpk.o petsc-uedge.o svrpy.$(SO) \
	  Makefile.svr petscmodinterface2.mod petscmodinterfaces.mod \
	  petscmod.mod snesmonitormod.mod svr_p.F90 svr_p.o svrpymodule.c \
	  svrpymodule.o svr.scalars uoa.mod

# bootstrap to make Forthon play nicely
dist-local:
	for i in $(DISTSOURCES); do \
	  if test ! -f $$i; then \
	    ln -s $(srcdir)/$$i .; \
	  fi \
	done
	if test ! -e ../localrules.py; then \
	  ln -s `(cd $(srcdir)/..; pwd -P)`/localrules.py ../; \
	fi
