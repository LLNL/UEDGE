## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am,v 7.0 2018/02/28 18:32:48 meyer8 Exp $
##
## UEDGE project
##
## ############################################################

BUILDDIR	= ../build
FORTHON_INCS	= @FORTHON_INCS@
FC_PIC_FLAG	= @FC_PIC_FLAG@
SO	= so

if HAVE_PETSC
PETSC_F90	= petscMod.F90
endif


SOURCES = \
	blasext.m \
	brent.m \
	comutil.m \
	misc.m \
	mnbrak.m \
	mpivarsmod.F \
	$(PETSC_F90)

DISTSOURCES = \
	handlers.c \
	handlers_CRAY.c \
	handlers_SP.c \
	vector.c \
	com.h \
	llnltyps.h \
	dsum.f \
	dummy_py.f \
	error.f \
	getmsg.f \
	ssum.f \
	com.v

DISTSOURCES += petscMod.F90

EXTRA_DIST = \
	$(SOURCES) \
	$(DISTSOURCES)

INCLUDES = \
	-I$(srcdir) \
	-I$(top_srcdir)

.m.F:
	$(top_srcdir)/pyscripts/convert1.py -i $(srcdir) -o $$PWD `basename $<`

all-local: dist-local $(BUILDDIR)/compy.$(SO)

$(BUILDDIR)/compy.$(SO): \
	blasext.F brent.F comutil.F misc.F mnbrak.F \
	dsum.f dummy_py.f error.f getmsg.f ssum.f \
	com.v handlers.c mpivarsmod.F $(PETSC_F90)
	echo Makefile CC = $(CC); \
	source $(top_srcdir)/debugenv.sh; \
	$(FORTHON_CMD) $(INCLUDES) --builddir $(BUILDDIR) -I.. \
	  -f blasext.F com brent.F comutil.F misc.F mnbrak.F \
	  dsum.f dummy_py.f error.f getmsg.f ssum.f handlers.c \
	  mpivarsmod.F $(PETSC_F90)
	if test -f compy.so -a $(BUILDDIR) != .; then \
	  mv compy.so $(BUILDDIR) ; \
	fi

clean-local:
	rm -f *.F forthonf2c.h mpi_nim.mod handlers.o *.mod
	for i in $(DISTSOURCES); do \
	  if test -L $$i; then \
	    rm -f $$i; \
	  fi \
	done
	- cd $(BUILDDIR) && rm -f blasext.o brent.o compymodule.o \
	  comutil.o misc.o mnbrak.o dsum.o dummy_py.o error.o getmsg.o \
	  ssum.o compymodule.c com_p.o com_p.F90 com.scalars petscMod*.o \
	  comgeo.mod compy.so dim.mod fitdata.mod Forthon.c forthonf2c.h \
	  Forthon.h Forthon.o limiter.mod linkbbb.mod Makefile.com \
	  multicharge.mod noggeo.mod petsc_snes_param.mod \
	  rz_grid_global.mod rz_grid_info.mod share.mod timespl.mod \
	  timing.mod xpoint_indices.mod

uninstall-local:
	rm -f $(prefix)/lib/compy.$(SO)

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

