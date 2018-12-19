## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am,v 7.0 2018/02/28 18:32:52 meyer8 Exp $
##
## UEDGE project
##
## ############################################################

BUILDDIR	= ../build
FORTHON_INCS	= @FORTHON_INCS@
SO	= so

SOURCES = \
	flxcomp.m \
	flxdriv.m \
	flxread.m \
	flxwrit.m \
	dvshdf5.m

DISTSOURCES = \
	flx.v

EXTRA_DIST = \
	$(SOURCES) \
	$(DISTSOURCES)

INCLUDES = \
	-I$(srcdir) \
	-I$(top_builddir)/bbb \
	-I$(top_builddir)/api \
	-I$(top_builddir)/aph \
	-I$(top_builddir)/com \
	-I$(top_builddir) \
	$(FCIOWRAPPERS_INC) \
        $(HDF5_MODINC)

.m.F:
	$(top_srcdir)/pyscripts/convert1.py -i $(srcdir) -o \
	  $$PWD `basename $<`

all-local: dist-local $(BUILDDIR)/flxpy.$(SO)

../com/com.scalars:
	cd ../com && make com.scalars

$(BUILDDIR)/flxpy.$(SO): \
  	flxcomp.F flxdriv.F flxread.F flxwrit.F dvshdf5.F \
  	flx.v ../com/com.v
	echo Makefile CC = $(CC); \
	source $(top_srcdir)/debugenv.sh; \
	$(FORTHON_CMD) $(INCLUDES) --builddir $(BUILDDIR) \
	  -d com -f flxcomp.F flx flxdriv.F flxread.F flxwrit.F \
	  dvshdf5.F
	if test -f flxpy.so -a $(BUILDDIR) != .; then \
	  mv flxpy.so $(BUILDDIR) ; \
        fi

clean-local:
	rm -f *.F forthonf2c.h *.mod
	for i in $(DISTSOURCES); do \
	  if test -L $$i; then \
	    rm -f $$i; \
	  fi \
	done
	-cd $(BUILDDIR) && rm -f flxcomp.o flxdriv.o flxread.o  flxwrit.o dvshdf5.o \
	  aeqdskv.mod dimflx.mod flxin.mod flx_p.F90 flx_p.o flxpymodule.c \
	  flxpymodule.o flxpy.so flx.scalars inpf0.mod inpf.mod Makefile.flx \
	  neqdskv.mod polflx.mod workdn.mod

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
