## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am,v 7.0 2018/02/28 18:32:43 meyer8 Exp $
##
## UEDGE project
##
## ############################################################

BUILDDIR	= ../build
FORTHON_INCS	= @FORTHON_INCS@
SO	= so

SOURCES = \
	aphrates.m \
	aphread.m

DISTSOURCES = \
	aph.v

EXTRA_DIST = \
	$(SOURCES) \
	$(DISTSOURCES)

INCLUDES = \
	-I$(srcdir) \
	-I$(top_builddir)/com \
	-I$(top_srcdir)

.m.F:
	$(top_srcdir)/pyscripts/convert1.py -i $(srcdir) -o \
	  $$PWD `basename $<`

all-local: dist-local $(BUILDDIR)/aphpy.$(SO)


$(BUILDDIR)/aphpy.$(SO): aphrates.F aphread.F aph.v \
	$(BUILDDIR)/compymodule.o
	echo Makefile CC = $(CC); \
	source $(top_srcdir)/debugenv.sh; \
	$(FORTHON_CMD) $(INCLUDES) --builddir $(BUILDDIR) \
	  -f aphrates.F aph aphread.F
	if test -f aphpy.so -a $(BUILDDIR) != .; then \
	  mv aphpy.so $(BUILDDIR) ; \
	fi

clean-local:
	rm -f $(BUILDDIR)/aphpy.$(SO) aphrates.F aphread.F forthonf2c.h *.mod
	for i in $(DISTSOURCES); do \
	  if test -L $$i; then \
	    rm -f $$i; \
	  fi \
	done
	-cd $(BUILDDIR) && rm -f aphrates.F aphread.F aph_p.F90 aph_p.o \
	  aphpymodule.c aphpymodule.o aphrates.o aphread.o aph.scalars \
	  aphwrk.mod data_input.mod ionization_energy.mod Makefile.aph \
	  physical_constants.mod rtdata.mod rtdegas.mod

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
