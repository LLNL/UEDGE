## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am,v 7.0 2018/02/28 18:59:58 meyer8 Exp $
##
## UEDGE project
##
## ############################################################

BUILDDIR	= ../build
FORTHON_INCS	= @FORTHON_INCS@
SO	= so

SOURCES = \
	wdf.m

DISTSOURCES = \
	wdf.v

EXTRA_DIST = \
	$(SOURCES) \
	$(DISTSOURCES)

INCLUDES = \
	-I$(srcdir) \
	-I$(top_builddir)/svr \
	-I$(top_builddir)/grd \
	-I$(top_builddir)/flx \
	-I$(top_builddir)/bbb \
	-I$(top_builddir)/api \
	-I$(top_builddir)/aph \
	-I$(top_builddir)/com \
	-I$(top_srcdir)

.m.F:
	$(top_srcdir)/pyscripts/convert1.py -i $(srcdir) -o \
	  $$PWD `basename $<`

all-local: dist-local $(BUILDDIR)/wdfpy.$(SO)


$(BUILDDIR)/wdfpy.$(SO): wdf.F wdf.v
	echo Makefile CC = $(CC); \
	source $(top_srcdir)/debugenv.sh; \
	$(FORTHON_CMD) $(INCLUDES) --builddir $(BUILDDIR) wdf
	if test -f wdfpy.$(SO) -a $(BUILDDIR) != .; then \
	  mv wdfpy.$(SO) $(BUILDDIR) ; \
	fi

clean-local:
	rm -f $(BUILDDIR)/wdfpy.$(SO) *.F forthonf2c.h *.mod
	for i in $(DISTSOURCES); do \
	  if test -L $$i; then \
	    rm -f $$i; \
	  fi \
	done
	-cd $(BUILDDIR) && rm -f wdf.o auxw.mod degas1.mod degas2.mod \
	  dimwdf.mod eqdsk.mod linkgrd.mod Makefile.wdf nodew.mod \
	  options.mod wdf_p.F90 wdf_p.o wdfpymodule.c wdfpymodule.o \
	  wdf.scalars

uninstall-local:
	rm -f $(prefi$(BUILDDIR)/wdfpy.$(SO)

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
