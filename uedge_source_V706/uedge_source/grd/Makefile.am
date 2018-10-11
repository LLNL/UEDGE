## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am,v 7.0 2018/02/28 18:32:53 meyer8 Exp $
##
## UEDGE project
##
## ############################################################

BUILDDIR	= ../build
FORTHON_INCS	= @FORTHON_INCS@
SO	= so

SOURCES = \
	grdcomp.m \
	grddriv.m \
	grdinit.m \
	grdread.m \
	grdwrit.m

DISTSOURCES = \
	grd.v

EXTRA_DIST = \
	$(SOURCES) \
	$(DISTSOURCES)


INCLUDES = \
	-I$(srcdir) \
	-I$(top_builddir)/flx \
	-I$(top_builddir)/bbb \
	-I$(top_builddir)/api \
	-I$(top_builddir)/aph \
	-I$(top_builddir)/com \
	-I$(top_srcdir)

.m.F:
	$(top_srcdir)/pyscripts/convert1.py -i $(srcdir) -o \
	  $$PWD `basename $<`

all-local: dist-local $(BUILDDIR)/grdpy.$(SO)

../com/com.scalars :
	cd ../com ; make com.scalars


$(BUILDDIR)/grdpy.$(SO): grdcomp.F grddriv.F grdinit.F grdread.F grdwrit.F \
  	grd.v ../com/com.v
	echo Makefile CC = $(CC); \
	source $(top_srcdir)/debugenv.sh; \
	$(FORTHON_CMD) $(INCLUDES) --builddir $(BUILDDIR) \
	  -d com -f grdcomp.F grd grddriv.F grdinit.F \
	  grdread.F grdwrit.F
	if test -f grdpy.so -a $(BUILDDIR) != .; then \
          mv grdpy.so $(BUILDDIR); \
        fi

clean-local:
	rm -f *.F forthonf2c.h *.mod
	for i in $(DISTSOURCES); do \
	  if test -L $$i; then \
	    rm -f $$i; \
	  fi \
	done
	-cd $(BUILDDIR) && rm -f grdcomp.o grddriv.o grdinit.o grdread.o \
	  grdwrit.o grdpy.$(SO) aeqdskg.mod analgrd.mod argfc.mod \
	  curves.mod dimensions.mod dnull_temp.mod efitg.mod grd_p.F90 \
	  grd_p.o grdpymodule.c grdpymodule.o grd.scalars inmesh.mod \
	  linkco.mod Makefile.grd mmod.mod refinex.mod spline.mod \
	  system.mod torannulus.mod transfm.mod transit.mod uegrid.mod \
	  xmesh.mod

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
