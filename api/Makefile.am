## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am,v 7.0 2018/02/28 18:32:44 meyer8 Exp $
##
## UEDGE project
##
## ############################################################

BUILDDIR	= ../build
FORTHON_INCS	= @FORTHON_INCS@
SO	= so

SOURCES = \
	apifcn.m \
	apip93.m \
	apisorc.m \
	fimp.m \
	fmombal.m \
	inelrates.m \
	sputt.m

DISTSOURCES = \
	api.h \
	api.v

EXTRA_DIST = \
	$(SOURCES) \
	$(DISTSOURCES) babel.make


INCLUDES = \
	-I$(srcdir) \
	-I$(top_builddir)/aph \
	-I$(top_builddir)/com \
	-I$(top_srcdir)

.m.F:
	$(top_srcdir)/pyscripts/convert1.py -i $(srcdir) -o \
	  $$PWD `basename $<`

all-local: dist-local $(BUILDDIR)/apipy.$(SO)


../com/com.scalars:
	(cd ../com ; make com.scalars )

$(BUILDDIR)/apipy.$(SO): \
  	apifcn.F apip93.F apisorc.F fimp.F fmombal.F \
  	inelrates.F sputt.F api.v ../com/com.v $(BUILDDIR)/compymodule.o
	echo Makefile CC = $(CC); \
	source $(top_srcdir)/debugenv.sh; \
	$(FORTHON_CMD) $(INCLUDES) --builddir $(BUILDDIR) \
	  --noimplicitnone \
	  -d com -f apifcn.F api apip93.F apisorc.F fimp.F \
	  fmombal.F inelrates.F sputt.F
	if test -f apipy.so -a $(BUILDDIR) != .; then \
	  mv apipy.so $(BUILDDIR) ; \
	fi

clean-local:
	rm -f apipy.$(SO) *.F forthonf2c.h *.mod
	for i in $(DISTSOURCES); do \
	  if test -L $$i; then \
	    rm -f $$i; \
	  fi \
	done
	-cd $(BUILDDIR) && rm -f apifcn.o apip93.o apisorc.o fimp.o \
	  fmombal.o inelrates.o sputt.o api_p.F90 api_p.o apipymodule.c \
	  apipymodule.o apipy.so api.scalars cyield.mod emissivities.mod \
	  impurity_charge.mod impurity_source.mod impurity_transport.mod \
	  imslwrk.mod input.mod Makefile.api normalization_constants.mod \
	  p93dat.mod physical_constants2.mod pixels.mod radiation.mod \
	  reduced_ion_constants.mod reduced_ion_variables.mod \
	  sources_at_walls.mod

uninstall-local:
	rm -f $(prefix)/lib/apipy.$(SO)

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
