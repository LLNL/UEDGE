## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am,v 7.0 2018/02/28 18:59:55 meyer8 Exp $
##
## UEDGE project
##
## ############################################################

BUILDDIR	= ../build
FORTHON_INCS	= @FORTHON_INCS@
SO	= so

SOURCES = \
	ncl_driver.m

DISTSOURCES = \
	nclass_mod.f \
	rarray_copy.f \
	rarray_sum.f \
	rarray_zero.f \
	u_erf.f \
	u_lu_backsub.f \
	u_lu_decomp.f \
	write_mod.f \
	ncl.v

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

all-local: dist-local $(BUILDDIR)/nclpy.$(SO)

../com/com.scalars :
	cd ../com ; make com.scalars


$(BUILDDIR)/nclpy.$(SO): ncl_driver.F nclass_mod.f  rarray_copy.f  rarray_sum.f \
	   rarray_zero.f  u_erf.f  u_lu_backsub.f  u_lu_decomp.f  write_mod.f \
	ncl.v ../com/com.v
	echo Makefile CC = $(CC); \
	source $(top_srcdir)/debugenv.sh; \
	$(FORTHON_CMD) $(INCLUDES) --builddir $(BUILDDIR) \
	  -d com -f nclass_mod.f ncl rarray_copy.f  rarray_sum.f \
	   rarray_zero.f  u_erf.f  u_lu_backsub.f  u_lu_decomp.f  write_mod.f ncl_driver.F
	if test -f nclpy.so -a $(BUILDDIR) != .; then \
          mv nclpy.so $(BUILDDIR); \
        fi

clean-local:
	rm -f *.F forthonf2c.h *.mod
	for i in $(DISTSOURCES); do \
	  if test -L $$i; then \
	    rm -f $$i; \
	  fi \
	done
	-cd $(BUILDDIR) && rm -f ncl_driver.o *.o

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
