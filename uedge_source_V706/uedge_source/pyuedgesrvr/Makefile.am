##
## $Id: Makefile.am,v 7.0 2018/02/28 18:59:56 meyer8 Exp $
##
##

ABS_TOP_SRCDIR	= @ABS_TOP_SRCDIR@

# Babel variables
BABEL_SHARED		= @BABEL_SHARED@
BABEL_SHARED_LIBTOOL	= @BABEL_SHARED_LIBTOOL@
BABEL_SHARED_CONFIG	= @BABEL_SHARED_CONFIG@
BABEL_SHARED_CC	        = @BABEL_SHARED_CC@
BABEL_SHARED_CFLAGS	= @BABEL_SHARED_CFLAGS@
ICONV_LIBS		= @ICONV_LIBS@

# This is so that automake does not swallow the include
@MAKEINCLUDE@ babel.make

# Must link splicer-block files into place first
SPLICER_FILES	= ue/Uedge_Impl.py

SUBDIRS		= ue

STAGE1SRCS	= $(IORSRCS) $(LAUNCHSRCS)
STAGE1OBJS	= $(STAGE1SRCS:.c=.lo)
STAGE2SRCS	= $(SKELSRCS)
STAGE2OBJS	= $(STAGE2SRCS:.c=.lo)
CC		= $(BABEL_SHARED_CC)
CFLAGS		= $(BABEL_SHARED_CFLAGS)

INCLUDES = \
	`$(BABEL_SHARED_CONFIG) --includes-py` \
	`$(BABEL_SHARED_CONFIG) --includes-xml2` \
	`$(BABEL_SHARED_CONFIG) --includes`

# This is how babel is called.
SIDL_FILES		= $(top_srcdir)/uedge.sidl $(FACETSIFC_INCDIR)/FacetsIfc.sidl $(FACETSIFC_INCDIR)/Facets0dIfc.sidl
BABEL_SHARED_FLAGS	= --server=python --strict-provenance

all: babel-stamp python-stamp libUedge1.la libUedge2.la libUedge.scl	\
  ue/Uedge.so facets/FacetsIfc.so facets/Facets0dIfc.so

ue/Uedge.so: python-stamp
facets/FacetsIfc.so: python-stamp
facets/Facets0dIfc.so: python-stamp

$(IORSRCS) $(SKELHDRS) $(SKELSRCS) $(LAUNCHSRCS): babel-stamp
$(STAGE1SRCS) $(STAGE2SRCS) setup.py: babel-stamp
facets/facets_FacetsIfc_Module.c facets/facets_Facets0dIfc_Module.c ue/ue_Uedge_Module.c: babel-stamp
babel.make: babel-stamp
babel-stamp: $(SIDL_FILES)
	@rm -f babel-temp
	@touch -f babel-temp
	if test "$(abs_top_srcdir)" != "$(abs_top_builddir)"; then \
	  for i in $(SPLICER_FILES); do \
            rm -f $$i; \
	    ln -s $(abs_top_srcdir)/pyuedgesrvr/$$i $$i; \
	  done; \
	fi
	chmod u+w $(srcdir)/ue/Uedge_Impl.py
	$(BABEL_SHARED) $(BABEL_SHARED_FLAGS) --vpath=$(srcdir) $(SIDL_FILES)
	@echo NOTICE babel.make remade.  May need to remake again.
	@mv -f babel-temp $@

python-stamp: setup.py facets/facets_FacetsIfc_Module.c facets/facets_Facets0dIfc_Module.c ue/ue_Uedge_Module.c
	( \
	  BABEL_SHARED_INCDIRS=`$(BABEL_SHARED_CONFIG) --includes | sed -e 's/-I/--include-dirs=/g'`; \
	  BABEL_SHARED_PYINCDIRS=`$(BABEL_SHARED_CONFIG) --includes-py | sed -e 's/-I/--include-dirs=/g'` ; \
	  BABEL_SHARED_LIBDIRS=`$(BABEL_SHARED_CONFIG) --libdir` ;\
	  python setup.py \
	    $$BABEL_SHARED_INCDIRS $$BABEL_SHARED_PYINCDIRS --include-dirs=. \
	    --library-dirs="$$BABEL_SHARED_LIBDIRS" \
	    build_ext --inplace && touch python-stamp \
	)

libUedge.scl:
	echo '<?xml version="1.0" ?>' > $@
	echo '<scl>' >> $@
	echo '  <library uri="'`pwd`'/libUedge1.la" scope="local" resolution="lazy" >' >> $@
	echo '    <class name="ue.Uedge" desc="ior/impl" />' >> $@
	echo '  </library>' >> $@
	echo '  <library uri="'`pwd`'/libUedge2.la" scope="local" resolution="lazy" >' >> $@
	echo '    <class name="ue.Uedge" desc="python/impl" />' >> $@
	echo '  </library>' >> $@
	echo '</scl>' >> $@

libUedge1.la: $(STAGE1OBJS)
	$(BABEL_SHARED_LIBTOOL) --tag=CC --mode=link $(CC) -rpath `pwd`/.libs -module -no-undefined  -o $@ $^ `$(BABEL_SHARED_CONFIG) --libs-c-client` $(ICONV_LIBS)
libUedge2.la : $(STAGE2OBJS)
	$(BABEL_SHARED_LIBTOOL) --tag=CC --mode=link $(CC) -rpath `pwd`/.libs -module -no-undefined -o $@ $^ `$(BABEL_SHARED_CONFIG) --libs-python-client` `$(BABEL_SHARED_CONFIG) --libs-python`

$(IORSRCS:.c=.lo) : $(IORSRCS) $(IORHDRS)
$(SKELSRCS:.c=.lo) : $(IORHDRS) $(SKELHDRS)
$(LAUNCHSRCS:.c=.lo) : $(IORHDRS)

.c.lo:
	$(BABEL_SHARED_LIBTOOL) --tag=CC --mode=compile $(CC) -c $(CFLAGS) \
		$(INCLUDES) -o $@ $<

install-data-local:
	$(INSTALL) -d $(prefix)/include
	$(INSTALL) ue_Uedge_IOR.h $(prefix)/include
	$(INSTALL) facets_FacetsIfc_IOR.h $(prefix)/include
	$(INSTALL) facets_Facets0dIfc_IOR.h $(prefix)/include
	cp -f libUedge.scl libUedge.scl.bak
	sed -e "s%`pwd`%$(prefix)/lib%g" <libUedge.scl.bak >libUedge.scl
	$(INSTALL) -d $(prefix)/lib
	$(INSTALL) -d $(prefix)/lib/facets
	$(BABEL_SHARED_LIBTOOL) --mode=install $(INSTALL) libUedge1.la $(prefix)/lib
	$(BABEL_SHARED_LIBTOOL) --mode=install $(INSTALL) libUedge2.la $(prefix)/lib
	$(INSTALL) libUedge.scl $(prefix)/lib
	$(INSTALL) facets/__init__.py $(prefix)/lib/facets
	$(INSTALL) facets/FacetsIfc.so $(prefix)/lib/facets
	$(INSTALL) facets/Facets0dIfc.so $(prefix)/lib/facets
	mv libUedge.scl.bak libUedge.scl

uninstall-local:
	$(BABEL_SHARED_LIBTOOL) --mode=uninstall /bin/rm -f $(prefix)/lib/libUedge1.la
	$(BABEL_SHARED_LIBTOOL) --mode=uninstall /bin/rm -f $(prefix)/lib/libUedge2.la
	rm -f $(prefix)/lib/libUedge.scl
	rm -f $(prefix)/lib/facets/__init__.py
	rm -f $(prefix)/lib/facets/FacetsIfc.so
	rm -f $(prefix)/lib/facets/Facets0dIfc.so

# clean-local does not remove generated sources
# dependencies should make sure those are removed
clean-local:
	rm -rf build .libs
	rm -f *.o *.lo *.la *~
	rm -f facets/*.so facets/*.o facets/*.pyc
	rm -f ue/*.so ue/*.o ue/*.pyc
	rm -f libUedge.scl

distclean-local:
	rm -f $(IORSRCS) $(SKELHDRS) $(SKELSRCS) $(LAUNCHSRCS)
	rm -f $(STAGE1SRCS) $(STAGE2SRCS) $(IORHDRS) $(PYMOD_HDRS)
	rm -f  *.scl setup.py babel.make babel-stamp
	rm -f ue/babel.make ue/*_Module.c facets/babel.make facets/*_Module.c
	rm -f facets/__init__.py ue/__init__.py

realclean: clean-local

SUFFIXES: .lo la
.PHONY: clean realclean

