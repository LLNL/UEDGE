## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am,v 7.0 2018/02/28 18:59:56 meyer8 Exp $
##
## UEDGE project
##
## ############################################################

# Babel variables
BABEL_SHARED           = @BABEL_SHARED@
BABEL_SHARED_CONFIG    = @BABEL_SHARED_CONFIG@
BABEL_SHARED_LIBTOOL   = @BABEL_SHARED_LIBTOOL@
ICONV_LIBS      = @ICONV_LIBS@

# There is no way to use babel.make.  It has to be generated
# prior to invocation of automake, but the location of babel
# may not be known by then.  So once babel is found, we can
# "make babel.make" to determine the following variables, which
# we then cut and paste into the makefile.
PYMOD_SRCS = ue_Uedge_Module.c
PYTHONADMIN = __init__.py
PYTHONSRC = Uedge_Impl.py
# Curiously, none of the above seem to need to be known anyway?

EXTRA_DIST =					\
	Uedge_Impl.py

# For some reason one needs write privs to Uedge_Impl.py for this to
# function properly, this is a hack to fix this problem.
# jrc 16oct08: The following has to be done up one directory,
# as that is where babel.make is created along with these files
#all-local: dist-local
#	chmod u+w $(srcdir)/Uedge_Impl.py
#
#dist-local:
#	for i in $(EXTRA_DIST); do		\
#	  if test ! -f $$i; then		\
#	    ln -s $(srcdir)/$$i .;		\
#	  fi					\
#	done

#clean-local:
#	for i in $(EXTRA_DIST); do		\
#	  if test -L $$i; then			\
#	    rm -f $$i;				\
#	  fi					\
#	done

install-data-local:
	$(INSTALL) -d $(prefix)/lib/ue
	$(INSTALL) __init__.py $(prefix)/lib/ue
	$(INSTALL) Uedge_Impl.py $(prefix)/lib/ue
	$(INSTALL) Uedge_Aux.py $(prefix)/lib/ue
	if test -f Uedge.so; then \
	  $(INSTALL) Uedge.so $(prefix)/lib/ue; \
	fi

uninstall-local:
	rm -f $(prefix)/lib/ue/Uedge.so
	rm -f $(prefix)/lib/ue/__init__.py
	rm -f $(prefix)/lib/ue/Uedge_Impl.py
	rm -f $(prefix)/lib/ue/Uedge_Aux.py
