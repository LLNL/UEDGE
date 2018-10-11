## ##########################################################################
##
## File:    Makefile.am
##
## Purpose:   Make sure scripts also distributed.
##
## Version:   $Id: Makefile.am,v 7.0 2018/02/28 18:43:48 meyer8 Exp $
##
## ##########################################################################

dist_bin_SCRIPTS = \
    convert.py convert1.py filelists.py parallel.py \
    rdcontdt.py rdinitdt.py uedge.py uedgeutils.py \
    uefacets.py uefcifc.py ueparallel.py uewall.py

all-local: dist-local

clean-local:
	rm -f *.pyc
	for i in $(dist_bin_SCRIPTS); do	\
	 if test -L $$i; then			\
	  rm -f $$i;				\
	 fi					\
	done

dist-local:
	for i in $(dist_bin_SCRIPTS); do	\
	 if test ! -f $$i; then		\
	  ln -s $(srcdir)/$$i .;		\
	 fi					\
	done

