##########################################################################
#
# Makefile for com
#   $Id: Makefile,v 7.0 2018/02/28 18:32:43 meyer8 Exp $
#
##########################################################################


# TOP and CONFIG are set just to ennable running make starting in
# this dir. they are overridden by make-variables passed from above.
# in progress -- must set up all passflags as in Pyuedge/Makefile.
# perhaps better to add a file Pyuedge/Makefile.pass
TOP = ../..
CONFIG  = $(TOP)/Pyuedge/config

include $(CONFIG)/Makefile.tmpl

##PASSFLAGS = TOP=../$(TOP)
PASSFLAGS = TOP=$(TOP)
 
PASSOTHERFLAGS =

all::
	@cd $(ARCH); \
	$(MAKE) -f $(CONFIG)/module_make $(PASSFLAGS) $(PASSOTHERFLAGS)

clean::
	@cd $(ARCH); \
	$(MAKE) -f $(CONFIG)/module_make $(PASSFLAGS) clean
