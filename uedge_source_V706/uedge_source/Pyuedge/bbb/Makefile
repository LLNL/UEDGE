##########################################################################
#
# Makefile for bbb
#   $Id: Makefile,v 7.0 2018/02/28 18:32:43 meyer8 Exp $
#
##########################################################################


include $(CONFIG)/Makefile.tmpl

OTHER_VFILES = ../../com/$(ARCH)/com.v
OTHER_SCALARS = -d ../../com/$(ARCH)/com.scalars
OTHER_DFILES = ../../com/$(ARCH)/com.d

#PASSFLAGS = TOP=../$(TOP)
PASSFLAGS = TOP=$(TOP)

PASSOTHERFLAGS = OTHER_SCALARS="$(OTHER_SCALARS)" \
        OTHER_VFILES="$(OTHER_VFILES)" OTHER_DFILES="$(OTHER_DFILES)"

all::
	@cd $(ARCH); \
	$(MAKE) -f $(CONFIG)/module_make $(PASSFLAGS) $(PASSOTHERFLAGS)

clean::
	@cd $(ARCH); \
	$(MAKE) -f $(CONFIG)/module_make $(PASSFLAGS) clean
