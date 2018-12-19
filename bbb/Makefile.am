## ############################################################
##
## File:  Makefile.am
##
## Purpose:  Run automake in top directory to make Makefile.in
##           in all subdirectories
##
## $Id: Makefile.am,v 7.1 2018/09/25 15:52:24 meyer8 Exp $
##
## UEDGE project
##
## ############################################################

BUILDDIR	= ../build
FORTHON_INCS	= @FORTHON_INCS@
SO	= so

SOURCES = \
	boundary.m \
	convert.m \
	geometry.m \
	griddubl.m \
	oderhs.m \
	odesetup.m \
	odesolve.m \
	potencur.m \
    ext_neutrals.m \
	turbulence.m

DISTSOURCES = \
	bbb.h \
	bbb.v \
	bbb_basis.v

EXTRA_DIST = \
	$(SOURCES) \
	$(DISTSOURCES)

INCLUDES = \
	-I$(srcdir) \
	-I$(top_builddir)/api \
	-I$(top_builddir)/aph \
	-I$(top_builddir)/com \
	-I$(top_builddir)/svr \
	-I$(top_srcdir) \
	-I../base \
	$(BABEL_INC)

.m.F:
	$(top_srcdir)/pyscripts/convert1.py -i $(srcdir) -o \
	  $$PWD `basename $<`

all-local: dist-local $(BUILDDIR)/bbbpy.$(SO)

../com/com.scalars:
	cd ../com && make com.scalars

$(BUILDDIR)/bbbpy.$(SO): \
  	boundary.F convert.F geometry.F griddubl.F oderhs.F \
  	odesetup.F odesolve.F potencur.F turbulence.F \
  	ext_neutrals.F bbb.v bbb_basis.v ../com/com.v
	echo Makefile CC = $(CC); \
	source $(top_srcdir)/debugenv.sh; \
	$(FORTHON_CMD) $(INCLUDES) --builddir $(BUILDDIR) \
	  --macros $(top_builddir)/com/com.v -d com -f boundary.F \
	  bbb convert.F geometry.F griddubl.F oderhs.F odesetup.F \
	  odesolve.F potencur.F turbulence.F ext_neutrals.F
	if test -f bbbpy.so -a $(BUILDDIR) != .; then \
          mv bbbpy.so $(BUILDDIR) ; \
        fi

clean-local:
	rm -f $(BUILDDIR)/bbbpy.$(SO) *.F forthonf2c.h *.mod
	for i in $(DISTSOURCES); do \
	  if test -L $$i; then \
	    rm -f $$i; \
	  fi \
	done
	-cd $(BUILDDIR) && rm -f boundary.o convert.o geometry.o \
	  griddubl.o oderhs.o odesetup.o odesolve.o ext_neutrals.o \
	  potencur.o turbulence.o aux.mod bbb_p.F90 bbb_p.o bbbpymodule.c \
	  bbbpymodule.o bbb.scalars bcond.mod bfield.mod cdv.mod cfric.mod \
	  coefeq.mod comflo.mod compla.mod comtra.mod condition_number.mod \
	  conduc.mod constraints.mod decomp.mod err_msg_out.mod fixsrc.mod \
	  global_input.mod global_vars.mod gradients.mod grid.mod \
	  ident_vars.mod ilutv.mod imprad.mod impurity_source_flux.mod \
	  indexes.mod indices_domain_dcg.mod indices_domain_dcl.mod \
	  indices_loc_glob_map.mod interp.mod jacaux.mod \
	  jacobian_clipping.mod jacobian_csc.mod jacobian_full.mod \
	  jacobian.mod jacobian_part.mod jacreorder.mod jac_work_arrays.mod \
	  locflux.mod lsode.mod Makefile.bbb math_problem_size.mod \
	  mcn_bkgd.mod mcn_dim.mod mcn_sources.mod mcn_test.mod \
	  model_choice.mod newtaux.mod nonzero_diagonals.mod npes_mpi.mod \
	  oldpla.mod opt_input.mod parallv.mod petsc_lib.mod phyvar.mod \
	  postproc.mod poten.mod preconditioning.mod rccoef.mod \
	  reduced_ion_interface.mod rhsides.mod rz_cell_info.mod \
	  save_terms.mod selec.mod solver_work_arrays.mod stat.mod \
	  temporary_work_arrays.mod timary.mod time_dep_nwt.mod \
	  turbulence_comm.mod turbulence_diagnostics.mod turbulence.mod \
	  ueint.mod uepar.mod variable_perturbation.mod volsrc.mod \
	  wkspace.mod ynorm.mod zag_output.mod

uninstall-local:
	rm -f $(prefix)/lib/bbbpy.$(SO)

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

if HAVE_PETSC
odesolve.lo: petsc-uedge.lo
odesolve.o: petsc-uedge.o
endif
