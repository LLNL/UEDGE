! The routines below are for debugging of OMP implementation
subroutine WriteArrayReal(array,s,iu)
    implicit none
    real:: array(*)
    integer:: i,s,iu
    do i=1,s
        write(iu,*) array(i)
    enddo
end subroutine WriteArrayReal

subroutine WriteArrayInteger(array,s,iu)
    implicit none
    integer:: array(*)
    integer:: i,s,iu
    do i=1,s
        write(iu,*) array(i)
    enddo
end subroutine WriteArrayInteger
!JG 03/2020 Subroutine for debugging purpose: print out all the threadprivate variable.
!JG This routine was created with the OMPDebugger.py script, which also contains methods to analyze the outpout of this routine.
!JG This routine can be called within the jacobian constructor to determine which variables behave differenlty in serial vs openmp jacobian calculations
      subroutine DebugHelper(FileName)

      Use Bcond
      Use Cfric
      Use Coefeq
      Use Comflo
      Use Comgeo
      Use Compla
      Use Comtra
      Use Conduc
      Use Dim
      Use Gradients
      Use Imprad
      Use Imslwrk
      Use Indices_domain_dcl
      Use Jacobian_csc
      Use Locflux
      Use MCN_sources
      Use PNC_params
      Use Poten
      Use Rccoef
      Use Reduced_ion_interface
      Use Rhsides
      Use Save_terms
      Use Selec
      Use Time_dep_nwt
      Use Timespl
      Use Timing
      Use UEpar
      Use Wkspace
      implicit none
      integer:: iunit
      character(len = *) ::  filename
      open (newunit = iunit, file = trim(filename))
      write(iunit,*) "alfe"
      call WriteArrayReal(alfe,size(alfe),iunit)
      write(iunit,*) "alfneo"
      call WriteArrayReal(alfneo,size(alfneo),iunit)
      write(iunit,*) "bcel"
      call WriteArrayReal(bcel,size(bcel),iunit)
      write(iunit,*) "bcer"
      call WriteArrayReal(bcer,size(bcer),iunit)
      write(iunit,*) "bcil"
      call WriteArrayReal(bcil,size(bcil),iunit)
      write(iunit,*) "bcir"
      call WriteArrayReal(bcir,size(bcir),iunit)
      write(iunit,*) "betai"
      call WriteArrayReal(betai,size(betai),iunit)
      write(iunit,*) "betap"
      call WriteArrayReal(betap,size(betap),iunit)
      write(iunit,*) "cfneut"
      write(iunit,*) cfneut
      write(iunit,*) "cfneutdiv"
      write(iunit,*) cfneutdiv
      write(iunit,*) "cfvcsx"
      call WriteArrayReal(cfvcsx,size(cfvcsx),iunit)
      write(iunit,*) "cfvcsy"
      call WriteArrayReal(cfvcsy,size(cfvcsy),iunit)
      write(iunit,*) "cfvgpx"
      call WriteArrayReal(cfvgpx,size(cfvgpx),iunit)
      write(iunit,*) "cfvgpy"
      call WriteArrayReal(cfvgpy,size(cfvgpy),iunit)
      write(iunit,*) "cmneut"
      write(iunit,*) cmneut
      write(iunit,*) "cmneutdiv"
      write(iunit,*) cmneutdiv
      write(iunit,*) "coef1"
      write(iunit,*) coef1
      write(iunit,*) "coll_fe"
      call WriteArrayReal(coll_fe,size(coll_fe),iunit)
      write(iunit,*) "coll_fi"
      call WriteArrayReal(coll_fi,size(coll_fi),iunit)
      write(iunit,*) "conx"
      call WriteArrayReal(conx,size(conx),iunit)
      write(iunit,*) "conxe"
      call WriteArrayReal(conxe,size(conxe),iunit)
      write(iunit,*) "conxg"
      call WriteArrayReal(conxg,size(conxg),iunit)
      write(iunit,*) "conxge"
      call WriteArrayReal(conxge,size(conxge),iunit)
      write(iunit,*) "conxi"
      call WriteArrayReal(conxi,size(conxi),iunit)
      write(iunit,*) "cony"
      call WriteArrayReal(cony,size(cony),iunit)
      write(iunit,*) "conye"
      call WriteArrayReal(conye,size(conye),iunit)
      write(iunit,*) "conyg"
      call WriteArrayReal(conyg,size(conyg),iunit)
      write(iunit,*) "conyge"
      call WriteArrayReal(conyge,size(conyge),iunit)
      write(iunit,*) "conyi"
      call WriteArrayReal(conyi,size(conyi),iunit)
      write(iunit,*) "cs"
      write(iunit,*) cs
      write(iunit,*) "ctaue"
      call WriteArrayReal(ctaue,size(ctaue),iunit)
      write(iunit,*) "ctaui"
      call WriteArrayReal(ctaui,size(ctaui),iunit)
      write(iunit,*) "dclass_e"
      call WriteArrayReal(dclass_e,size(dclass_e),iunit)
      write(iunit,*) "dclass_i"
      call WriteArrayReal(dclass_i,size(dclass_i),iunit)
      write(iunit,*) "den"
      call WriteArrayReal(den,size(den),iunit)
      write(iunit,*) "dif2_use"
      call WriteArrayReal(dif2_use,size(dif2_use),iunit)
      write(iunit,*) "dif_use"
      call WriteArrayReal(dif_use,size(dif_use),iunit)
      write(iunit,*) "diffusivwrk"
      call WriteArrayReal(diffusivwrk,size(diffusivwrk),iunit)
      write(iunit,*) "difp_use"
      call WriteArrayReal(difp_use,size(difp_use),iunit)
      write(iunit,*) "dp1"
      write(iunit,*) dp1
      write(iunit,*) "dphi_iy1"
      call WriteArrayReal(dphi_iy1,size(dphi_iy1),iunit)
      write(iunit,*) "dtold"
      write(iunit,*) dtold
      write(iunit,*) "dtreal"
      write(iunit,*) dtreal
      write(iunit,*) "dutm_use"
      call WriteArrayReal(dutm_use,size(dutm_use),iunit)
      write(iunit,*) "dztot"
      call WriteArrayReal(dztot,size(dztot),iunit)
      write(iunit,*) "eeli"
      call WriteArrayReal(eeli,size(eeli),iunit)
      write(iunit,*) "eqp"
      call WriteArrayReal(eqp,size(eqp),iunit)
      write(iunit,*) "eqpg"
      call WriteArrayReal(eqpg,size(eqpg),iunit)
      write(iunit,*) "erliz"
      call WriteArrayReal(erliz,size(erliz),iunit)
      write(iunit,*) "erlrc"
      call WriteArrayReal(erlrc,size(erlrc),iunit)
      write(iunit,*) "eta1"
      call WriteArrayReal(eta1,size(eta1),iunit)
      write(iunit,*) "ex"
      call WriteArrayReal(ex,size(ex),iunit)
      write(iunit,*) "ey"
      call WriteArrayReal(ey,size(ey),iunit)
      write(iunit,*) "fcdif"
      write(iunit,*) fcdif
      write(iunit,*) "fdiaxlb"
      call WriteArrayReal(fdiaxlb,size(fdiaxlb),iunit)
      write(iunit,*) "fdiaxrb"
      call WriteArrayReal(fdiaxrb,size(fdiaxrb),iunit)
      write(iunit,*) "feex"
      call WriteArrayReal(feex,size(feex),iunit)
      write(iunit,*) "feexy"
      call WriteArrayReal(feexy,size(feexy),iunit)
      write(iunit,*) "feey"
      call WriteArrayReal(feey,size(feey),iunit)
      write(iunit,*) "feey4ord"
      call WriteArrayReal(feey4ord,size(feey4ord),iunit)
      write(iunit,*) "feeycbo"
      call WriteArrayReal(feeycbo,size(feeycbo),iunit)
      write(iunit,*) "fegx"
      call WriteArrayReal(fegx,size(fegx),iunit)
      write(iunit,*) "fegy"
      call WriteArrayReal(fegy,size(fegy),iunit)
      write(iunit,*) "feix"
      call WriteArrayReal(feix,size(feix),iunit)
      write(iunit,*) "feixy"
      call WriteArrayReal(feixy,size(feixy),iunit)
      write(iunit,*) "feiy"
      call WriteArrayReal(feiy,size(feiy),iunit)
      write(iunit,*) "feiy4ord"
      call WriteArrayReal(feiy4ord,size(feiy4ord),iunit)
      write(iunit,*) "feiycbo"
      call WriteArrayReal(feiycbo,size(feiycbo),iunit)
      write(iunit,*) "flngx"
      call WriteArrayReal(flngx,size(flngx),iunit)
      write(iunit,*) "flngxy"
      call WriteArrayReal(flngxy,size(flngxy),iunit)
      write(iunit,*) "flngy"
      call WriteArrayReal(flngy,size(flngy),iunit)
      write(iunit,*) "flox"
      call WriteArrayReal(flox,size(flox),iunit)
      write(iunit,*) "floxe"
      call WriteArrayReal(floxe,size(floxe),iunit)
      write(iunit,*) "floxebgt"
      call WriteArrayReal(floxebgt,size(floxebgt),iunit)
      write(iunit,*) "floxg"
      call WriteArrayReal(floxg,size(floxg),iunit)
      write(iunit,*) "floxge"
      call WriteArrayReal(floxge,size(floxge),iunit)
      write(iunit,*) "floxi"
      call WriteArrayReal(floxi,size(floxi),iunit)
      write(iunit,*) "floxibgt"
      call WriteArrayReal(floxibgt,size(floxibgt),iunit)
      write(iunit,*) "floy"
      call WriteArrayReal(floy,size(floy),iunit)
      write(iunit,*) "floye"
      call WriteArrayReal(floye,size(floye),iunit)
      write(iunit,*) "floyg"
      call WriteArrayReal(floyg,size(floyg),iunit)
      write(iunit,*) "floyge"
      call WriteArrayReal(floyge,size(floyge),iunit)
      write(iunit,*) "floyi"
      call WriteArrayReal(floyi,size(floyi),iunit)
      write(iunit,*) "fmihxpt"
      call WriteArrayReal(fmihxpt,size(fmihxpt),iunit)
      write(iunit,*) "fmity"
      call WriteArrayReal(fmity,size(fmity),iunit)
      write(iunit,*) "fmivxpt"
      call WriteArrayReal(fmivxpt,size(fmivxpt),iunit)
      write(iunit,*) "fmix"
      call WriteArrayReal(fmix,size(fmix),iunit)
      write(iunit,*) "fmixy"
      call WriteArrayReal(fmixy,size(fmixy),iunit)
      write(iunit,*) "fmiy"
      call WriteArrayReal(fmiy,size(fmiy),iunit)
      write(iunit,*) "fngx"
      call WriteArrayReal(fngx,size(fngx),iunit)
      write(iunit,*) "fngx4ord"
      call WriteArrayReal(fngx4ord,size(fngx4ord),iunit)
      write(iunit,*) "fngxy"
      call WriteArrayReal(fngxy,size(fngxy),iunit)
      write(iunit,*) "fngy"
      call WriteArrayReal(fngy,size(fngy),iunit)
      write(iunit,*) "fngy4ord"
      call WriteArrayReal(fngy4ord,size(fngy4ord),iunit)
      write(iunit,*) "fngysi"
      call WriteArrayReal(fngysi,size(fngysi),iunit)
      write(iunit,*) "fngyso"
      call WriteArrayReal(fngyso,size(fngyso),iunit)
      write(iunit,*) "fnix"
      call WriteArrayReal(fnix,size(fnix),iunit)
      write(iunit,*) "fnixcb"
      call WriteArrayReal(fnixcb,size(fnixcb),iunit)
      write(iunit,*) "fniy"
      call WriteArrayReal(fniy,size(fniy),iunit)
      write(iunit,*) "fniy4ord"
      call WriteArrayReal(fniy4ord,size(fniy4ord),iunit)
      write(iunit,*) "fniycb"
      call WriteArrayReal(fniycb,size(fniycb),iunit)
      write(iunit,*) "fniycbo"
      call WriteArrayReal(fniycbo,size(fniycbo),iunit)
      write(iunit,*) "fq2"
      call WriteArrayReal(fq2,size(fq2),iunit)
      write(iunit,*) "fq2d"
      call WriteArrayReal(fq2d,size(fq2d),iunit)
      write(iunit,*) "fqp"
      call WriteArrayReal(fqp,size(fqp),iunit)
      write(iunit,*) "fqpsatlb"
      call WriteArrayReal(fqpsatlb,size(fqpsatlb),iunit)
      write(iunit,*) "fqpsatrb"
      call WriteArrayReal(fqpsatrb,size(fqpsatrb),iunit)
      write(iunit,*) "fqx"
      call WriteArrayReal(fqx,size(fqx),iunit)
      write(iunit,*) "fqxb"
      call WriteArrayReal(fqxb,size(fqxb),iunit)
      write(iunit,*) "fqy"
      call WriteArrayReal(fqy,size(fqy),iunit)
      write(iunit,*) "fqya"
      call WriteArrayReal(fqya,size(fqya),iunit)
      write(iunit,*) "fqyae"
      call WriteArrayReal(fqyae,size(fqyae),iunit)
      write(iunit,*) "fqyai"
      call WriteArrayReal(fqyai,size(fqyai),iunit)
      write(iunit,*) "fqyao"
      call WriteArrayReal(fqyao,size(fqyao),iunit)
      write(iunit,*) "fqyb"
      call WriteArrayReal(fqyb,size(fqyb),iunit)
      write(iunit,*) "fqyd"
      call WriteArrayReal(fqyd,size(fqyd),iunit)
      write(iunit,*) "fqydt"
      call WriteArrayReal(fqydt,size(fqydt),iunit)
      write(iunit,*) "fqydti"
      call WriteArrayReal(fqydti,size(fqydti),iunit)
      write(iunit,*) "fqygp"
      call WriteArrayReal(fqygp,size(fqygp),iunit)
      write(iunit,*) "fqym"
      call WriteArrayReal(fqym,size(fqym),iunit)
      write(iunit,*) "fqymi"
      call WriteArrayReal(fqymi,size(fqymi),iunit)
      write(iunit,*) "fqyn"
      call WriteArrayReal(fqyn,size(fqyn),iunit)
      write(iunit,*) "frice"
      call WriteArrayReal(frice,size(frice),iunit)
      write(iunit,*) "frici"
      call WriteArrayReal(frici,size(frici),iunit)
      write(iunit,*) "fricnrl"
      call WriteArrayReal(fricnrl,size(fricnrl),iunit)
      write(iunit,*) "fxe"
      write(iunit,*) fxe
      write(iunit,*) "fxi"
      write(iunit,*) fxi
      write(iunit,*) "ghxpt"
      write(iunit,*) ghxpt
      write(iunit,*) "gpex"
      call WriteArrayReal(gpex,size(gpex),iunit)
      write(iunit,*) "gpey"
      call WriteArrayReal(gpey,size(gpey),iunit)
      write(iunit,*) "gpix"
      call WriteArrayReal(gpix,size(gpix),iunit)
      write(iunit,*) "gpiy"
      call WriteArrayReal(gpiy,size(gpiy),iunit)
      write(iunit,*) "gprx"
      call WriteArrayReal(gprx,size(gprx),iunit)
      write(iunit,*) "gpry"
      call WriteArrayReal(gpry,size(gpry),iunit)
      write(iunit,*) "gradp"
      call WriteArrayReal(gradp,size(gradp),iunit)
      write(iunit,*) "gradt"
      call WriteArrayReal(gradt,size(gradt),iunit)
      write(iunit,*) "gtex"
      call WriteArrayReal(gtex,size(gtex),iunit)
      write(iunit,*) "gtey"
      call WriteArrayReal(gtey,size(gtey),iunit)
      write(iunit,*) "gtix"
      call WriteArrayReal(gtix,size(gtix),iunit)
      write(iunit,*) "gtiy"
      call WriteArrayReal(gtiy,size(gtiy),iunit)
      write(iunit,*) "gvxpt"
      write(iunit,*) gvxpt
      write(iunit,*) "hcxe"
      call WriteArrayReal(hcxe,size(hcxe),iunit)
      write(iunit,*) "hcxg"
      call WriteArrayReal(hcxg,size(hcxg),iunit)
      write(iunit,*) "hcxi"
      call WriteArrayReal(hcxi,size(hcxi),iunit)
      write(iunit,*) "hcxij"
      call WriteArrayReal(hcxij,size(hcxij),iunit)
      write(iunit,*) "hcxineo"
      call WriteArrayReal(hcxineo,size(hcxineo),iunit)
      write(iunit,*) "hcxn"
      call WriteArrayReal(hcxn,size(hcxn),iunit)
      write(iunit,*) "hcye"
      call WriteArrayReal(hcye,size(hcye),iunit)
      write(iunit,*) "hcyg"
      call WriteArrayReal(hcyg,size(hcyg),iunit)
      write(iunit,*) "hcyi"
      call WriteArrayReal(hcyi,size(hcyi),iunit)
      write(iunit,*) "hcyij"
      call WriteArrayReal(hcyij,size(hcyij),iunit)
      write(iunit,*) "hcyn"
      call WriteArrayReal(hcyn,size(hcyn),iunit)
      write(iunit,*) "i1"
      write(iunit,*) i1
      write(iunit,*) "i2"
      write(iunit,*) i2
      write(iunit,*) "i2p"
      write(iunit,*) i2p
      write(iunit,*) "i3"
      write(iunit,*) i3
      write(iunit,*) "i4"
      write(iunit,*) i4
      write(iunit,*) "i5"
      write(iunit,*) i5
      write(iunit,*) "i5m"
      write(iunit,*) i5m
      write(iunit,*) "i6"
      write(iunit,*) i6
      write(iunit,*) "i7"
      write(iunit,*) i7
      write(iunit,*) "i8"
      write(iunit,*) i8
      write(iunit,*) "impradloc"
      call WriteArrayReal(impradloc,size(impradloc),iunit)
      write(iunit,*) "ispwrbcl"
      write(iunit,*) ispwrbcl
      write(iunit,*) "iwalli"
      call WriteArrayReal(iwalli,size(iwalli),iunit)
      write(iunit,*) "iwallo"
      call WriteArrayReal(iwallo,size(iwallo),iunit)
      write(iunit,*) "ixf6"
      write(iunit,*) ixf6
      write(iunit,*) "ixs1"
      write(iunit,*) ixs1
      write(iunit,*) "iyf6"
      write(iunit,*) iyf6
      write(iunit,*) "iys1"
      write(iunit,*) iys1
      write(iunit,*) "j1"
      write(iunit,*) j1
      write(iunit,*) "j1p"
      write(iunit,*) j1p
      write(iunit,*) "j2"
      write(iunit,*) j2
      write(iunit,*) "j2p"
      write(iunit,*) j2p
      write(iunit,*) "j3"
      write(iunit,*) j3
      write(iunit,*) "j4"
      write(iunit,*) j4
      write(iunit,*) "j5"
      write(iunit,*) j5
      write(iunit,*) "j5m"
      write(iunit,*) j5m
      write(iunit,*) "j5p"
      write(iunit,*) j5p
      write(iunit,*) "j6"
      write(iunit,*) j6
      write(iunit,*) "j6p"
      write(iunit,*) j6p
      write(iunit,*) "j7"
      write(iunit,*) j7
      write(iunit,*) "j8"
      write(iunit,*) j8
      write(iunit,*) "k2neo"
      call WriteArrayReal(k2neo,size(k2neo),iunit)
      write(iunit,*) "kappal"
      call WriteArrayReal(kappal,size(kappal),iunit)
      write(iunit,*) "kappar"
      call WriteArrayReal(kappar,size(kappar),iunit)
      write(iunit,*) "kincorlb"
      call WriteArrayReal(kincorlb,size(kincorlb),iunit)
      write(iunit,*) "kincorrb"
      call WriteArrayReal(kincorrb,size(kincorrb),iunit)
      write(iunit,*) "ktneo"
      call WriteArrayReal(ktneo,size(ktneo),iunit)
      write(iunit,*) "kxbohm"
      call WriteArrayReal(kxbohm,size(kxbohm),iunit)
      write(iunit,*) "kxe_use"
      call WriteArrayReal(kxe_use,size(kxe_use),iunit)
      write(iunit,*) "kxi_use"
      call WriteArrayReal(kxi_use,size(kxi_use),iunit)
      write(iunit,*) "kybohm"
      call WriteArrayReal(kybohm,size(kybohm),iunit)
      write(iunit,*) "kye_use"
      call WriteArrayReal(kye_use,size(kye_use),iunit)
      write(iunit,*) "kyi_use"
      call WriteArrayReal(kyi_use,size(kyi_use),iunit)
      write(iunit,*) "lng"
      call WriteArrayReal(lng,size(lng),iunit)
      write(iunit,*) "loglambda"
      call WriteArrayReal(loglambda,size(loglambda),iunit)
      write(iunit,*) "msor"
      call WriteArrayReal(msor,size(msor),iunit)
      write(iunit,*) "msorold"
      call WriteArrayReal(msorold,size(msorold),iunit)
      write(iunit,*) "msorxr"
      call WriteArrayReal(msorxr,size(msorxr),iunit)
      write(iunit,*) "msorxrold"
      call WriteArrayReal(msorxrold,size(msorxrold),iunit)
      write(iunit,*) "na"
      call WriteArrayReal(na,size(na),iunit)
      write(iunit,*) "ncrhs"
      write(iunit,*) ncrhs
      write(iunit,*) "ne"
      call WriteArrayReal(ne,size(ne),iunit)
      write(iunit,*) "netap"
      call WriteArrayReal(netap,size(netap),iunit)
      write(iunit,*) "ney0"
      call WriteArrayReal(ney0,size(ney0),iunit)
      write(iunit,*) "ney1"
      call WriteArrayReal(ney1,size(ney1),iunit)
      write(iunit,*) "nfsp"
      write(iunit,*) nfsp
      write(iunit,*) "ng"
      call WriteArrayReal(ng,size(ng),iunit)
      write(iunit,*) "ngy0"
      call WriteArrayReal(ngy0,size(ngy0),iunit)
      write(iunit,*) "ngy1"
      call WriteArrayReal(ngy1,size(ngy1),iunit)
      write(iunit,*) "ni"
      call WriteArrayReal(ni,size(ni),iunit)
      write(iunit,*) "nit"
      call WriteArrayReal(nit,size(nit),iunit)
      write(iunit,*) "nity0"
      call WriteArrayReal(nity0,size(nity0),iunit)
      write(iunit,*) "nity1"
      call WriteArrayReal(nity1,size(nity1),iunit)
      write(iunit,*) "nixpt"
      call WriteArrayReal(nixpt,size(nixpt),iunit)
      write(iunit,*) "niy0"
      call WriteArrayReal(niy0,size(niy0),iunit)
      write(iunit,*) "niy0s"
      call WriteArrayReal(niy0s,size(niy0s),iunit)
      write(iunit,*) "niy1"
      call WriteArrayReal(niy1,size(niy1),iunit)
      write(iunit,*) "niy1s"
      call WriteArrayReal(niy1s,size(niy1s),iunit)
      write(iunit,*) "nm"
      call WriteArrayReal(nm,size(nm),iunit)
      write(iunit,*) "nratio"
      call WriteArrayReal(nratio,size(nratio),iunit)
      write(iunit,*) "ntau"
      call WriteArrayReal(ntau,size(ntau),iunit)
      write(iunit,*) "nucx"
      call WriteArrayReal(nucx,size(nucx),iunit)
      write(iunit,*) "nucxi"
      call WriteArrayReal(nucxi,size(nucxi),iunit)
      write(iunit,*) "nuelg"
      call WriteArrayReal(nuelg,size(nuelg),iunit)
      write(iunit,*) "nueli"
      call WriteArrayReal(nueli,size(nueli),iunit)
      write(iunit,*) "nuii"
      call WriteArrayReal(nuii,size(nuii),iunit)
      write(iunit,*) "nuiistar"
      call WriteArrayReal(nuiistar,size(nuiistar),iunit)
      write(iunit,*) "nuix"
      call WriteArrayReal(nuix,size(nuix),iunit)
      write(iunit,*) "nuiz"
      call WriteArrayReal(nuiz,size(nuiz),iunit)
      write(iunit,*) "nurc"
      call WriteArrayReal(nurc,size(nurc),iunit)
      write(iunit,*) "nuvl"
      call WriteArrayReal(nuvl,size(nuvl),iunit)
      write(iunit,*) "nz2"
      call WriteArrayReal(nz2,size(nz2),iunit)
      write(iunit,*) "nzloc"
      call WriteArrayReal(nzloc,size(nzloc),iunit)
      write(iunit,*) "openbox"
      write(iunit,*) openbox
      write(iunit,*) "parvis"
      call WriteArrayReal(parvis,size(parvis),iunit)
      write(iunit,*) "pg"
      call WriteArrayReal(pg,size(pg),iunit)
      write(iunit,*) "pgy0"
      call WriteArrayReal(pgy0,size(pgy0),iunit)
      write(iunit,*) "pgy1"
      call WriteArrayReal(pgy1,size(pgy1),iunit)
      write(iunit,*) "phi"
      call WriteArrayReal(phi,size(phi),iunit)
      write(iunit,*) "phiv"
      call WriteArrayReal(phiv,size(phiv),iunit)
      write(iunit,*) "phiy0"
      call WriteArrayReal(phiy0,size(phiy0),iunit)
      write(iunit,*) "phiy0s"
      call WriteArrayReal(phiy0s,size(phiy0s),iunit)
      write(iunit,*) "phiy1"
      call WriteArrayReal(phiy1,size(phiy1),iunit)
      write(iunit,*) "phiy1s"
      call WriteArrayReal(phiy1s,size(phiy1s),iunit)
      write(iunit,*) "pr"
      call WriteArrayReal(pr,size(pr),iunit)
      write(iunit,*) "prad"
      call WriteArrayReal(prad,size(prad),iunit)
      write(iunit,*) "pradc"
      call WriteArrayReal(pradc,size(pradc),iunit)
      write(iunit,*) "pradcff"
      call WriteArrayReal(pradcff,size(pradcff),iunit)
      write(iunit,*) "pradhyd"
      call WriteArrayReal(pradhyd,size(pradhyd),iunit)
      write(iunit,*) "pradz"
      call WriteArrayReal(pradz,size(pradz),iunit)
      write(iunit,*) "pradzc"
      call WriteArrayReal(pradzc,size(pradzc),iunit)
      write(iunit,*) "pre"
      call WriteArrayReal(pre,size(pre),iunit)
      write(iunit,*) "prev"
      call WriteArrayReal(prev,size(prev),iunit)
      write(iunit,*) "pri"
      call WriteArrayReal(pri,size(pri),iunit)
      write(iunit,*) "priv"
      call WriteArrayReal(priv,size(priv),iunit)
      write(iunit,*) "priy0"
      call WriteArrayReal(priy0,size(priy0),iunit)
      write(iunit,*) "priy1"
      call WriteArrayReal(priy1,size(priy1),iunit)
      write(iunit,*) "prtv"
      call WriteArrayReal(prtv,size(prtv),iunit)
      write(iunit,*) "psor"
      call WriteArrayReal(psor,size(psor),iunit)
      write(iunit,*) "psor_tmpov"
      call WriteArrayReal(psor_tmpov,size(psor_tmpov),iunit)
      write(iunit,*) "psorbgg"
      call WriteArrayReal(psorbgg,size(psorbgg),iunit)
      write(iunit,*) "psorbgz"
      call WriteArrayReal(psorbgz,size(psorbgz),iunit)
      write(iunit,*) "psorc"
      call WriteArrayReal(psorc,size(psorc),iunit)
      write(iunit,*) "psorcxg"
      call WriteArrayReal(psorcxg,size(psorcxg),iunit)
      write(iunit,*) "psorcxgc"
      call WriteArrayReal(psorcxgc,size(psorcxgc),iunit)
      write(iunit,*) "psordis"
      call WriteArrayReal(psordis,size(psordis),iunit)
      write(iunit,*) "psorg"
      call WriteArrayReal(psorg,size(psorg),iunit)
      write(iunit,*) "psorgc"
      call WriteArrayReal(psorgc,size(psorgc),iunit)
      write(iunit,*) "psori"
      call WriteArrayReal(psori,size(psori),iunit)
      write(iunit,*) "psorold"
      call WriteArrayReal(psorold,size(psorold),iunit)
      write(iunit,*) "psorrg"
      call WriteArrayReal(psorrg,size(psorrg),iunit)
      write(iunit,*) "psorrgc"
      call WriteArrayReal(psorrgc,size(psorrgc),iunit)
      write(iunit,*) "psorxr"
      call WriteArrayReal(psorxr,size(psorxr),iunit)
      write(iunit,*) "psorxrc"
      call WriteArrayReal(psorxrc,size(psorxrc),iunit)
      write(iunit,*) "psorxrold"
      call WriteArrayReal(psorxrold,size(psorxrold),iunit)
      write(iunit,*) "pwrebkg"
      call WriteArrayReal(pwrebkg,size(pwrebkg),iunit)
      write(iunit,*) "pwribkg"
      call WriteArrayReal(pwribkg,size(pwribkg),iunit)
      write(iunit,*) "pwrze"
      call WriteArrayReal(pwrze,size(pwrze),iunit)
      write(iunit,*) "pwrzec"
      call WriteArrayReal(pwrzec,size(pwrzec),iunit)
      write(iunit,*) "q2cd"
      call WriteArrayReal(q2cd,size(q2cd),iunit)
      write(iunit,*) "qipar"
      call WriteArrayReal(qipar,size(qipar),iunit)
      write(iunit,*) "resco"
      call WriteArrayReal(resco,size(resco),iunit)
      write(iunit,*) "resee"
      call WriteArrayReal(resee,size(resee),iunit)
      write(iunit,*) "reseg"
      call WriteArrayReal(reseg,size(reseg),iunit)
      write(iunit,*) "resei"
      call WriteArrayReal(resei,size(resei),iunit)
      write(iunit,*) "resmo"
      call WriteArrayReal(resmo,size(resmo),iunit)
      write(iunit,*) "resng"
      call WriteArrayReal(resng,size(resng),iunit)
      write(iunit,*) "resphi"
      call WriteArrayReal(resphi,size(resphi),iunit)
      write(iunit,*) "rtau"
      call WriteArrayReal(rtau,size(rtau),iunit)
      write(iunit,*) "rtaue"
      call WriteArrayReal(rtaue,size(rtaue),iunit)
      write(iunit,*) "rtaux"
      call WriteArrayReal(rtaux,size(rtaux),iunit)
      write(iunit,*) "rtauy"
      call WriteArrayReal(rtauy,size(rtauy),iunit)
      write(iunit,*) "seec"
      call WriteArrayReal(seec,size(seec),iunit)
      write(iunit,*) "seev"
      call WriteArrayReal(seev,size(seev),iunit)
      write(iunit,*) "seg_ue"
      call WriteArrayReal(seg_ue,size(seg_ue),iunit)
      write(iunit,*) "seic"
      call WriteArrayReal(seic,size(seic),iunit)
      write(iunit,*) "seiv"
      call WriteArrayReal(seiv,size(seiv),iunit)
      write(iunit,*) "smoc"
      call WriteArrayReal(smoc,size(smoc),iunit)
      write(iunit,*) "smov"
      call WriteArrayReal(smov,size(smov),iunit)
      write(iunit,*) "sng_ue"
      call WriteArrayReal(sng_ue,size(sng_ue),iunit)
      write(iunit,*) "snic"
      call WriteArrayReal(snic,size(snic),iunit)
      write(iunit,*) "sniv"
      call WriteArrayReal(sniv,size(sniv),iunit)
      write(iunit,*) "sputflxlb"
      call WriteArrayReal(sputflxlb,size(sputflxlb),iunit)
      write(iunit,*) "sputflxpf"
      call WriteArrayReal(sputflxpf,size(sputflxpf),iunit)
      write(iunit,*) "sputflxrb"
      call WriteArrayReal(sputflxrb,size(sputflxrb),iunit)
      write(iunit,*) "sputflxw"
      call WriteArrayReal(sputflxw,size(sputflxw),iunit)
      write(iunit,*) "sxyxpt"
      write(iunit,*) sxyxpt
      write(iunit,*) "te"
      call WriteArrayReal(te,size(te),iunit)
      write(iunit,*) "tempa"
      call WriteArrayReal(tempa,size(tempa),iunit)
      write(iunit,*) "tev"
      call WriteArrayReal(tev,size(tev),iunit)
      write(iunit,*) "tey0"
      call WriteArrayReal(tey0,size(tey0),iunit)
      write(iunit,*) "tey1"
      call WriteArrayReal(tey1,size(tey1),iunit)
      write(iunit,*) "tg"
      call WriteArrayReal(tg,size(tg),iunit)
      write(iunit,*) "tgy0"
      call WriteArrayReal(tgy0,size(tgy0),iunit)
      write(iunit,*) "tgy1"
      call WriteArrayReal(tgy1,size(tgy1),iunit)
      write(iunit,*) "ti"
      call WriteArrayReal(ti,size(ti),iunit)
      write(iunit,*) "tiv"
      call WriteArrayReal(tiv,size(tiv),iunit)
      write(iunit,*) "tiy0"
      call WriteArrayReal(tiy0,size(tiy0),iunit)
      write(iunit,*) "tiy0s"
      call WriteArrayReal(tiy0s,size(tiy0s),iunit)
      write(iunit,*) "tiy1"
      call WriteArrayReal(tiy1,size(tiy1),iunit)
      write(iunit,*) "tiy1s"
      call WriteArrayReal(tiy1s,size(tiy1s),iunit)
      write(iunit,*) "totb2val"
      write(iunit,*) totb2val
      write(iunit,*) "totfeexl"
      call WriteArrayReal(totfeexl,size(totfeexl),iunit)
      write(iunit,*) "totfeexr"
      call WriteArrayReal(totfeexr,size(totfeexr),iunit)
      write(iunit,*) "totfeixl"
      call WriteArrayReal(totfeixl,size(totfeixl),iunit)
      write(iunit,*) "totfeixr"
      call WriteArrayReal(totfeixr,size(totfeixr),iunit)
      write(iunit,*) "travis"
      call WriteArrayReal(travis,size(travis),iunit)
      write(iunit,*) "trax_use"
      call WriteArrayReal(trax_use,size(trax_use),iunit)
      write(iunit,*) "tray_use"
      call WriteArrayReal(tray_use,size(tray_use),iunit)
      write(iunit,*) "ttimpfe"
      write(iunit,*) ttimpfe
      write(iunit,*) "ttimpjf"
      write(iunit,*) ttimpjf
      write(iunit,*) "ttngfd2"
      write(iunit,*) ttngfd2
      write(iunit,*) "ttngfxy"
      write(iunit,*) ttngfxy
      write(iunit,*) "ttngxlog"
      write(iunit,*) ttngxlog
      write(iunit,*) "ttngylog"
      write(iunit,*) ttngylog
      write(iunit,*) "ttnpg"
      write(iunit,*) ttnpg
      write(iunit,*) "ttotfe"
      write(iunit,*) ttotfe
      write(iunit,*) "ttotjf"
      write(iunit,*) ttotjf
      write(iunit,*) "up"
      call WriteArrayReal(up,size(up),iunit)
      write(iunit,*) "upe"
      call WriteArrayReal(upe,size(upe),iunit)
      write(iunit,*) "upi"
      call WriteArrayReal(upi,size(upi),iunit)
      write(iunit,*) "upxpt"
      call WriteArrayReal(upxpt,size(upxpt),iunit)
      write(iunit,*) "uu"
      call WriteArrayReal(uu,size(uu),iunit)
      write(iunit,*) "uug"
      call WriteArrayReal(uug,size(uug),iunit)
      write(iunit,*) "uup"
      call WriteArrayReal(uup,size(uup),iunit)
      write(iunit,*) "uz"
      call WriteArrayReal(uz,size(uz),iunit)
      write(iunit,*) "v2"
      call WriteArrayReal(v2,size(v2),iunit)
      write(iunit,*) "v2cb"
      call WriteArrayReal(v2cb,size(v2cb),iunit)
      write(iunit,*) "v2cd"
      call WriteArrayReal(v2cd,size(v2cd),iunit)
      write(iunit,*) "v2ce"
      call WriteArrayReal(v2ce,size(v2ce),iunit)
      write(iunit,*) "v2dd"
      call WriteArrayReal(v2dd,size(v2dd),iunit)
      write(iunit,*) "v2rd"
      call WriteArrayReal(v2rd,size(v2rd),iunit)
      write(iunit,*) "v2xgp"
      call WriteArrayReal(v2xgp,size(v2xgp),iunit)
      write(iunit,*) "ve2cb"
      call WriteArrayReal(ve2cb,size(ve2cb),iunit)
      write(iunit,*) "ve2cd"
      call WriteArrayReal(ve2cd,size(ve2cd),iunit)
      write(iunit,*) "vex"
      call WriteArrayReal(vex,size(vex),iunit)
      write(iunit,*) "vey"
      call WriteArrayReal(vey,size(vey),iunit)
      write(iunit,*) "veycb"
      call WriteArrayReal(veycb,size(veycb),iunit)
      write(iunit,*) "veycp"
      call WriteArrayReal(veycp,size(veycp),iunit)
      write(iunit,*) "visvol_q"
      call WriteArrayReal(visvol_q,size(visvol_q),iunit)
      write(iunit,*) "visvol_v"
      call WriteArrayReal(visvol_v,size(visvol_v),iunit)
      write(iunit,*) "visx"
      call WriteArrayReal(visx,size(visx),iunit)
      write(iunit,*) "visxneo"
      call WriteArrayReal(visxneo,size(visxneo),iunit)
      write(iunit,*) "visy"
      call WriteArrayReal(visy,size(visy),iunit)
      write(iunit,*) "visyxpt"
      call WriteArrayReal(visyxpt,size(visyxpt),iunit)
      write(iunit,*) "vsoree"
      call WriteArrayReal(vsoree,size(vsoree),iunit)
      write(iunit,*) "vsoreec"
      call WriteArrayReal(vsoreec,size(vsoreec),iunit)
      write(iunit,*) "vy"
      call WriteArrayReal(vy,size(vy),iunit)
      write(iunit,*) "vy_cft"
      call WriteArrayReal(vy_cft,size(vy_cft),iunit)
      write(iunit,*) "vy_use"
      call WriteArrayReal(vy_use,size(vy_use),iunit)
      write(iunit,*) "vyavis"
      call WriteArrayReal(vyavis,size(vyavis),iunit)
      write(iunit,*) "vycb"
      call WriteArrayReal(vycb,size(vycb),iunit)
      write(iunit,*) "vyce"
      call WriteArrayReal(vyce,size(vyce),iunit)
      write(iunit,*) "vycf"
      call WriteArrayReal(vycf,size(vycf),iunit)
      write(iunit,*) "vycp"
      call WriteArrayReal(vycp,size(vycp),iunit)
      write(iunit,*) "vycr"
      call WriteArrayReal(vycr,size(vycr),iunit)
      write(iunit,*) "vydd"
      call WriteArrayReal(vydd,size(vydd),iunit)
      write(iunit,*) "vyg"
      call WriteArrayReal(vyg,size(vyg),iunit)
      write(iunit,*) "vygp"
      call WriteArrayReal(vygp,size(vygp),iunit)
      write(iunit,*) "vygtan"
      call WriteArrayReal(vygtan,size(vygtan),iunit)
      write(iunit,*) "vyhxpt"
      call WriteArrayReal(vyhxpt,size(vyhxpt),iunit)
      write(iunit,*) "vyrd"
      call WriteArrayReal(vyrd,size(vyrd),iunit)
      write(iunit,*) "vytan"
      call WriteArrayReal(vytan,size(vytan),iunit)
      write(iunit,*) "vyte_cft"
      call WriteArrayReal(vyte_cft,size(vyte_cft),iunit)
      write(iunit,*) "vyti_cft"
      call WriteArrayReal(vyti_cft,size(vyti_cft),iunit)
      write(iunit,*) "vyvxpt"
      call WriteArrayReal(vyvxpt,size(vyvxpt),iunit)
      write(iunit,*) "w"
      call WriteArrayReal(w,size(w),iunit)
      write(iunit,*) "w0"
      call WriteArrayReal(w0,size(w0),iunit)
      write(iunit,*) "w1"
      call WriteArrayReal(w1,size(w1),iunit)
      write(iunit,*) "w2"
      call WriteArrayReal(w2,size(w2),iunit)
      write(iunit,*) "w3"
      call WriteArrayReal(w3,size(w3),iunit)
      write(iunit,*) "wjdote"
      call WriteArrayReal(wjdote,size(wjdote),iunit)
      write(iunit,*) "wvh"
      call WriteArrayReal(wvh,size(wvh),iunit)
      write(iunit,*) "xcnearlb"
      write(iunit,*) xcnearlb
      write(iunit,*) "xcnearrb"
      write(iunit,*) xcnearrb
      write(iunit,*) "yld_carbi"
      call WriteArrayReal(yld_carbi,size(yld_carbi),iunit)
      write(iunit,*) "yld_carbo"
      call WriteArrayReal(yld_carbo,size(yld_carbo),iunit)
      write(iunit,*) "yldot_pert"
      call WriteArrayReal(yldot_pert,size(yldot_pert),iunit)
      write(iunit,*) "zcoef"
      write(iunit,*) zcoef
      write(iunit,*) "zeff"
      call WriteArrayReal(zeff,size(zeff),iunit)
      write(iunit,*) "znot"
      call WriteArrayReal(znot,size(znot),iunit)
      close(iunit)
      end subroutine DebugHelper


#ifdef _OPENMP
      subroutine OmpCopyPointerbbb2
use OmpCopybbb
!call OmpCopyPointerctaue
!call OmpCopyPointerctaui
!call OmpCopyPointercfvcsx
!call OmpCopyPointercfvcsy
!call OmpCopyPointercfvgpx
!call OmpCopyPointercfvgpy
!call OmpCopyPointeriwalli
!call OmpCopyPointeriwallo
!call OmpCopyPointerfngysi
!call OmpCopyPointerfngyso
!call OmpCopyPointeryld_carbi
!call OmpCopyPointeryld_carbo
!call OmpCopyPointerfqpsatlb
!call OmpCopyPointerfqpsatrb
!call OmpCopyPointertotfeexl
!call OmpCopyPointertotfeexr
!call OmpCopyPointertotfeixl
!call OmpCopyPointertotfeixr
!call OmpCopyPointersputflxlb
!call OmpCopyPointersputflxrb
!call OmpCopyPointersputflxw
!call OmpCopyPointersputflxpf
call OmpCopyPointerparvis
call OmpCopyPointertravis
!call OmpCopyPointerdiffusivwrk
!call OmpCopyPointercoll_fe
!call OmpCopyPointercoll_fi
!call OmpCopyPointerni
!call OmpCopyPointernm
!call OmpCopyPointernz2
!call OmpCopyPointeruu
!call OmpCopyPointeruup


!call OmpCopyPointerup
!call OmpCopyPointerupi
!call OmpCopyPointeruz
!
!
!call OmpCopyPointerv2
!call OmpCopyPointerv2xgp
!call OmpCopyPointerv2ce
!call OmpCopyPointerv2cb
!call OmpCopyPointerve2cb
!call OmpCopyPointerv2cd
!call OmpCopyPointerve2cd
!call OmpCopyPointerq2cd
!call OmpCopyPointerv2rd
!call OmpCopyPointerv2dd
!call OmpCopyPointervy
!call OmpCopyPointervygp
!call OmpCopyPointervytan
!call OmpCopyPointervygtan
!call OmpCopyPointervyce
!call OmpCopyPointervycb
!call OmpCopyPointerveycb
!call OmpCopyPointervycp
!call OmpCopyPointerveycp
!call OmpCopyPointervyrd
!call OmpCopyPointervydd
!call OmpCopyPointervyavis
!call OmpCopyPointervex
!call OmpCopyPointerupe
!call OmpCopyPointervey
!
!!!!!!!!!!!!!!!!!!!!!!
!call OmpCopyPointervycf
!call OmpCopyPointervycr
!call OmpCopyPointerte
!call OmpCopyPointerti
!call OmpCopyPointerng
!call OmpCopyPointerlng
!call OmpCopyPointeruug
!call OmpCopyPointervyg
!call OmpCopyPointertg
!call OmpCopyPointertev
!call OmpCopyPointertiv
!call OmpCopyPointerniy0
!call OmpCopyPointerniy1
!call OmpCopyPointerniy0s
!call OmpCopyPointerniy1s
!call OmpCopyPointerney0
!call OmpCopyPointerney1
!call OmpCopyPointernity0
!call OmpCopyPointernity1
!call OmpCopyPointertey0
!call OmpCopyPointertey1
!call OmpCopyPointertiy0
!call OmpCopyPointertiy1
!call OmpCopyPointertiy0s
!call OmpCopyPointertiy1s
!call OmpCopyPointertgy0
!
!!call OmpCopyPointertgy1
!
!write(*,*) '#Master::: Calling routine to copy variable ngy0'
!
!call OmpCopyPointerngy0
!
!write(*,*) '#Master::: Calling routine to copy variable ngy1'
!
!call OmpCopyPointerngy1
!
!write(*,*) '#Master::: Calling routine to copy variable pgy0'
!
!call OmpCopyPointerpgy0
!
!write(*,*) '#Master::: Calling routine to copy variable pgy1'
!
!call OmpCopyPointerpgy1
!
!write(*,*) '#Master::: Calling routine to copy variable pg'
!
!call OmpCopyPointerpg
!
!write(*,*) '#Master::: Calling routine to copy variable phiy0'
!
!call OmpCopyPointerphiy0
!
!write(*,*) '#Master::: Calling routine to copy variable phiy1'
!
!call OmpCopyPointerphiy1
!
!write(*,*) '#Master::: Calling routine to copy variable phiy0s'
!
!call OmpCopyPointerphiy0s
!
!write(*,*) '#Master::: Calling routine to copy variable phiy1s'
!
!call OmpCopyPointerphiy1s
!
!write(*,*) '#Master::: Calling routine to copy variable pr'
!
!call OmpCopyPointerpr
!
!write(*,*) '#Master::: Calling routine to copy variable prev'
!
!call OmpCopyPointerprev
!
!write(*,*) '#Master::: Calling routine to copy variable prtv'
!
!call OmpCopyPointerprtv
!
!write(*,*) '#Master::: Calling routine to copy variable pri'
!
!call OmpCopyPointerpri
!
!write(*,*) '#Master::: Calling routine to copy variable priv'
!
!call OmpCopyPointerpriv
!
!write(*,*) '#Master::: Calling routine to copy variable priy0'
!
!call OmpCopyPointerpriy0
!
!write(*,*) '#Master::: Calling routine to copy variable priy1'
!
!call OmpCopyPointerpriy1
!
!write(*,*) '#Master::: Calling routine to copy variable pre'
!
!call OmpCopyPointerpre
!
!write(*,*) '#Master::: Calling routine to copy variable ne'
!
!call OmpCopyPointerne
!
!write(*,*) '#Master::: Calling routine to copy variable nit'
!
!call OmpCopyPointernit
!
!write(*,*) '#Master::: Calling routine to copy variable phi'
!
!call OmpCopyPointerphi
!
!write(*,*) '#Master::: Calling routine to copy variable phiv'
!
!call OmpCopyPointerphiv
!
!write(*,*) '#Master::: Calling routine to copy variable zeff'
!
!call OmpCopyPointerzeff
!
!write(*,*) '#Master::: Calling routine to copy variable loglambda'
!
!call OmpCopyPointerloglambda
!
!write(*,*) '#Master::: Calling routine to copy variable netap'
!
!call OmpCopyPointernetap
!call OmpCopyPointerznot
!call OmpCopyPointerupxpt
!call OmpCopyPointernixpt
!call OmpCopyPointervisyxpt
!call OmpCopyPointervyhxpt
!call OmpCopyPointervyvxpt
!call OmpCopyPointerfmihxpt
!call OmpCopyPointerfmivxpt
!call OmpCopyPointerrtaux
!call OmpCopyPointerrtauy
!call OmpCopyPointerrtau
!call OmpCopyPointerbetap
!call OmpCopyPointerfqp
!call OmpCopyPointerfq2
!call OmpCopyPointerfqx
!call OmpCopyPointerfqxb
!call OmpCopyPointerfdiaxlb
!call OmpCopyPointerfdiaxrb
!call OmpCopyPointerfloxebgt
!call OmpCopyPointerfloxibgt
!call OmpCopyPointerfqy
!
!write(*,*) '#Master::: Calling routine to copy variable fqyb'
!
!call OmpCopyPointerfqyb
!
!write(*,*) '#Master::: Calling routine to copy variable fqyn'
!
!call OmpCopyPointerfqyn
!
!write(*,*) '#Master::: Calling routine to copy variable fqym'
!
!call OmpCopyPointerfqym
!
!write(*,*) '#Master::: Calling routine to copy variable fqymi'
!
!call OmpCopyPointerfqymi
!
!write(*,*) '#Master::: Calling routine to copy variable fqya'
!
!call OmpCopyPointerfqya
!
!write(*,*) '#Master::: Calling routine to copy variable fqydt'
!
!call OmpCopyPointerfqydt
!
!write(*,*) '#Master::: Calling routine to copy variable fqydti'
!
!call OmpCopyPointerfqydti
!
!write(*,*) '#Master::: Calling routine to copy variable fqyao'
!
!call OmpCopyPointerfqyao
!
!write(*,*) '#Master::: Calling routine to copy variable fqyae'
!
!call OmpCopyPointerfqyae
!
!write(*,*) '#Master::: Calling routine to copy variable fqyai'
!
!call OmpCopyPointerfqyai
!
!write(*,*) '#Master::: Calling routine to copy variable fqyd'
!
!call OmpCopyPointerfqyd
!
!write(*,*) '#Master::: Calling routine to copy variable fqygp'
!
!call OmpCopyPointerfqygp
!
!write(*,*) '#Master::: Calling routine to copy variable fq2d'
!
!call OmpCopyPointerfq2d
!
!write(*,*) '#Master::: Calling routine to copy variable fnix'
!
!call OmpCopyPointerfnix
!call OmpCopyPointerfnixcb
!call OmpCopyPointerfniy
!call OmpCopyPointerfniy4ord
!call OmpCopyPointerfniycb
!call OmpCopyPointerfmix
!call OmpCopyPointerfmiy
!call OmpCopyPointerfmixy
!call OmpCopyPointerfmity
!call OmpCopyPointerfeex
!call OmpCopyPointerfeey
!call OmpCopyPointerfeexy
!call OmpCopyPointerfeey4ord
!call OmpCopyPointerfeix
!call OmpCopyPointerfeiy
!call OmpCopyPointerfegx
!call OmpCopyPointerfegy
!call OmpCopyPointerqipar
!call OmpCopyPointerfniycbo
!call OmpCopyPointerfeiycbo
!call OmpCopyPointerfeeycbo
!call OmpCopyPointerfeixy
!call OmpCopyPointerfeiy4ord
!call OmpCopyPointerfngx
!call OmpCopyPointerfngx4ord
!call OmpCopyPointerflngx
!call OmpCopyPointerfngy
!call OmpCopyPointerfngy4ord
!call OmpCopyPointerflngy
!call OmpCopyPointerfngxy
!call OmpCopyPointerflngxy
!call OmpCopyPointerbcel
!call OmpCopyPointerbcer
!call OmpCopyPointerbcil
!call OmpCopyPointerbcir
!call OmpCopyPointerkappal
!call OmpCopyPointerkappar
!call OmpCopyPointerdphi_iy1
!call OmpCopyPointerkincorlb
!call OmpCopyPointerkincorrb
!call OmpCopyPointerex
!call OmpCopyPointerey
!call OmpCopyPointergpix
!call OmpCopyPointergpiy
!call OmpCopyPointergpex
!call OmpCopyPointergpey
!call OmpCopyPointergprx
!call OmpCopyPointergpry
!call OmpCopyPointergtex
!call OmpCopyPointergtey
!!
!!write(*,*) '#Master::: Calling routine to copy variable gtix'
!!
!!call OmpCopyPointergtix
!!
!!write(*,*) '#Master::: Calling routine to copy variable gtiy'
!!
!!call OmpCopyPointergtiy
!!
!!write(*,*) '#Master::: Calling routine to copy variable frice'
!!
!!call OmpCopyPointerfrice
!!call OmpCopyPointerfrici
!!call OmpCopyPointerfricnrl
!!call OmpCopyPointeralfe
!!call OmpCopyPointerbetai
!!call OmpCopyPointerw
!!call OmpCopyPointerw0
!!call OmpCopyPointerw1
!!
!!write(*,*) '#Master::: Calling routine to copy variable w2'
!!
!!call OmpCopyPointerw2
!!
!!write(*,*) '#Master::: Calling routine to copy variable w3'
!!
!!call OmpCopyPointerw3
!!
!!write(*,*) '#Master::: Calling routine to copy variable wvh'
!!
!!call OmpCopyPointerwvh
!!
!!write(*,*) '#Master::: Calling routine to copy variable flox'
!!
!!call OmpCopyPointerflox
!!
!!write(*,*) '#Master::: Calling routine to copy variable floy'
!!
!!call OmpCopyPointerfloy
!!
!!write(*,*) '#Master::: Calling routine to copy variable conx'
!!
!!call OmpCopyPointerconx
!!
!!write(*,*) '#Master::: Calling routine to copy variable cony'
!!
!!call OmpCopyPointercony
!!
!!write(*,*) '#Master::: Calling routine to copy variable floxe'
!!
!!call OmpCopyPointerfloxe
!!
!!write(*,*) '#Master::: Calling routine to copy variable floye'
!!
!!call OmpCopyPointerfloye
!!
!!write(*,*) '#Master::: Calling routine to copy variable floxi'
!!
!!call OmpCopyPointerfloxi
!!
!!write(*,*) '#Master::: Calling routine to copy variable floyi'
!!
!!call OmpCopyPointerfloyi
!!
!!write(*,*) '#Master::: Calling routine to copy variable floxg'
!!
!!call OmpCopyPointerfloxg
!!
!!write(*,*) '#Master::: Calling routine to copy variable floyg'
!!
!!call OmpCopyPointerfloyg
!!
!!write(*,*) '#Master::: Calling routine to copy variable conxe'
!!
!!call OmpCopyPointerconxe
!!
!!write(*,*) '#Master::: Calling routine to copy variable conye'
!!
!!call OmpCopyPointerconye
!!
!!write(*,*) '#Master::: Calling routine to copy variable conxi'
!!
!!call OmpCopyPointerconxi
!!
!!write(*,*) '#Master::: Calling routine to copy variable conyi'
!!
!!call OmpCopyPointerconyi
!!
!!write(*,*) '#Master::: Calling routine to copy variable conxg'
!!
!!call OmpCopyPointerconxg
!!
!!write(*,*) '#Master::: Calling routine to copy variable conyg'
!!
!!call OmpCopyPointerconyg
!!
!!write(*,*) '#Master::: Calling routine to copy variable floxge'
!!
!!call OmpCopyPointerfloxge
!!
!!write(*,*) '#Master::: Calling routine to copy variable floyge'
!!
!!call OmpCopyPointerfloyge
!!
!!write(*,*) '#Master::: Calling routine to copy variable conxge'
!!
!!call OmpCopyPointerconxge
!!
!!write(*,*) '#Master::: Calling routine to copy variable conyge'
!!
!!call OmpCopyPointerconyge
!!
!!write(*,*) '#Master::: Calling routine to copy variable visx'
!!
!!call OmpCopyPointervisx
!!
!!write(*,*) '#Master::: Calling routine to copy variable visy'
!!
!!call OmpCopyPointervisy
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxe'
!!
!!call OmpCopyPointerhcxe
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcye'
!!
!!call OmpCopyPointerhcye
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxij'
!!
!!call OmpCopyPointerhcxij
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcyij'
!!
!!call OmpCopyPointerhcyij
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxg'
!!
!!call OmpCopyPointerhcxg
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcyg'
!!
!!call OmpCopyPointerhcyg
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxi'
!!
!!call OmpCopyPointerhcxi
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxineo'
!!
!!call OmpCopyPointerhcxineo
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcyi'
!!
!!call OmpCopyPointerhcyi
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcxn'
!!
!!call OmpCopyPointerhcxn
!!
!!write(*,*) '#Master::: Calling routine to copy variable hcyn'
!!
!!call OmpCopyPointerhcyn
!!
!!write(*,*) '#Master::: Calling routine to copy variable kxbohm'
!!
!!call OmpCopyPointerkxbohm
!!
!!write(*,*) '#Master::: Calling routine to copy variable kybohm'
!!
!!call OmpCopyPointerkybohm
!!
!!write(*,*) '#Master::: Calling routine to copy variable dif_use'
!!
!!call OmpCopyPointerdif_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable difp_use'
!!
!!call OmpCopyPointerdifp_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable dif2_use'
!!
!!call OmpCopyPointerdif2_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable tray_use'
!!
!!call OmpCopyPointertray_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable trax_use'
!!
!!call OmpCopyPointertrax_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable kye_use'
!!
!!call OmpCopyPointerkye_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable kyi_use'
!!
!!call OmpCopyPointerkyi_use
!!call OmpCopyPointerkxe_use
!!call OmpCopyPointerkxi_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable dutm_use'
!!
!!call OmpCopyPointerdutm_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable vy_use'
!!
!!call OmpCopyPointervy_use
!!
!!write(*,*) '#Master::: Calling routine to copy variable vy_cft'
!!
!!call OmpCopyPointervy_cft
!!
!!write(*,*) '#Master::: Calling routine to copy variable vyte_cft'
!!
!!call OmpCopyPointervyte_cft
!!
!!write(*,*) '#Master::: Calling routine to copy variable vyti_cft'
!!
!!call OmpCopyPointervyti_cft
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuiz'
!!
!!call OmpCopyPointernuiz
!!
!!write(*,*) '#Master::: Calling routine to copy variable nucx'
!!
!!call OmpCopyPointernucx
!!
!!write(*,*) '#Master::: Calling routine to copy variable nucxi'
!!
!!call OmpCopyPointernucxi
!!
!!write(*,*) '#Master::: Calling routine to copy variable nueli'
!!
!!call OmpCopyPointernueli
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuelg'
!!
!!call OmpCopyPointernuelg
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuix'
!!
!!call OmpCopyPointernuix
!!
!!write(*,*) '#Master::: Calling routine to copy variable nurc'
!!
!!call OmpCopyPointernurc
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuvl'
!!
!!call OmpCopyPointernuvl
!!
!!write(*,*) '#Master::: Calling routine to copy variable eqp'
!!
!!call OmpCopyPointereqp
!!
!!write(*,*) '#Master::: Calling routine to copy variable eqpg'
!!
!!call OmpCopyPointereqpg
!!
!!write(*,*) '#Master::: Calling routine to copy variable eeli'
!!
!!call OmpCopyPointereeli
!!
!!write(*,*) '#Master::: Calling routine to copy variable pradhyd'
!!
!!call OmpCopyPointerpradhyd
!!
!!write(*,*) '#Master::: Calling routine to copy variable eta1'
!!
!!call OmpCopyPointereta1
!!
!!write(*,*) '#Master::: Calling routine to copy variable rtaue'
!!
!!call OmpCopyPointerrtaue
!!
!!write(*,*) '#Master::: Calling routine to copy variable dclass_e'
!!
!!call OmpCopyPointerdclass_e
!!
!!write(*,*) '#Master::: Calling routine to copy variable dclass_i'
!!
!!call OmpCopyPointerdclass_i
!!
!!write(*,*) '#Master::: Calling routine to copy variable visxneo'
!!
!!call OmpCopyPointervisxneo
!!
!!write(*,*) '#Master::: Calling routine to copy variable visvol_v'
!!
!!call OmpCopyPointervisvol_v
!!
!!write(*,*) '#Master::: Calling routine to copy variable visvol_q'
!!
!!call OmpCopyPointervisvol_q
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuii'
!!
!!call OmpCopyPointernuii
!!
!!write(*,*) '#Master::: Calling routine to copy variable nuiistar'
!!
!!call OmpCopyPointernuiistar
!!
!!write(*,*) '#Master::: Calling routine to copy variable alfneo'
!!
!!call OmpCopyPointeralfneo
!!
!!write(*,*) '#Master::: Calling routine to copy variable k2neo'
!!
!!call OmpCopyPointerk2neo
!!
!!write(*,*) '#Master::: Calling routine to copy variable ktneo'
!!
!!call OmpCopyPointerktneo
!!
!!write(*,*) '#Master::: Calling routine to copy variable snic'
!!
!!call OmpCopyPointersnic
!!
!!write(*,*) '#Master::: Calling routine to copy variable sniv'
!!
!!call OmpCopyPointersniv
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorc'
!!
!!call OmpCopyPointerpsorc
!!call OmpCopyPointerpsor
!!call OmpCopyPointerpsorxrc
!!call OmpCopyPointerpsorxr
!!
!!write(*,*) '#Master::: Calling routine to copy variable psor_tmpov'
!!
!!call OmpCopyPointerpsor_tmpov
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorgc'
!!
!!call OmpCopyPointerpsorgc
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorg'
!!
!!call OmpCopyPointerpsorg
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorrgc'
!!
!!call OmpCopyPointerpsorrgc
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorrg'
!!
!!call OmpCopyPointerpsorrg
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorcxgc'
!!
!!call OmpCopyPointerpsorcxgc
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorcxg'
!!
!!call OmpCopyPointerpsorcxg
!!
!!write(*,*) '#Master::: Calling routine to copy variable psori'
!!
!!call OmpCopyPointerpsori
!!
!!write(*,*) '#Master::: Calling routine to copy variable psordis'
!!
!!call OmpCopyPointerpsordis
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorbgg'
!!
!!call OmpCopyPointerpsorbgg
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorbgz'
!!
!!call OmpCopyPointerpsorbgz
!!
!!write(*,*) '#Master::: Calling routine to copy variable erliz'
!!
!!call OmpCopyPointererliz
!!
!!write(*,*) '#Master::: Calling routine to copy variable erlrc'
!!
!!call OmpCopyPointererlrc
!!
!!write(*,*) '#Master::: Calling routine to copy variable vsoreec'
!!
!!call OmpCopyPointervsoreec
!!
!!write(*,*) '#Master::: Calling routine to copy variable vsoree'
!!
!!call OmpCopyPointervsoree
!!
!!write(*,*) '#Master::: Calling routine to copy variable pwrebkg'
!!
!!call OmpCopyPointerpwrebkg
!!
!!write(*,*) '#Master::: Calling routine to copy variable pwribkg'
!!
!!call OmpCopyPointerpwribkg
!!
!!write(*,*) '#Master::: Calling routine to copy variable wjdote'
!!
!!call OmpCopyPointerwjdote
!!
!!write(*,*) '#Master::: Calling routine to copy variable smoc'
!!
!!call OmpCopyPointersmoc
!!
!!write(*,*) '#Master::: Calling routine to copy variable smov'
!!
!!call OmpCopyPointersmov
!!
!!write(*,*) '#Master::: Calling routine to copy variable msor'
!!
!!call OmpCopyPointermsor
!!
!!write(*,*) '#Master::: Calling routine to copy variable msorxr'
!!
!!call OmpCopyPointermsorxr
!!
!!write(*,*) '#Master::: Calling routine to copy variable seec'
!!
!!call OmpCopyPointerseec
!!
!!write(*,*) '#Master::: Calling routine to copy variable seev'
!!
!!call OmpCopyPointerseev
!!
!!write(*,*) '#Master::: Calling routine to copy variable seic'
!!
!!call OmpCopyPointerseic
!!
!!write(*,*) '#Master::: Calling routine to copy variable seiv'
!!
!!call OmpCopyPointerseiv
!!
!!
!!
!!call OmpCopyPointerresco
!!
!! to copy variable resng'
!!
!!call OmpCopyPointerresng
!!
!!write(*,*) '#Master::: Calling routine to copy variable reseg'
!!
!!call OmpCopyPointerreseg
!!
!!write(*,*) '#Master::: Calling routine to copy variable resmo'
!!
!!call OmpCopyPointerresmo
!!
!!write(*,*) '#Master::: Calling routine to copy variable resee'
!!
!!call OmpCopyPointerresee
!!
!!write(*,*) '#Master::: Calling routine to copy variable resei'
!!
!!call OmpCopyPointerresei
!!
!!write(*,*) '#Master::: Calling routine to copy variable resphi'
!!
!!call OmpCopyPointerresphi
!!
!!write(*,*) '#Master::: Calling routine to copy variable sng_ue'
!!
!!call OmpCopyPointersng_ue
!!
!!write(*,*) '#Master::: Calling routine to copy variable seg_ue'
!!
!!call OmpCopyPointerseg_ue
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorold'
!!
!!call OmpCopyPointerpsorold
!!
!!write(*,*) '#Master::: Calling routine to copy variable psorxrold'
!!
!!call OmpCopyPointerpsorxrold
!!
!!write(*,*) '#Master::: Calling routine to copy variable msorold'
!!
!!call OmpCopyPointermsorold
!!
!!write(*,*) '#Master::: Calling routine to copy variable msorxrold'
!!
!!call OmpCopyPointermsorxrold
!!
!!write(*,*) '#Master::: Calling routine to copy variable yldot_pert'
!!
!!call OmpCopyPointeryldot_pert
!!call OmpCopyPointernzloc
!!call OmpCopyPointerimpradloc
!!call OmpCopyPointerpwrzec
!!call OmpCopyPointerpwrze
!!call OmpCopyPointerpradc
!!call OmpCopyPointerpradcff
!!call OmpCopyPointerprad
!!call OmpCopyPointerpradzc
!!call OmpCopyPointerpradz
!!call OmpCopyPointerna
!!call OmpCopyPointerntau
!!call OmpCopyPointernratio
!!call OmpCopyPointertempa
!!call OmpCopyPointerden
!!call OmpCopyPointergradp
!!call OmpCopyPointergradt
!!call OmpCopyPointerdztot
end subroutine OmpCopyPointerbbb2


subroutine ompCopyConvert
      use OmpCopybbb
      implicit none
      call OmpCopyPointerne
      call OmpCopyPointernit
      call OmpCopyPointernm
      call OmpCopyPointernz2
      call OmpCopyPointerni
      call OmpCopyPointerng
      call OmpCopyPointerlng
      call OmpCopyPointerte
      call OmpCopyPointertg
      call OmpCopyPointerti
      call OmpCopyPointerphi
      call OmpCopyPointerup
      call OmpCopyPointerpr
      call OmpCopyPointerzeff
      call OmpCopyPointerpri
      call OmpCopyPointerpre
      call OmpCopyPointerznot
      call OmpCopyPointertg
      call OmpCopyPointerpg
      call OmpCopyPointergprx
      call OmpCopyPointerney0
      call OmpCopyPointerney1
      call OmpCopyPointernity0
      call OmpCopyPointernity1
      call OmpCopyPointergpry
      call OmpCopyPointergpix
      call OmpCopyPointerniy0
      call OmpCopyPointerniy1
      call OmpCopyPointerniy0s
      call OmpCopyPointerniy1s
      call OmpCopyPointerpriy0
      call OmpCopyPointerpriy1
      call OmpCopyPointergpiy
      call OmpCopyPointertey0
      call OmpCopyPointertey1
      call OmpCopyPointertiy0
      call OmpCopyPointertiy1
      call OmpCopyPointerphiy0
      call OmpCopyPointerphiy1
      call OmpCopyPointertiy0s
      call OmpCopyPointertiy1s
      call OmpCopyPointerphiy0s
      call OmpCopyPointerphiy1s
      call OmpCopyPointerngy0
      call OmpCopyPointerngy1
      call OmpCopyPointertgy0
      call OmpCopyPointertgy1
      call OmpCopyPointerpgy0
      call OmpCopyPointerpgy1
      call OmpCopyPointergpex
      call OmpCopyPointergtex
      call OmpCopyPointergtix
      call OmpCopyPointerex
      call OmpCopyPointergpey
      call OmpCopyPointergtey
      call OmpCopyPointergtiy
      call OmpCopyPointerey
      call OmpCopyPointerphiv
      call OmpCopyPointertiv
      call OmpCopyPointertev
      call OmpCopyPointerprev
      call OmpCopyPointerprtv
      call OmpCopyPointerpriv
      end subroutine ompCopyConvert
#endif
subroutine jac_write(filename,neq, jac, jaccol, jacrow)

!  This function serves to output the jacobian for viewing purposes
      character(*)::filename
      integer neq,j,k
      real jac(*)
      integer jacrow(*), jaccol(*)

      open(UNIT=88,FILE=trim(filename),STATUS='REPLACE')

      do j=1,neq
        do k=jacrow(j),jacrow(j+1)-1
          write(88,*)j,jacrow(j),jaccol(k),jac(k)
        end do
      end do
      close(88)
      end subroutine jac_write

subroutine EvalDumpJac(FileName,neq,yl,yldot00)

    Use ParallelSettings,only:OMPParallelJac,MPIParallelJac
    Use OMPJacSettings,only:iidebugprint,ivdebugprint,DebugJac,ForceSerialCheck,CheckJac,DumpFullJac, DumpJac
    Use Cdv,only:exmain_aborted
    implicit none
    ! ... Input arguments:
    character(*)::filename
    integer,intent(in):: neq      !      total number of equations (all grid points)
    real,intent(inout)   :: yl(*)          ! dependent variables
    real,intent(in)   :: yldot00(neq+2) ! right-hand sides evaluated at yl

    real   :: jaccopy(neq*neq)     ! nonzero Jacobian elements
    integer:: jacopy(neq*neq)   ! col indices of nonzero Jacobian elements
    integer:: iacopy(neq+1)   ! pointers to beginning of each row in jac,ja
    real wk(neq)     ! work space available to this subroutine
    write(*,*) '--- Evaluating full jacobian (serial) for analysis of bandwidth'
    call jac_calc (neq, 0, yl, yldot00, neq, neq, wk,neq*neq, jaccopy, jacopy, iacopy)
    write(*,*) '--- Dumping jacobian for analysis of bandwidth'
    call jac_write(FileName,neq, jaccopy, jacopy, iacopy)
end subroutine EvalDumpJac
