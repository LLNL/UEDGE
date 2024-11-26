        SUBROUTINE right_boundary(neq, yl, yldot)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: neq
        REAL, INTENT(OUT) :: yl(neq), yldot(neq)        

          Use(Dim)      # nx,ny,nhsp,nzspt,nzsp,nisp,ngsp,nusp,nxpt
          Use(Share)    # nxpt,nxc,geometry,cutlo,islimon,ix_lim,iy_lims
                        # isudsym
          Use(Xpoint_indices) # ixlb,ixpt1,ixpt2,ixrb,iysptrx1,iysptrx2
          Use(Math_problem_size)   # neqmx 
          Use(Phyvar)
          Use(UEpar)    # isnewpot,r0slab,cslim,dcslim,csfaclb,csfacrb,csfacti,
                        # isnion,isupon,isteon,istion,isngon,isnionxy,isuponxy,
                        # isteonxy,istionxy,isngonxy,isphionxy, ismolcrm
          Use(Aux)      # ixmp
          Use(Coefeq)   # fac2sp,cf2ef,exjbdry
          Use(Bcond)    # iflux,ncore,tcoree,tcorei,tbmin,nbmin,ngbmin,
                        # tepltl,tipltl,tepltr,tipltr,
                        # istewc,istiwc,istepfc,istipfc,
                        # tewalli,tiwalli,tewallo,tiwallo,isextrnp,
                        # isextrnpf,isextrtpf,isextrngc,isextrnw,isextrtw,
                        # iflcore,pcoree,pcorei,ifluxni,ckinfl,isupss,
                        # isnwconiix,isnwconoix,nwalli,nwallo,iscpli,iscplo,
                        # fngysi,fngyso,albedoo,albedoi,matwallo,matwalli,
                        # sinphi,isfixlb,nib,teb,tib,nibprof,tebprof,tibprof,
                        # engbsr,epsbs,rlimiter,ngcore,isngcore,isutcore,
                        # ixfixnc,incixc,isupcore,isfixrb,chemsputi,chemsputo
                        # islbcn,islbcu,islbce,islbci,islbcg,isexunif
                        # fchemygwi,fchemylb,fphysylb,fchemygwo,fchemyrb,fphysyrb
                        # xcnearlb,xcnearrb,openbox,fqpsatlb,fqpsatrb
                        # cfueb,ikapmod,cfvytanbc
          Use(Parallv)  # nxg,nyg
          Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,j1,j2,j3,j4,j5,j6,j7,xlinc
          Use(Comgeo)   # gx,gy,gyf,sx,sy,xcwi,xcwo,yylb,rrv,sygytotc,isixcore
          Use(Compla)   # mi, mg
          Use(Comflo)   # fqx,fqy,fnix,fniy,feex,feey,feix,feiy,fngx,fngy
                        # fdiaxlb, fdiaxrb
          Use(Conduc)   # visx
          Use(Indexes)
          Use(Ynorm)    # temp0,nnorm,ennorm
          Use(Poten)    # newbcl,newbcr,bcee,bcei,rsigpl,bcel,bcer,bcil,bcir,
                        # kappal,kappar,bctype,phi0l,phi0r,isfdiax
          Use(Rccoef)   # recylb,recyrb,alblb,albrb,recycw,sputtr,
                        # recycm,recyce,recycmlb,recycmrb,recyllim,recyrlim
          Use(Bfield)   # rbfbt,btot
          Use(Imprad)   # isimpon
          Use(Impurity_source_flux)   # fnzysi,fnzyso
          Use(Gradients)  # ey
          Use(RZ_grid_info)      # rm
          Use(Indices_domain_dcl)   #ixmxbcl,ixmnbcl,iymxbcl,iymnbcl,ispwrbcl
          Use(Interp)
          Use(Jacaux)   # yldot_diag
          Use(Npes_mpi) # npes
          Use(Indices_domain_dcg) # ndomain,ispwrbc
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

          Use(MCN_dim)
          Use(MCN_sources) # edisspl, edisspr, cmntgpl, cmntgpl
          Use(Utilities)


c...  local scalars
          real totfeix, totfeex, kfeix, cosphi,
     .     ueb, nbound, tbound, ut0, sumb, feeytotc, feiytotc,
     .     r_major, fniytotc, fng_chem, vbound, eng_sput, flx_incid,
     .     yld_chm, t0p, zflux_chm, fqytotc, flux_inc,
     .     totfnex, totfnix, fqpsate, qpfac, aq, arglgphi, faceel,
     .     faceel2, csfac, lambdae, uztotc, uztotc1, uztotc2,
     .     fngytotc, fmiytotc, sytotc, f_cgpld, sfeeytotc, sfeiytotc,
     .     vxa, ta0, flxa
          integer ii,isphion2, nzsp_rt, jz
          real hflux, zflux
          integer ifld, ihyd, iimp, ix_fl_bc, ixc1, igsp2
          real dif_imp_flux, fng_alb, fngyw, nharmave
          real upbdry, upbdry1, upbdry2, uugoal, fniy_recy, lengg, xtotc
          integer ixt, ixt1, ixt2, ixt3, jx, ixc, ierr
          integer ixtl, ixtl1, ixtr,ixtr1
          #Former Aux module variables
          integer ix,iy,igsp,iv,iv1,iv2,iv3,iv4,ix1,ix2,ix3,ix4
          real osmw
          real t0

          real yld96, kappa










c********************************************************************
c...  First, check if isfixrb=2 for using symmetry BC at ix = nx+1
c*******************************************************************
      if (i6 .ge. nx+1 .and. isfixrb(1) .gt. 0) then  #begin symmetry BC at nx+1
       do iy = j2, j5
c...  First do the ion density
         do ifld = 1, nisp
            if(isnionxy(nx+1,iy,ifld) .eq. 1) then
               iv1 = idxn(nx+1,iy,ifld)
               yldot(iv1) = nurlxn *
     .              (nib(ifld)*nibprof(iy)-ni(nx+1,iy,ifld))/n0(ifld)
               if(isfixrb(1).eq.2) yldot(iv1) = nurlxn * (1/n0(ifld)) *
     .              (ni(nx,iy,ifld) - ni(nx+1,iy,ifld))
            endif
         enddo

c...  Now do the parallel velocity
         do ifld = 1, nusp
            if(isuponxy(nx+1,iy,ifld) .eq. 1) then
               iv2 = idxu(nx+1,iy,ifld)
               yldot(iv2) = nurlxu*(up(nx,iy,ifld) - up(nx+1,iy,ifld))/
     .                                                         vpnorm
               if(isfixrb(1).eq.2) yldot(iv2) = nurlxu*
     .                            (0. - up(nx+1,iy,ifld))/vpnorm
               if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
                  cs = sqrt( (te(nx+1,iy)+ti(nx+1,iy))/mi(ifld) )
                  yldot(iv2) = nurlxu*(cs -up(nx+1,iy,ifld))/vpnorm
               endif
            endif
         enddo

c...  now do the gas and temperatures
         if(isteonxy(nx+1,iy) .eq. 1) then
           iv1 = idxte(nx+1,iy)
           yldot(iv1) = nurlxe * ne(nx+1,iy) *
     .                     (teb*ev*tebprof(iy) - te(nx+1,iy))/ennorm
           if(isfixrb(1).eq.2) yldot(iv1) = nurlxe * ne(nx+1,iy) *
     .                               (te(nx,iy) - te(nx+1,iy))/ennorm
           if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
              yldot(iv1) = - nurlxe*(feex(nx,iy)/sx(nx,iy) - bcee*
     .                           ne(nx+1,iy)*vex(nx,iy)*te(nx+1,iy))/
     .                                 (vpnorm*ennorm)
           endif
         endif
         if(istionxy(nx+1,iy) .eq. 1) then
            iv2 = idxti(nx+1,iy)
            yldot(iv2) = nurlxi * ne(nx+1,iy) *
     .                     (tib*ev*tibprof(iy) - ti(nx+1,iy))/ennorm
            if(isfixrb(1).eq.2) yldot(iv2) = nurlxi * ne(nx+1,iy) *
     .                               (ti(nx,iy) - ti(nx+1,iy))/ennorm
            if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
               yldot(iv2) = -nurlxi*
     .          ( feix(nx,iy) - bcei*ti(nx+1,iy)*fac2sp*fnix(nx,iy,1) ) /
     .                                         (vpnorm*ennorm*sx(nx,iy))
            endif
         endif
         do igsp = 1, nhgsp # not valid for ngsp > nhgsp; only on hydrog. gas
            if (isngonxy(nx+1,iy,igsp) .eq. 1) then
               iv = idxg(nx+1,iy,igsp)
               yldot(iv) = nurlxg * (ngbackg(igsp) - 
     .                                     ng(nx+1,iy,igsp)) / n0g(igsp)
               if(isfixrb(1).eq.2) yldot(iv) = nurlxg * 
     .                     (ng(nx,iy,igsp) - ng(nx+1,iy,igsp))/n0g(igsp)
               if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
                  flux_inc = fac2sp*fnix(nx,iy,1)
                  if (ishymol.eq.1 .and. igsp.eq.2) then
                    flxa= ismolcrm*(1-albrb(iy,1,nxpt))* onesided_maxwellian(
     .                  engbsr*tg(nx,iy,1), ng(nx,iy,1), mg(1), 
     .                  sx(nx, iy), engbsr*tgmin*ev
     .              )
                    if (isupgon(1) .eq. 1) then  # two atoms for one molecule
                      flux_inc = 0.5*( fnix(nx,iy,1) + fnix(nx,iy,2) -flxa) 
                    else
                      flux_inc = 0.5*( fnix(nx,iy,1) + fngx(nx,iy,1) -flxa) 
                    endif
                  endif
                  osmw = onesided_maxwellian(
     .                  engbsr*tg(nx,iy,1), ng(nx,iy,igsp), mg(igsp),
     .                  isoldalbarea*sx(nx,iy) + (1-isoldalbarea)*sxnp(nx,iy),
     .                  engbsr*tgmin*ev
     .            )      
                  yldot(iv) = -nurlxg * ( fngx(nx,iy,igsp)
     .                      + fngxrb_use(iy,igsp,1)
     .                      - fngxsrb(iy,igsp,1) 
     .                      + recyrb(iy,igsp,1)*flux_inc
     .                      - (1-albrb(iy,igsp,nxpt))*osmw
     .                  ) / (vpnorm*n0g(igsp)*sx(nx,iy))
               endif
            endif
         enddo      # end of igsp loop over gas


c ... Neutral temperature - test if tg eqn is on, then set BC
	 do igsp = 1, ngsp
           if (istgonxy(nx+1,iy,igsp) == 1) then
             iv = idxtg(nx+1,iy,igsp)
             yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(nx+1,iy,igsp))/(temp0*ev)
             if(isfixrb(1)==2) then #just above applies if isfixrb=1
               yldot(iv)=nurlxg*(tg(nx,iy,igsp)-tg(nx+1,iy,igsp))/
     .                                                      (temp0*ev)
             endif
             if(isfixrb(1)==2 .and. yyrb(iy,1) > rlimiter) then
               yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(nx+1,iy,igsp))/(temp0*ev)
             endif
           endif
         enddo

         if (isphionxy(nx+1,iy) .eq. 1) then
            iv = idxphi(nx+1,iy)
            yldot(iv) = nurlxp*(phi(nx,iy) - phi(nx+1,iy))/temp0
         endif

       enddo         # end of main loop starting with iy = j2, j5

      endif         # end of ix = nx+1, isfixrb.ne.0 boundary conditions

c...  If isfixrb=2, check if i2,i5 range for yldot in pandf includes ixpt1(1);
c...  the overwrite ix=ixpt1(1) up eqn with up(ixpt1(1)) --> 0 eqn.
      if (isfixrb(1) .eq. 2) then
         if (i2.le.ixpt1(1) .and. i5.ge.ixpt1(1) .and. j2.le.iysptrx1(1)) then  
           do ifld = 1, nusp
             if (isupcore(ifld).eq.0) then
               do iy = 0+1-iymnbcl, iysptrx1(1)
                 if(isuponxy(ixpt1(1),iy,ifld)==1) then
                   iv = idxu(ixpt1(1),iy,ifld)
                   yldot(iv) = nurlxu*(0.-up(ixpt1(1),iy,ifld))/vpnorm
                 endif
               enddo
             endif
           enddo
         endif 
      endif                             #end symm. BC at ix=nx+1 for isfixrb=2

c************************************************************************
c     begin standard divertor plate conditions for right boundaries
c************************************************************************
      do jx = 1, nxpt # loop over nxpt mesh regions
      if (xcnearrb .or. openbox) then
c     First, the density equations --
        do ifld = 1, nisp
          if ( i6.ge.(ixrb(jx)+1-isextrnp) .and. isfixrb(jx)==0) then
          do iy = j2, j5

            ixt  = ixrb(jx) + 1    # analog of ix=nx+1
            ixt1 = ixm1(ixt,iy)    # analog of ix=nx
            ixt2 = ixm1(ixt1,iy)   # analog of ix=nx-1
            ixt3 = ixm1(ixt2,iy)   # analog of ix=nx-2

        if (isnionxy(ixt,iy,ifld)==1) then
              iv1 = idxn(ixt,iy,ifld)
              if (isupgon(1)==1 .and. zi(ifld)==0.0) then   ## neutrals
                if (recyrb(iy,1,jx) .gt. 0.) then           # recycling
                  osmw = onesided_maxwellian(
     .                  tg(ixt1,iy,1), ni(ixt1,iy,ifld), mi(ifld), 
     .                  isoldalbarea*sx(ixt1,iy) + (1-isoldalbarea)*sxnp(ixt1,iy),
     .                  tgmin*ev
     .            )
                  yldot(iv1) = nurlxg * (fnix(ixt1,iy,ifld) 
     .                  + recyrb(iy,1,jx)*fnix(ixt1,iy,1)
     .                  + fngxrb_use(iy,1,jx)
     .                  - (1-albrb(iy,1,jx))*osmw
     .                  - fngxsrb(iy,1,jx) 
     .              ) / (vpnorm*n0(ifld)*sx(ixt1,iy))
                elseif (recyrb(iy,1,jx) <=  0. .and. 
     .                  recyrb(iy,1,jx) >= -1.) then   # recyrb is albedo
                  osmw = onesided_maxwellian(
     .                  tg(ixt1,iy,1), ni(ixt,iy,ifld), mi(1),
     .                  sx(ixt1,iy), tgmin*ev
     .            )
                  yldot(iv1) = nurlxg * ( fnix(ixt1,iy,ifld) 
     .                      - (1+recyrb(iy,1,jx))*osmw
     .                 ) / (vpnorm*n0(ifld)*sx(ixt1,iy))
                elseif (recyrb(iy,1,jx) < -1. .and.
     .                  recyrb(iy,1,jx) > -2.) then    #fix density ngrfix
                  yldot(iv1)=nurlxg*(ngrfix - ni(ixt,iy,ifld))/n0(ifld)
                elseif (recyrb(iy,1,jx) <= -2.) then   #zero gradient
                  yldot(iv1) = nurlxn*(ni(ixt1,iy,ifld)-ni(ixt,iy,ifld))/
     .                                                          n0(ifld)
                endif # end if-test on recyrb
              else                                     ## ions
                if (isextrnp==0) then                  # zero x-gradient
                  yldot(iv1) = nurlxn*(ni(ixt1,iy,ifld)-ni(ixt,iy,ifld))/
     .                                                          n0(ifld)
                else                                   # extrapolation
                  nbound =  ni(ixt1,iy,ifld) + gxf(ixt2,iy)*
     .                   (ni(ixt1,iy,ifld)-ni(ixt2,iy,ifld))/gxf(ixt1,iy)
                  nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                    ni(ixt,iy,ifld)-1) ) ) + 0.2*ni(ixt,iy,ifld)
                  yldot(iv1) = nurlxn*(nbound-ni(ixt,iy,ifld))/n0(ifld)
                endif # end if-test on isextrnp
              endif   # end if-test on isupgon and zi
            endif     # end if-test on isnionxy  
          enddo       # end do-loop in iy
          endif       # end if-test on i6 and isfixrb and isnion
        enddo         # end do-loop on ifld


c     Now do the parallel velocity and other variables --
      if ( (i6 .ge. (ixrb(jx)+1)) .and. isfixrb(jx)==0 ) then
      do iy = j2, j5 # begin big do-loop on iy

        ixt  = ixrb(jx) + 1    # analog of ix=nx+1
        ixt1 = ixm1(ixt,iy)    # analog of ix=nx
        ixt2 = ixm1(ixt1,iy)   # analog of ix=nx-1
        ixt3 = ixm1(ixt2,iy)   # analog of ix=nx-2

        kfeix = 0.
        sumb = 0.
        do ifld = 1, nfsp  # set up generalized Bohm condition
          upi(ixt,iy,ifld) = upi(ixt1,iy,ifld) # not set before; upi not var
          upi(ixt,iy,1   ) =  up(ixt ,iy,1   ) # need to keep up(,,1) as var
          ueb = cfueb*( cf2ef*v2ce(ixt1,iy,ifld)*rbfbt(ixt,iy) -
     .            vytan(ixt1,iy,ifld) ) / rrv(ixt1,iy) 
          sumb = sumb + ni(ixt,iy,ifld)*zi(ifld)**2*te(ixt,iy) /
     .                (mi(ifld)*(upi(ixt,iy,ifld)+ueb)**2 - ti(ixt,iy))
        enddo # end do-loop on ifld for Bohm condition
        sumb = sqrt(abs(sumb/ne(ixt,iy)))

c       Next, the momentum equations --
        do ifld = 1, nusp
          if (isuponxy(ixt,iy,ifld)==1) then
            iv2 = idxu(ixt1,iy,ifld)  #ixt1 ~ nx
	    iv = idxu(ixt,iy,ifld)    #ixt ~ nx+1
            cs = csfacrb(ifld,jx)*sqrt( (te(ixt,iy)+
     .                                  csfacti*ti(ixt,iy))/mi(ifld) )
            if (isupgon(1)==1 .and. zi(ifld)==0.0) then  ## neutrals
              if (recycmrb(iy,1,jx) > -9.9) then  # backscatter with recycm
                yldot(iv2) = -nurlxu*(recycmrb(iy,1,jx)*up(ixt1,iy,1) + 
     .                                       up(ixt1,iy,ifld))/vpnorm
              elseif (recycmrb(iy,1,jx) <= -9.9 .and. 
     .                      recycmrb(iy,1,jx) > -10.1) then # zero x-gradient
                yldot(iv2) = nurlxu*(up(ixt2,iy,ifld) -
     .                                          up(ixt1,iy,ifld))/vpnorm
              else  #neutral thermal flux to wall if recycm < -10.1
c...              if up > 0, leave unchanged; if up<0, big reduction
                osmw = onesided_maxwellian(
     .                  tg(ixt,iy,1), 0.5*(nm(ixt1,iy,ifld)+nm(ixt,iy,ifld)),
     .                  mi(ifld), sx(ixt1,iy), tgmin*ev
     .          )
                yldot(iv2) = -nurlxu*( fmix(ixt1,iy,ifld) 
     .                  - up(ixt,iy,ifld)*cgmompl*osmw
     .              ) / (vpnorm*fnorm(ifld)*sx(ixt1,iy))
              endif
              yldot(iv) =nurlxu*(up(ixt1,iy,ifld)-up(ixt,iy,ifld))/vpnorm
            else                                         ## ions
              ueb = cfueb*( cf2ef*v2ce(ixt1,iy,ifld)*rbfbt(ixt,iy) -
     .                vytan(ixt1,iy,ifld) ) / rrv(ixt1,iy) 
              yldot(iv2) = -nurlxu*(sumb - 1.)  # multispecies Bohm
              if (isbohmms==0) then          # simple Bohm condition
                yldot(iv2) = nurlxu * (cs-ueb-up(ixt1,iy,ifld))/vpnorm
              endif
              if(isupss(ifld)==1 .and. up(ixt2,iy,ifld)+ueb .gt. cs)
                                             # dup/dx=0 if supersonic
     .          yldot(iv2) = nurlxu*(up(ixt2,iy,ifld)-up(ixt1,iy,ifld))/
     .                                                           vpnorm
              if (isupss(ifld)==-1) then     # slip boundary conditions
                yldot(iv2) = nurlxu*(up(ixt2,iy,ifld)-up(ixt1,iy,ifld))/
     .                                                           vpnorm
              elseif (isupss(ifld)==-2) then # extrap. + no pos. uu
                vbound = up(ixt2,iy,ifld) - gx(ixt2,iy)*
     .                   (up(ixt3,iy,ifld)-up(ixt2,iy,ifld))/gx(ixt1,iy) 
                vbound = max(vbound, -ueb)   # forces uu & fnix >= 0
                yldot(iv2) = nurlxu*(vbound - up(ixt1,iy,ifld))/vpnorm
              elseif (isupss(ifld)==-3) then # modified Bohm condition
                vbound = -ueb +2*cs*uu(ixt2,iy,ifld)/
     .                             (uu(ixt2,iy,ifld)+rrv(ixt,iy)*cs)
                vbound = max(vbound, -ueb)   # forces uu & fnix >= 0
                yldot(iv2) = nurlxu * (vbound-up(ixt1,iy,ifld))/vpnorm
              endif # end if-test on isupss

c Finally set unused up(ixt,,) = up(ixt1,,); note ixt~nx+1
              yldot(iv) =nurlxu*(up(ixt1,iy,ifld)-up(ixt,iy,ifld))/vpnorm
            endif # end if-test on isupgon
          endif # end if-test on isupon
          if (zi(ifld)==0.0) then 
            kfeix = kfeix - cftiexclg*cfvcsx(ifld)*0.5*sx(ixt1,iy)
     .                     *visx(ixt1,iy,ifld)*gx(ixt1,iy)
     .           *( up(ixt1,iy,ifld)**2 - up(ixt2,iy,ifld)**2 ) 
          else
            kfeix = kfeix - cfvcsx(ifld)*0.5*sx(ixt1,iy)
     .                     *visx(ixt1,iy,ifld)*gx(ixt1,iy)
     .           *( up(ixt1,iy,ifld)**2 - up(ixt2,iy,ifld)**2 )
          endif
        enddo # end do-loop on ifld

c       Next, the potential equation --
        if(isphion+isphiofft .eq. 1) then
cc           fqpsatrb(iy,jx) = qe*isfdiax*( ne(ixt,iy)*v2ce(ixt1,iy,1)*
cc     .                           rbfbt(ixt,iy)*sx(ixt1,iy) + fdiaxrb(iy,jx) )
cc           do ifld = 1, nfsp   # note fqp,fqpsat are poloidal proj. of || curr
cc            fqpsatrb(iy,jx) = fqpsatrb(iy,jx) + qe*zi(ifld)*ni(ixt,iy,ifld)*
cc     .                        upi(ixt1,iy,ifld)*sx(ixt1,iy)*rrv(ixt1,iy)  
cc           enddo
           lambdae = 2e16*(te(ixt,iy)/ev)**2/ne(ixt,iy)  #approx mfp
           kincorrb(iy,jx) = 1./(1 + cfkincor*(lambdae/lcone(ixt,iy))*
     .                                       abs(ev*phi(ixt,iy)/te(ixt,iy)))
           fqpsate = qe*ne(ixt,iy)*sqrt(te(ixt,iy)/(2*pi*me))*
     .                               kincorrb(iy,jx)*sx(ixt1,iy)*rrv(ixt1,iy)
c          NOTE: by definition, fqpsate is always > 0
           if (ikapmod==0) then
cc           if (fqp(ixt1,iy) > 0.) then #limit to saturation current
cc             fqp(ixt1,iy)= ( abs((fqp(ixt1,iy)*fqpsatrb(iy,jx)))**exjbdry/
cc     .                  (abs(fqp(ixt1,iy))**exjbdry +
cc     .                   abs(fqpsatrb(iy,jx))**exjbdry) )**(1/exjbdry)
cc           endif
              if ( fqpsatrb(iy,jx)-(1.-gamsec)*fqp(ixt1,iy) > 0) then  # force +ve log argument
                arglgphi=(((fqpsatrb(iy,jx)-(1.-gamsec)*fqp(ixt1,iy))
     .                             /fqpsate)**2 + exp(-kappamx)**2)**(0.5)
              else
                arglgphi = exp(-kappamx)
              endif
              if (iskaprex.eq.0) kappar(iy,jx) = - log(arglgphi)
cccTDR            if ((isnewpot==1) .and. ((iy==1) .or. (iy==ny))) then
           elseif (ikapmod==1) then
              kappar(iy,jx) = kappa( fqpsatrb(iy,jx), fqpsate, 
     .                                       (1.-gamsec)*fqp(ixt1,iy) )
           endif
           if (isphionxy(ixt,iy) .eq. 1) then
              iv = idxphi(ixt,iy)
              if ((isnewpot==100) .and. ((iy==1) .or. (iy==ny))) then
                 continue
              elseif (isphirbc==1)
                 yldot(iv) = - nurlxp * (phi(ixt,iy) - phi0r(iy,jx))/temp0 
              else
                 yldot(iv) = -nurlxp*(1.-bctype(iy))*
     .              ( phi(ixt,iy)-kappar(iy,jx)*te(ixt,iy)/ev
     .                                     -phi0r(iy,jx) )/temp0
     .                       -nurlxp*bctype(iy)*(1.-gamsec)*fqp(ixt1,iy)
     .                                 / (fqpsatrb(iy,jx)+cutlo)
              endif
           endif
	else    # corresponding to (isphion+isphiofft .ne. 1)
           fqpsate = 0.
	   kappar(iy,jx) = 3.  # estimate for e*phi/Te; may be used for sputtering
        endif
        isphion2 = isphion + isphiofft
        bcer(iy,jx) = (1-newbcr(jx)*isphion2) * bcee
     .              + newbcr(jx)*isphion2 * (2. + kappar(iy,jx))
        if (iskaprex.eq.1) bcer(iy,jx) = 2. + kappar(iy,jx)
        bcir(iy,jx) = (1-newbcr(jx)*isphion2) * bcei
     .              + newbcr(jx)*isphion2 * (2.5)

        t0 = te(ixt,iy)/ev
        f_cgpld = .5*(1. - cos(pi*(t0 - temin)/(.3 - temin)))
                                # ramp cgpld effect down at low guard cell
                                # Te; rampdown might not be necesasry if
                                # guard cell were allowed to access negative
                                # temperatures (ETM 7 Jan 2014)
        if(t0 < temin) f_cgpld = 0.
        if(t0 > 0.3) f_cgpld = 1.

c   Do the electron temp Eqn -----------------------------------
        if (isteonxy(ixt,iy) == 1) then
	  iv1 = idxte(ixt,iy)
          if (ibctepr == 1) then
             faceel =  bcer(iy,jx)*(fqpsate/qe)*exp(-kappar(iy,jx))
             faceel2 = bcer(iy,jx)*(fqpsate/qe)*exp(-kappamx+2) 
             totfeexr(iy,jx) = feex(ixt1,iy) + cfeexdbo*( 
     .             2.5*fdiaxrb(iy,jx) + floxebgt(ixt1,iy) )*te(ixt,iy)
	     totfnex = ne(ixt,iy)*vex(ixt1,iy)*sx(ixt1 ,iy)
             if (isphion+isphiofft==1) then
              yldot(iv1) = nurlxe*(totfeexr(iy,jx) 
     .                       -faceel*te(ixt,iy)
     .                       -faceel2*(te(ixt,iy)-te(ixt1,iy))
     .                       -cmneut*fnix(ixt1,iy,1)*recycp(1)*eedisspr*ev
     .                                  )/(sx(ixt1,iy)*vpnorm*ennorm)
              else
                osmw = onesided_maxwellian(
     .                  tg(ixt1,iy,1), ng(ixt1,iy,1), mg(1),
     .                  sx(ixt1,iy), tgmin*ev
     .          )
                yldot(iv1) = nurlxe*(totfeexr(iy,jx)
     .                  - totfnex*te(ixt,iy)*bcer(iy,jx)
     .                  - cgpld*f_cgpld*osmw*0.5*ediss*ev
     .                  - cmneut*fnix(ixt1,iy,1)*recycp(1)*eedisspr*ev
     .              ) / (sx(ixt1,iy)*vpnorm*ennorm)
             endif
          elseif (ibctepr .eq. 0) then
             yldot(iv1) = nurlxe*(tepltr*ev-te(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          elseif (ibctepr .eq. 2) then
             yldot(iv1) = nurlxe*(te(ixt1,iy)-te(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          endif 
        endif  # end loop for isteon=1

c  Do the ion temp Eqn ---------------------------
        if (istionxy(ixt,iy) == 1) then
          iv2 = idxti(ixt,iy)
          if (ibctipr == 1) then
c            totfeix is thermal + kinetic + viscous flux
            totfeixr(iy,jx) = feix(ixt1,iy) + ckinfl*kfeix 
            totfnix = 0.
            do ifld = 1, nfsp
            # the condition on zi must include use of v2cd which is not properly set to zero in pandf1
            if(zi(ifld)>1e-10) then
              totfeixr(iy,jx) = totfeixr(iy,jx) + cfeixdbo*(
     .                    2.5* ni(ixt,iy,ifld)*v2cd(ixt1,iy,ifld)*
     .                          sx(ixt1,iy)*rbfbt(ixt,iy)+
     .                           floxibgt(ixt1,iy,ifld) )*ti(ixt,iy)
               totfnix = totfnix + fnix(ixt1,iy,ifld)
               endif
            enddo
            if (isupgon(1)==1) then
c               Different boundary conditions for neutral momentum equation
                osmw = onesided_maxwellian( 
     .              tg(ixt1,iy,1), ng(ixt1,iy,1), mg(1), sx(ixt1,iy), tgmin*ev
     .          ) 
                yldot(iv2) = nurlxi*(totfeixr(iy,jx)
     .                      - totfnix*bcir(iy,jx)*ti(ixt,iy)
     .                      + cftiexclg*( 
     .                          - cfneut*fnix(ixt1,iy,iigsp)*bcen*tg(ixt,iy,1)
     .                          - (cgengpl*2.*tg(ixt,iy,1) 
     .                          - cgpld*eion*ev)*osmw
     .                          - cmneut*fnix(ixt1,iy,1)*recycp(1)
     .                              * cmntgpr*(ti(ixt,iy)-eidisspr*ev)
     .                      ) 
     .              ) / (vpnorm*ennorm*sx(ixt1,iy))
            else
                 yldot(iv2) = nurlxi*(totfeixr(iy,jx) 
     .                        -totfnix*bcir(iy,jx)*ti(ixt,iy)+
     .            cftiexclg*( -cmneut*fnix(ixt1,iy,1)*recycp(1)*
     .                                 cmntgpr*(ti(ixt,iy)-eidisspr*ev)    
     .                                   ) )/(vpnorm*ennorm*sx(ixt1,iy))
            endif # end test for isupgon(1)==1

          elseif (ibctipr .eq. 0) then
             yldot(iv2) = nurlxi*(tipltr*ev - ti(ixt,iy))
     .                                  * 1.5 * ne(ixt,iy)/ennorm
          elseif (ibctipr .eq. 2) then
             yldot(iv2) = nurlxi*(ti(ixt1,iy) - ti(ixt,iy))
     .                                  * 1.5 * ne(ixt,iy)/ennorm
          endif # end loop for bctipr
        endif  # end loop for istion=1

c       Next, the hydrogenic gas equations --
        do igsp = 1, nhgsp  # imp gas below
           if (isngonxy(ixt,iy,igsp) .eq. 1) then
             iv = idxg(ixt,iy,igsp)
             if (recyrb(iy,igsp,jx) .gt. 0.) then  # normal recycling
               flux_inc = fac2sp*fnix(ixt1,iy,1)
               if (ishymol.eq.1 .and. igsp.eq.2) then
                flxa= ismolcrm*(1-albrb(iy,1,jx))*onesided_maxwellian(
     .              tg(ixt1,iy,1), ng(ixt1,iy,1), mg(1), sx(ixt1,iy), tgmin*ev)
                 if (isupgon(1) .eq. 1) then  # two atoms for one molecule
                   flux_inc = 0.5*( fnix(ixt1,iy,1) +fnix(ixt1,iy,2)-flxa) 
                 else
                   flux_inc = 0.5*( fnix(ixt1,iy,1) +fngx(ixt1,iy,1)-flxa) 
                 endif
               endif
               osmw = onesided_maxwellian(
     .              tg(ixt1,iy,igsp), 1.0, mg(igsp),
     .              isoldalbarea*sx(ixt1,iy) + (1-isoldalbarea)*sxnp(ixt1,iy),
     .              tgmin*ev
     .         )
               yldot(iv) = nurlxg *  ( fngx(ixt1,iy,igsp)
     .                  + fngxrb_use(iy,igsp,jx)
     .                  - fngxsrb(iy,igsp,jx) 
     .                  + recyrb(iy,igsp,jx)*flux_inc
     .                  - (1-albrb(iy,igsp,jx))*osmw*ng(ixt1,iy,igsp)
     .              ) / (vpnorm*n0g(igsp)*sx(ixt1,iy))
             elseif( (recyrb(iy,igsp,jx) <=  0.) 
     .              .and. (recyrb(iy,igsp,jx) >= -1.) 
     .       ) then  # recyrb is albedo
                osmw = onesided_maxwellian(
     .              tg(ixt,iy,igsp), 1.0, mg(igsp), sx(ixt1,iy), tgmin*ev
     .          )
                yldot(iv) =  nurlxg*( fngx(ixt1,iy,igsp)
     .                  - (1+recyrb(iy,igsp,jx))*osmw*ng(ixt,iy,igsp)
     .                  ) / (osmw*n0g(igsp))
             elseif (recyrb(iy,igsp,jx) < -1. .and.
     .               recyrb(iy,igsp,jx) > -2.)  then  #fix density ngrfix
               yldot(iv)=nurlxg*(ngrfix - ng(ixt,iy,igsp))/n0g(igsp)
             elseif (recyrb(iy,igsp,jx) <= -2.) then   #zero gradient
               yldot(iv) = nurlxg*(ng(ixt1,iy,igsp) -
     .                                      ng(ixt,iy,igsp))/n0g(igsp)
             endif
             if (is1D_gbx.eq.1) yldot(iv) = nurlxg*(ng(ixt1,iy,igsp) -
     .                                      ng(ixt,iy,igsp))/n0g(igsp)
c   Special coding for Maxim
             if (isngrf==1) yldot(iv) = nurlxg*(ngrfix - 
     .                                    ng(ixt,iy,igsp))/n0g(igsp)
           endif # end if-test on isngon
        enddo # end do-loop on igsp



c       Finally, do impurity gas equations --
c       Set neutral gas impurity flux to be the sum of the fluxes of the
c       sputtered impurities plus recycled impurities from all charge states.
        if (isimpon .ge. 4 .and. isimpon .le. 7 .and. nzspt.ge.1) then
           nzsp_rt = nhsp
           do igsp = nhgsp+1, ngsp
	      jz = max(igsp - nhgsp, 1)   # identify impurity index
	      if (jz > 1) nzsp_rt = nzsp_rt + nzsp(jz-1) #prev index for fnix
              if (isngonxy(ixt,iy,igsp) .eq. 1) then
                 iv = idxg(ixt,iy,igsp)
                 hflux = 0.
                 do ihyd = 1, nhsp
                    if (zi(ihyd).gt.0.) hflux = hflux + fnix(ixt1,iy,ihyd)
                 enddo
                 zflux = 0.
                 do iimp = 1, nzsp(jz) # need to limit this loop per species
                   zflux = zflux + fnix(ixt1,iy,nzsp_rt+iimp)
                 enddo
                 sputflxrb(iy,igsp,jx) = 0.
                 zflux_chm = 0.
                 if (isph_sput(igsp) .ge. 1) then  # use fits for phys sput
                   do ifld = ipsputt_s, ipsputt_e 
                     eng_sput = ( 0.5*mi(ifld)*up(ixt1,iy,ifld)**2 + 
     .                             ti(ixt,iy) + zi(ifld)*
     .                             kappar(iy,jx)*te(ixt,iy) )/ev
                     if(zi(ifld)>0.) sputflxrb(iy,igsp,jx) = sputflxrb(iy,igsp,jx) +
     .                                             fnix(ixt1,iy,ifld)*
     .                            fphysyrb(igsp,jx)*yld96(matp,matt,eng_sput)
                   

                   enddo
                   if (isph_sput(igsp) .ge. 2) then  # add chem sput from ions
                      do ifld = 1,1  #(ipsputt_s, ipsputt_e) place holder for imp/imp sputt
                        eng_sput = ( 0.5*mi(ifld)*up(ixt1,iy,ifld)**2 + 
     .                               ti(ixt,iy) + zi(ifld)*
     .                               kappar(iy,jx)*te(ixt,iy) )/ev
                        flx_incid = abs(fnix(ixt1,iy,ifld))/sx(ixt1,iy)
                        call sputchem (isch_sput(igsp),eng_sput,tvplatrb(iy,jx),
     .                             flx_incid, yld_chm)
                        sputflxrb(iy,igsp,jx) = sputflxrb(iy,igsp,jx) + 
     .                      fchemyrb(igsp,jx)*fnix(ixt1,iy,ifld)*yld_chm
                      enddo
                   endif   #test on isph_sput >=2
                   if (isph_sput(igsp) .eq. 3) then # add chem sput from h. neut
                     do igsp2 = 1, 1+ishymol
                       t0p = max(tg(ixt1,iy,igsp2),tgmin*ev)
                       flx_incid=ng(ixt,iy,igsp2)*.25*
     .                                        sqrt(8*t0p/(pi*mg(igsp2)))
                       call sputchem (isch_sput(igsp),t0p/ev,tvplatrb(iy,jx),
     .                                               flx_incid, yld_chm)
                       zflux_chm = zflux_chm + flx_incid*
     .                             fchemyrb(igsp,jx)*yld_chm*sx(ixt1,iy)
                     enddo
                   endif    #test on isph_sput = 3
                 endif   #test on isph_sput > 1              
                 if (   (sputtrb(iy,igsp,jx) .ge. 0.)
     .                   .or. (abs(sputflxrb(iy,igsp,jx)).gt.0.)
     .            ) then
                    osmw = onesided_maxwellian(
     .                  cdifg(igsp)*tg(ixt1,iy,igsp), ng(ixt1,iy,igsp), mg(igsp),
     .                  isoldalbarea*sx(ixt1,iy) + (1-isoldalbarea)*sxnp(ixt1,iy),
     .                  tgmin*ev
     .              )
                    zflux = - sputtrb(iy,igsp,jx) * hflux  
     .                      - sputflxrb(iy,igsp,jx)
     .                      - recyrb(iy,igsp,jx) * zflux
     .                      + (1-albrb(iy,igsp,jx))*osmw
     .                      - zflux_chm 
     .                      + fngxsrb(iy,igsp,jx)
     .                      - fngxrb_use(iy,igsp,jx)
                    yldot(iv) = nurlxg * (fngx(ixt1,iy,igsp) - zflux) /
     .                         (n0(igsp) * vpnorm * sx(ixt1,iy))
                 elseif (sputtrb(iy,igsp,jx).ge.-9.9) then # neg. sputtrb ==> albedo
                    osmw = onesided_maxwellian(
     .                  cdifg(igsp)*tg(ixt1,iy,igsp), 1.0, mg(igsp),
     .                  isoldalbarea*sx(ixt1,iy) + (1-isoldalbarea)*sxnp(ixt1,iy),
     .                  tgmin*ev
     .              )
                    yldot(iv) =  nurlxg*( fngx(ixt1,iy,igsp) -
     .                  (1+sputtrb(iy,igsp,jx))*osmw*ng(ixt1,iy,igsp)
     .                  ) / (osmw*n0g(igsp))
                 else                                # sputtrb < -9.9 ==> fix dens
                    yldot(iv) = -nurlxg*(ng(ixt,iy,igsp)-ngplatrb(igsp,jx))/
     .                                                        n0g(igsp)    
                 endif  # end if-test on sputtrb
              endif  # end if-test on isngon
           enddo  # end do-loop in igsp
        endif  # end if-test on isimpon and nzspt

c ... Neutral temperature - test if tg eqn is on, then set BC
	do igsp = 1, ngsp
          if (istgonxy(ixt,iy,igsp) == 1) then
            iv = idxtg(ixt,iy,igsp)
            if (istgrb(igsp) == 0) then  #set tg=tgwall
              yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                 tg(ixt,iy,igsp))/(temp0*ev)
            elseif (istgrb(igsp) == 1) then  #extrapolate tg from interior
              tbound = tg(ixt1,iy,igsp) + gxf(ixt2,iy)*
     .               (tg(ixt1,iy,igsp)-tg(ixt2,iy,igsp))/gxf(ixt1,iy)
              tbound = max(tbound,0.5*temin*ev)
              yldot(iv) = nurlxg*(tbound - tg(ixt,iy,igsp))/(temp0*ev)
            elseif (istgrb(igsp) == 2)  #placeholder for gradient BC
        call xerrab("**INPUT ERROR: istgrb=2 grad opt not implemented")
            elseif (istgrb(igsp) == 3)  #Maxwell thermal flux to wall
                osmw = onesided_maxwellian( 
     .              cdifg(igsp)*tg(ixt1,iy,igsp), ng(ixt1,iy,igsp),
     .              mg(igsp), sx(ixt1,iy), tgmin*ev
     .          )
                yldot(iv) =  nurlxg*( fegx(ixt1,iy,igsp) 
     .                  - 2*cgengmpl*max(
     .                      cdifg(igsp)*tg(ixt1,iy,igsp), temin*ev
     .                  )*osmw) / (sx(ixt1,iy)*vpnorm*ennorm)
            elseif (istgrb(igsp) == 4) 
	      if (isupgon(igsp)==1) then
                if (recyrb(iy,igsp,jx) .gt. 0.) then
                  fng_alb = (1-albrb(iy,igsp,jx))*onesided_maxwellian(
     .                  tg(ixt1,iy,igsp), ng(ixt1,iy,igsp), mg(igsp),
     .                  sx(ixt1,iy), tgmin*ev)
                  yldot(iv) = nurlxg * ( fegx(ixt1,iy,igsp)
     .                                  -cfalbedo*fng_alb*max(tg(ixt1,iy,igsp),tgmin*ev)
     .                                  +recyrb(iy,igsp,jx)*(1.-cfdiss)
     .                                  *fnix(ixt1,iy,1)
     .                                  *recyce*cfalbedo
     .                                  *( kappar(iy,jx)*zi(1)*te(ixt,iy)
     .                                    +ti(ixt,iy) ) )/
     .                                  (vpnorm*ennorm*sx(ixt1,iy))
                elseif (recyrb(iy,igsp,jx) <=  0. .and.
     .                  recyrb(iy,igsp,jx) >= -1.) then   # recyrb is albedo
                  osmw = onesided_maxwellian(
     .                  tg(ixt1,iy,igsp), ng(ixt1,iy,igsp),
     .                  mg(igsp), sx(ixt1,iy), tgmin*ev
     .            )
                  yldot(iv) = nurlxg * ( fegx(ixt1,iy,igsp)
     .                      - cfalbedo*(1+recyrb(iy,igsp,jx))*osmw
     .                      * max(tg(ixt1,iy,igsp),tgmin*ev)
     .                  ) / (vpnorm*ennorm*sx(ixt1,iy))
                elseif (recyrb(iy,igsp,jx) < -1.) then  #..half Maxwellian
                  t0 = max(tg(ixt1,iy,igsp),tgmin*ev)
                  yldot(iv) = nurlxg * ( fegx(ixt1,iy,igsp)
     .                                  -cfalbedo*fnix(ixt1,iy,iigsp)*t0 )/
     .                                   (vpnorm*ennorm*sx(ixt1,iy))
                endif
              endif
            elseif (istgrb(igsp) == 5) then  #set tg=ti*cftgtipltr
              yldot(iv) = nurlxg*(ti(ixt,iy)*cftgtipltr(igsp) -
     .                                 tg(ixt,iy,igsp))/(temp0*ev)
            else
              call xerrab("**INPUT ERROR: istgrb set to unknown option")
            endif
          endif
        enddo

      enddo # end big do-loop on iy

c ...   Make electric field at plate uniform over last 2 cells --
        if (isexunif.eq.1) then
           do iy = j2, j5
	     if(isphionxy(ixt1,iy)==1) then
               iv = idxphi(ixt1,iy)  # next formula assumes uniform mesh
               yldot(iv)=-nurlxp*(3*phi(ixt1,iy)-2*phi(ixt,iy)-
     .                                            phi(ixt2,iy))/temp0
             endif
           enddo
        endif

      endif # end if-test on i6 and isfixrb

      endif # end if-test on xcnearrb and openbox
      enddo # end do-loop over nxpt mesh regions
c************************************************************************
c     end standard divertor plate conditions for right boundaries
c************************************************************************

        RETURN
        END SUBROUTINE right_boundary