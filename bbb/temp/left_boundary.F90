        SUBROUTINE left_boundary(neq, yl, yldot)
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
     .     ueb, nbound, tbound, ut0, sumb,
     .     fng_chem, vbound, eng_sput, flx_incid,
     .     yld_chm, t0p, zflux_chm, fqytotc, flux_inc,
     .     totfnex, totfnix, fqpsate, qpfac, arglgphi, faceel,
     .     faceel2, csfac, lambdae, 
     .     f_cgpld, 
     .     vxa, ta0, flxa
          integer ii,isphion2, nzsp_rt, jz
          real hflux, zflux
          integer ifld, ihyd, iimp, ix_fl_bc, ixc1, igsp2
          real fng_alb
          integer ixt, ixt1, ixt2, ixt3, jx, ixc, ierr
          integer ixtl, ixtl1, ixtr,ixtr1
          #Former Aux module variables
          integer ix,iy,igsp,iv,iv1,iv2,iv3,iv4,ix1,ix2,ix3,ix4,ixd
          real osmw
          real t0

          real yld96, kappa



c ====================================================================
c ======================== The ix = 0 boundary =======================
c ====================================================================

            ix = 0
            ixd = 1
c********************************************************************
c...  First check if ix=0 has fixed boundary values, no potential
c...  isfixlb=1 sets all profiles; isfixlb=2 sets reflection boundary 
c...  conditions
c********************************************************************
      if (i3 .le. 0 .and. isfixlb(1) .ne. 0) then
      do iy = j2, j5
         do ifld = 1 , nisp
            if(isnionxy(ix,iy,ifld) .eq. 1) then
               iv1 = idxn(ix,iy,ifld)
               yldot(iv1) = nurlxn *
     .                (nib(ifld)*nibprof(iy)-ni(ix,iy,ifld))/n0(ifld)
               if(isfixlb(1).eq.2) yldot(iv1) = nurlxn * (1/n0(ifld)) *
     .                               (ni(1,iy,ifld) - ni(ix,iy,ifld))
            endif
         enddo
         do ifld = 1, nusp
            if(isuponxy(ix,iy,ifld) .eq. 1) then
               iv2 = idxu(ix,iy,ifld)
               yldot(iv2) = nurlxu *
     .           (upb(ifld)*upbprof(iy) - up(ix,iy,ifld))/vpnorm
               if(isfixlb(1).eq.2) yldot(iv2) = nurlxu *
     .                           (0. - up(ix,iy,ifld))/vpnorm
               if(isfixlb(1).eq.2 .and. yylb(iy,1).gt.rlimiter) then
                  cs = sqrt( (te(ix,iy)+ti(ix,iy))/mi(ifld) )
                  yldot(iv2) = nurlxu*
     .                          (-cs -up(ix,iy,ifld))/vpnorm

               endif
            endif
         enddo


c...  now do the gas and temperatures
         if(isteonxy(ix,iy) .eq. 1) then
           iv1 = idxte(ix,iy)
           yldot(iv1) = nurlxe * ne(ix,iy) *
     .                     (teb*ev*tebprof(iy) - te(ix,iy))/ennorm
           if(isfixlb(1).eq.2) yldot(iv1) = nurlxe * ne(ix,iy) *
     .                               (te(ixd,iy) - te(ix,iy))/ennorm
           if(isfixlb(1).eq.2 .and. yylb(iy,1).gt.rlimiter) then
              yldot(iv1) = - nurlxe*(feex(ix,iy)/sx(ix,iy) - bcee*
     .                               ne(ix,iy)*vex(ix,iy)*te(ix,iy))/
     .                                 (vpnorm*ennorm)
           endif
         endif
         if(istionxy(ix,iy) .eq. 1) then
           iv2 = idxti(ix,iy)
            yldot(iv2) = nurlxi * ne(ix,iy) *
     .                     (tib*ev*tibprof(iy) - ti(ix,iy))/ennorm
            if(isfixlb(1).eq.2) yldot(iv2) = nurlxi * ne(ix,iy) *
     .                               (ti(ixd,iy) - ti(ix,iy))/ennorm
            if(isfixlb(1).eq.2 .and. yylb(iy,1).gt.rlimiter) then
               yldot(iv2) = -nurlxi*
     .          ( feix(ix,iy) - bcei*ti(ix,iy)*fac2sp*fnix(ix,iy,1) ) / 
     .                                       (vpnorm*ennorm*sx(ix,iy))
            endif
         endif
         do igsp = 1, ngsp
            if(isngonxy(ix,iy,igsp) .eq. 1) then
               iv = idxg(ix,iy,igsp)
               yldot(iv) = nurlxg * (ngbackg(igsp) - 
     .                                         ng(ix,iy,igsp))/n0g(igsp)
               if(isfixlb(1).eq.2) yldot(iv) = nurlxg * 
     .                        (ng(ixd,iy,igsp) - ng(ix,iy,igsp))/n0g(igsp)
               if(isfixlb(1).eq.2 .and. yylb(iy,1).gt.rlimiter) then
                  flux_inc = fac2sp*fnix(ix,iy,1)
                  if (ishymol.eq.1 .and. igsp.eq.2) then
                    ta0 = engbsr * max(tg(ixd,iy,1),temin*ev)
                    vxa = 0.25 * sqrt( 8*ta0/(pi*mg(1)) )
                    flxa = ismolcrm*(1-alblb(iy,1,1))*ng(ixd,iy,1)*vxa*sx(ix,iy)

                    if (isupgon(1) .eq. 1) then  # two atoms per molecule
                      flux_inc = 0.5*( fnix(ix,iy,1) + fnix(ix,iy,2) + flxa)
                    else
                      flux_inc = 0.5*( fnix(ix,iy,1) + fngx(ix,iy,1) + flxa) 
                    endif
                  endif
                  osmw = onesided_maxwellian(
     .                  tg(ixd,iy,igsp), ng(ixd,iy,igsp), mg(1), 
     .                  isoldalbarea*sx(ix,iy) + (1-isoldalbarea)*sxnp(ix,iy),
     .                  temin*ev
     .            )
                  yldot(iv) = -nurlxg * ( fngx(ix,iy,igsp)
     .                      - fngxlb_use(iy,igsp,1)
     .                      + fngxslb(iy,igsp,1) 
     .                      + recylb(iy,igsp,1)*flux_inc
     .                      + (1-alblb(iy,igsp,1))*osmw
     .                  ) / (vpnorm*n0g(igsp)*sx(ix,iy))
               endif
               if (is1D_gbx.eq.1) yldot(iv) = nurlxg*(ng(ixd,iy,igsp) -
     .                                    ng(ix,iy,igsp))/n0g(igsp)
            endif
         enddo
c ... Neutral temperature - test if tg eqn is on, then set BC
        do igsp = 1, ngsp
           if (istgonxy(ix,iy,igsp) == 1) then
             iv = idxtg(ix,iy,igsp)
             yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(ix,iy,igsp))/(temp0*ev)
             if(isfixlb(1)==2) then #just above applies if isfixlb=1
               yldot(iv)=nurlxg*(tg(ixd,iy,igsp)-tg(ix,iy,igsp))/(temp0*ev)
             endif
             if(isfixlb(1)==2 .and. yylb(iy,1) > rlimiter) then
               yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(ix,iy,igsp))/(temp0*ev)
             endif
           endif
         enddo

         if (isphionxy(ix,iy) .eq. 1) then
            iv = idxphi(ix,iy)
            yldot(iv) = nurlxp*(phi(ixd,iy) - phi(ix,iy))/temp0
         endif

      enddo

      endif         # end of ix = 0, isfixlb.ne.0 boundary conditions

c...  If isfixlb=2, check if i2,i5 range for yldot loop in pandf includes ixpt2(1)
c...  Then overwrite pandf value is up(ixpt2(1)) --> 0 eqn.
      if (isfixlb(1) .eq. 2) then
         if (i2.le.ixpt2(1) .and. i5.ge.ixpt2(1) .and. j2.le.iysptrx2(1)) then  
            do ifld = 1, nusp
               do iy = 0+1-iymnbcl, iysptrx2(1)
                 if(isuponxy(ixpt2(1),iy,ifld)==1) then
                     iv = idxu(ixpt2(1),iy,ifld)
                     yldot(iv) = nurlxu*(0.-up(ixpt2(1),iy,ifld))/vpnorm
                 endif
               enddo
            enddo
         endif 
      endif            # end of isfixlb=2, check for ix=ixpt2
c***********************************************************************

c************************************************************************
c************************************************************************
c     begin standard divertor plate conditions for left boundaries 
c************************************************************************
c************************************************************************
      do jx = 1, nxpt # loop over nxpt mesh regions
      if (xcnearlb .or. openbox) then
c     First, the density equations --
        do ifld = 1, nisp
          if ( i3.le.(ixlb(jx)+isextrnp) .and. isfixlb(jx)==0) then
          do iy = j2, j5
       
            ixt  = ixlb(jx)        # analog of ix=0
            ixt1 = ixp1(ixt,iy)    # analog of ix=1
            ixt2 = ixp1(ixt1,iy)   # analog of ix=2
            ixt3 = ixp1(ixt1,iy)   # analog of ix=2

            if(isnionxy(ixt,iy,ifld)==1) then
              iv1 = idxn(ixt,iy,ifld)
              if (isupgon(1)==1 .and. zi(ifld)==0.0) then   ## neutrals
                if (recylb(iy,1,jx) .gt. 0.) then           # recycling
                  osmw = onesided_maxwellian(
     .              tg(ixt1,iy,1), ni(ixt1,iy,ifld), mi(ifld),
     .              isoldalbarea*sx(ixt,iy) + (1-isoldalbarea)*sxnp(ixt,iy),
     .              tgmin*ev
     .            )
                  yldot(iv1) = -nurlxg * (fnix(ixt,iy,ifld) 
     .                      + recylb(iy,1,jx)*fnix(ixt,iy,1) 
     .                      - fngxlb_use(iy,1,jx)
     .                      + (1-alblb(iy,1,jx))*osmw
     .                      - fngxslb(iy,1,jx) 
     .                  ) / (vpnorm*n0(ifld)*sx(ixt,iy))
                elseif( (recylb(iy,1,jx) <=  0.)
     .              .and. (recylb(iy,1,jx) >= -1.) 
     .          ) then  # recylb is albedo
                  osmw = onesided_maxwellian(
     .                  tg(ixt,iy,1), ni(ixt,iy,ifld), mi(1), 
     .                  sx(ixt,iy), tgmin*ev
     .            )
                  yldot(iv1) = -nurlxg * ( fnix(ixt,iy,ifld) 
     .                  + (1+recylb(iy,1,jx))*osmw
     .              ) / (vpnorm*n0(ifld)*sx(ixt,iy))
                elseif (recylb(iy,1,jx) < -1. .and.
     .                 recylb(iy,1,jx) > -2. ) then  #fix density nglfix
                  yldot(iv1)=nurlxg*(nglfix - ni(ixt,iy,ifld))/n0(ifld)
                elseif (recylb(iy,1,jx) <= -2.) then  #zero gradient
                  yldot(iv1) =nurlxn*(ni(ixt1,iy,ifld)-ni(ixt,iy,ifld))/
     .                                                          n0(ifld)
                endif # end if-test on recylb
              else                                     ## ions
                if (isextrnp==0) then                  # zero x-gradient
                   yldot(iv1) = nurlxn*(ni(ixt1,iy,ifld)-ni(ixt,iy,ifld))/
     .                                                          n0(ifld)
                else                                   # extrapolation
                   nbound =  ni(ixt1,iy,ifld) - gxf(ixt1,iy)*
     .                   (ni(ixt2,iy,ifld)-ni(ixt1,iy,ifld))/gxf(ixt,iy)
                   nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                    ni(ixt,iy,ifld)-1) ) ) + 0.2*ni(ixt,iy,ifld)
                   yldot(iv1) = nurlxn*(nbound-ni(ixt,iy,ifld))/n0(ifld)
                endif # end if-test on isextrnp
              endif # end if-test on isupgon and zi
            endif   # end if-test on isnionxy
          enddo     # end do-loop on iy
          endif     # end if-test on i3 and isfixlb and isnion
        enddo       # end do-loop on ifld

c     Now do the parallel velocity and other variables --
      if ( (i3 .le. ixlb(jx)) .and. isfixlb(jx)==0 ) then
      do iy = j2, j5 # begin big do-loop on iy

        ixt  = ixlb(jx)        # analog of ix=0
        ixt1 = ixp1(ixt,iy)    # analog of ix=1
        ixt2 = ixp1(ixt1,iy)   # analog of ix=2


        kfeix = 0.
        sumb = 0.
        do ifld = 1, nfsp  # set up generalized Bohm condition
          ueb = cfueb*( cf2ef*v2ce(ixt,iy,ifld)*rbfbt(ixt,iy) -
     .            vytan(ixt,iy,ifld) ) / rrv(ixt,iy) 
          sumb = sumb + ni(ixt,iy,ifld)*zi(ifld)**2*te(ixt,iy) /
     .                (mi(ifld)*(upi(ixt,iy,ifld)+ueb)**2 - ti(ixt,iy))
        enddo # end do-loop on ifld for Bohm condition
        sumb = sqrt(abs(sumb/ne(ixt,iy)))



c       Next, the momentum equations --
        do ifld = 1, nusp
          if (isuponxy(ixt,iy,ifld)==1) then
            iv2 = idxu(ixt,iy,ifld)
            cs = csfaclb(ifld,jx)*sqrt( (te(ixt,iy)+
     .                                  csfacti*ti(ixt,iy))/mi(ifld) )
            if (isupgon(1)==1 .and. zi(ifld)==0.0) then  ## neutrals
              if (recycmlb(iy,1,jx) > -9.9) then  # backscatter with recycm
                yldot(iv2) = -nurlxu*(recycmlb(iy,1,jx)*up(ixt,iy,1) + 
     .                                       up(ixt,iy,ifld))/vpnorm
              elseif (recycmlb(iy,1,jx) <= -9.9 .and. 
     .                    recycmlb(iy,1,jx) > -10.1) then # zero x-gradient
                yldot(iv2) = nurlxu*(up(ixt1,iy,ifld) -
     .                                          up(ixt,iy,ifld))/vpnorm
          else  # neutral thermal flux to wall if recycm < -10.1
                osmw = onesided_maxwellian( 
     .                  tg(ixt1,iy,1), 0.5*(nm(ixt1,iy,ifld)+nm(ixt,iy,ifld)),
     .                  mi(ifld), sx(ixt,iy), tgmin*ev
     .          )
                yldot(iv2) = -nurlxu*( fmix(ixt1,iy,ifld) 
     .              + up(ixt,iy,ifld)*cgmompl*osmw
     .              ) / (vpnorm*fnorm(ifld)*sx(ixt,iy))
              endif
            else                                         ## ions
              ueb = cfueb*( cf2ef*v2ce(ixt,iy,ifld)*rbfbt(ixt,iy) -
     .                vytan(ixt,iy,ifld) ) / rrv(ixt,iy) 
              yldot(iv2) = -nurlxu*(sumb - 1.)  # multispecies Bohm
              if (isbohmms==0) then          # simple Bohm condition
                yldot(iv2) = nurlxu * (-cs-ueb-up(ixt,iy,ifld))/vpnorm
              endif
              if(isupss(ifld)==1 .and. up(ixt1,iy,ifld)+ueb .lt. -cs)
                                             # dup/dx=0 if supersonic
     .          yldot(iv2) = nurlxu*(up(ixt1,iy,ifld)-up(ixt,iy,ifld))/
     .                                                           vpnorm
              if (isupss(ifld)==-1) then     # slip boundary conditions
                yldot(iv2) = nurlxu*(up(ixt1,iy,ifld)-up(ixt,iy,ifld))/
     .                                                           vpnorm
              elseif (isupss(ifld)==-2) then # extrap. + no pos. uu
                vbound = up(ixt1,iy,ifld) - gx(ixt2,iy)*
     .                   (up(ixt2,iy,ifld)-up(ixt1,iy,ifld))/gx(ixt1,iy) 
                vbound = min(vbound, -ueb)   # forces uu & fnix =< 0
                yldot(iv2) = nurlxu*(vbound - up(ixt,iy,ifld))/vpnorm
              elseif (isupss(ifld)==-3) then # modified Bohm condition
                vbound = -ueb -2*cs*uu(ixt1,iy,ifld)/(uu(ixt1,iy,ifld)-cs)
                vbound = min(vbound, -ueb)   # forces uu & fnix =< 0
                yldot(iv2) = nurlxu * (vbound-up(ixt,iy,ifld))/vpnorm
              endif # end if-test on isupss
            endif # end if-test on isupgon
          endif # end if-test on isuponxy
          






          if (zi(ifld)==0.0) then 
            kfeix = kfeix - cftiexclg*cfvcsx(ifld)*0.5*sx(ixt,iy)
     .                     *visx(ixt1,iy,ifld)*gx(ixt1,iy)
     .           *( up(ixt1,iy,ifld)**2 - up(ixt,iy,ifld)**2 )
          else
            kfeix = kfeix - cfvcsx(ifld)*0.5*sx(ixt,iy)
     .                     *visx(ixt1,iy,ifld)*gx(ixt1,iy)
     .           *( up(ixt1,iy,ifld)**2 - up(ixt,iy,ifld)**2 ) 
          endif
        enddo # end do-loop on ifld

c       Next, the potential equation --
        if(isphion+isphiofft .eq. 1) then
cc           fqpsatlb(iy,jx) = - qe*isfdiax*( ne(ixt,iy)*v2ce(ixt,iy,1)*
cc     .                              rbfbt(ixt,iy)*sx(ixt,iy) + fdiaxlb(iy,jx) )
cc           do ifld = 1, nfsp   # note fqp,fqpsat are poloidal proj. of || curr
cc              fqpsatlb(iy,jx) = fqpsatlb(iy,jx) - qe*zi(ifld)*ni(ixt,iy,ifld)*
cc     .                                   upi(ixt,iy,ifld)*sx(ixt,iy)*rrv(ixt,iy)
cc           enddo
           lambdae = 2e16*(te(ixt,iy)/ev)**2/ne(ixt,iy)  #approx mfp
           kincorlb(iy,jx) = 1./(1 + cfkincor*(lambdae/lcone(ixt,iy))*
     .                                       abs(ev*phi(ixt,iy)/te(ixt,iy)))
           fqpsate = qe*ne(ixt,iy)*sqrt(te(ixt,iy)/(2*pi*me))*
     .                                    kincorlb(iy,jx)*sx(ixt,iy)*rrv(ixt,iy)
c          NOTE: by definition, fqpsate is always > 0
           if (ikapmod==0) then
cc              if (fqp(ixt,iy) < 0.) then #limit to saturation current
cc                fqp(ixt,iy)=-( abs((fqp(ixt,iy)*fqpsatlb(iy,jx)))**exjbdry/
cc     .                      (abs(fqp(ixt,iy))**exjbdry +
cc     .                       abs(fqpsatlb(iy,jx))**exjbdry) )**(1/exjbdry)
cc              endif
              if ( fqpsatlb(iy,jx)+(1.-gamsec)*fqp(ixt,iy) > 0) then  # force +ve log argument
                 arglgphi=( ((fqpsatlb(iy,jx)+(1.-gamsec)*fqp(ixt,iy))
     .                          /fqpsate)**2 + exp(-kappamx)**2 )**0.5
              else
                 arglgphi = exp(-kappamx)
              endif
              if (iskaplex .eq. 0) kappal(iy,jx) = - log(arglgphi)
              if (newbcl(jx).eq.0 .and. iskaplex.eq.0) kappal(iy,jx) = 3.0
           elseif (ikapmod==1) then
              kappal(iy,jx) = kappa(fqpsatlb(iy,jx), fqpsate, 
     .                                         -(1.-gamsec)*fqp(ixt,iy))
           endif

           if (isphionxy(ixt,iy) .eq. 1) then
              iv = idxphi(ixt,iy)
              if ((isnewpot==100) .and. ((iy==1) .or. (iy==ny))) then
                 continue
              elseif (isphilbc==1) then
                 yldot(iv) = - nurlxp * (phi(ixt,iy) - phi0l(iy,jx))/temp0 
              else
                 yldot(iv) = -nurlxp*(1.-bctype(iy)) *
     .            (phi(ixt,iy)-kappal(iy,jx)*te(ixt,iy)/ev-phi0l(iy,jx))/temp0
     .                    -nurlxp * bctype(iy) * (1.-gamsec)*fqp(ixt,iy) 
     .                                         / (fqpsatlb(iy,jx)+cutlo)
              endif
           endif
	else    # corresponding to (isphion+isphiofft .ne. 1)
           fqpsate = 0.
	   kappal(iy,jx) = 3.  # some default of e*phi/Te - used for sputtering
        endif
        isphion2 = isphion + isphiofft
        bcel(iy,jx) = (1-newbcl(jx)*isphion2) * bcee
     .              + newbcl(jx)*isphion2 * (2. + kappal(iy,jx))
        if (iskaplex .eq. 1) bcel(iy,jx) = (2. + kappal(iy,jx))
        bcil(iy,jx) = (1-newbcl(jx)*isphion2) * bcei
     .              + newbcl(jx)*isphion2 * (2.5)

        t0 = te(ixt,iy)/ev
        f_cgpld = .5*(1. - cos(pi*(t0 - temin)/(.3 - temin)))
                                # ramp cgpld effect down at low guard cell
                                # Te; rampdown might not be necessary if
                                # guard cell were allowed to access negative
                                # temperatures (ETM 7 Jan 2014)
        if(t0 < temin) f_cgpld = 0.
        if(t0 > 0.3) f_cgpld = 1.

c   Do the electron temp Eqn -----------------------------------
        if (isteonxy(ixt,iy) == 1) then
      iv1 = idxte(ixt,iy)
          if (ibctepl == 1) then
            faceel =  bcel(iy,jx)*(fqpsate/qe)*exp(-kappal(iy,jx))
            faceel2 = bcel(iy,jx)*(fqpsate/qe)*exp(-kappamx+2)
            totfeexl(iy,jx) = feex(ixt,iy) + cfeexdbo*( 
     .             2.5*fdiaxlb(iy,jx) + floxebgt(ixt,iy) )*te(ixt,iy)
        totfnex = ne(ixt,iy)*vex(ixt,iy)*sx(ixt,iy)
            if (isphion+isphiofft==1) then
              yldot(iv1) = -nurlxe*(totfeexl(iy,jx)
     .                  +faceel*te(ixt,iy)
     .                  +faceel2*(te(ixt,iy)-te(ixt1,iy))
     .                  -cmneut*fnix(ixt,iy,1)*recycp(1)*eedisspl*ev
     .              ) / (sx(ixt,iy)*vpnorm*ennorm)
            else
              osmw = onesided_maxwellian(
     .              tg(ixt1,iy,1), ng(ixt1,iy,1), mg(1), 
     .              sx(ixt,iy), tgmin*ev
     .        )
              yldot(iv1) = -nurlxe*(totfeexl(iy,jx) 
     .                  - totfnex*te(ixt,iy)*bcel(iy,jx)
     .                  + cgpld*f_cgpld*osmw*0.5*ediss*ev
     .                  - cmneut*fnix(ixt,iy,1)*recycp(1)*eedisspl*ev
     .              ) / (sx(ixt,iy)*vpnorm*ennorm)
            endif
          elseif (ibctepl .eq. 0) then
             yldot(iv1) = nurlxe*(tepltl*ev-te(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          elseif (ibctepl .eq. 2) then
             yldot(iv1) = nurlxe*(te(ixt1,iy)-te(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          endif # end loop for ibctepl
        endif  # end loop for isteon=1

c  Do the ion temp Eqn ---------------------------
        if (istionxy(ixt,iy) == 1) then
          iv2 = idxti(ixt,iy)
          if (ibctipl == 1) then
c            totfeix is thermal + kinetic + viscous flux
            totfeixl(iy,jx) = feix(ixt,iy) + ckinfl*kfeix 
            totfnix = 0.
            do ifld = 1, nfsp
            # the condition on zi must include use of v2cd which is not properly set to zero in pandf1
            if(zi(ifld)>1e-10) then
              totfeixl(iy,jx) = totfeixl(iy,jx) + cfeixdbo*(
     .                      2.5* ni(ixt,iy,ifld)*v2cd(ixt,iy,ifld)*
     .                                sx(ixt,iy)*rbfbt(ixt,iy)+
     .                             floxibgt(ixt,iy,ifld) )*ti(ixt,iy)
              totfnix = totfnix + fnix(ixt,iy,ifld)
              endif
            enddo
            if (isupgon(1)==1) then
c              Different boundary conditions for neutral momentum equation
                osmw = onesided_maxwellian(
     .                  tg(ixt1,iy,1), ng(ixt1,iy,1), mg(1), sx(ixt,iy), tgmin*ev
     .          )
                yldot(iv2) = -nurlxi*(totfeixl(iy,jx)
     .                      - totfnix*ti(ixt,iy)*bcil(iy,jx)
     .                      + cftiexclg*( 
     .                        - cfneut*fnix(ixt,iy,iigsp)*tg(ixt,iy,1)*bcen
     .                        + (cgengpl*2.*tg(ixt,iy,1) 
     .                        - cgpld*eion*ev)*f_cgpld*osmw
     .                        - cmneut*fnix(ixt,iy,1)*recycp(1)
     .                          *cmntgpl*(ti(ixt,iy)-eidisspl*ev)
     .                      ) 
     .              ) / (vpnorm*ennorm*sx(ixt,iy))
            else
                 yldot(iv2) = -nurlxi*(totfeixl(iy,jx) 
     .                         -totfnix*bcil(iy,jx)*ti(ixt,iy)+
     .             cftiexclg*( -cmneut*fnix(ixt,iy,1)*recycp(1)*
     .                                cmntgpl*(ti(ixt,iy)-eidisspl*ev)
     .                                  ) )/(vpnorm*ennorm*sx(ixt,iy))
            endif #end loop for isupgon(1)

          elseif (ibctipl .eq. 0) then
             yldot(iv2) = nurlxi*(tipltl*ev - ti(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          elseif (ibctipl .eq. 2) then
             yldot(iv2) = nurlxi*(ti(ixt1,iy) - ti(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          endif # end loop for ibctipl
        endif  # end loop for istion=1

c       Do hydrogenic gas equations --
        do igsp = 1, nhgsp  # imp gas below
           if (isngonxy(ixt,iy,igsp) .eq. 1) then
             iv = idxg(ixt,iy,igsp)
             if (recylb(iy,igsp,jx) .gt. 0.) then  # normal recycling
               flux_inc = fac2sp*fnix(ixt,iy,1)
               if (ishymol.eq.1 .and. igsp.eq.2) then
                 ta0 = max(tg(ixt1,iy,1), temin*ev)
                 vxa = 0.25 * sqrt( 8*ta0/(pi*mg(1)) )
                 flxa = ismolcrm*(1-alblb(iy,1,jx))*ng(ixt1,iy,1)*vxa*sx(ixt,iy)

                 if (isupgon(1) .eq. 1) then  # two atoms for one molecule
                   flux_inc = 0.5*( fnix(ixt,iy,1) + fnix(ixt,iy,2) +flxa) 
                 else
                   flux_inc = 0.5*( fnix(ixt,iy,1) + fngx(ixt,iy,1) +flxa) 
                 endif
               endif
               osmw = onesided_maxwellian(
     .                  tg(ixt1,iy,igsp), ng(ixt1,iy,igsp), mg(igsp),
     .                  isoldalbarea*sx(ixt,iy) + (1-isoldalbarea)*sxnp(ixt,iy),
     .                  tgmin*ev
     .         )
               yldot(iv) = -nurlxg * ( fngx(ixt,iy,igsp)  
     .                  - fngxlb_use(iy,igsp,jx)
     .                  - fngxslb(iy,igsp,jx) 
     .                  + recylb(iy,igsp,jx)*flux_inc
     .                  + (1-alblb(iy,igsp,jx))*osmw
     .              ) / (vpnorm*n0g(igsp)*sx(ixt,iy))
             elseif ( (recylb(iy,igsp,jx) <=  0.)
     .           .and. (recylb(iy,igsp,jx) >= -1.) 
     .      ) then # recylb is albedo
                osmw = onesided_maxwellian(
     .              tg(ixt,iy,igsp), 1.0, mg(igsp), sx(ixt,iy), tgmin*ev
     .          )
                yldot(iv) = -nurlxg*( fngx(ixt,iy,igsp) +
     .              (1+recylb(iy,igsp,jx))*osmw*ng(ixt,iy,igsp)
     .              ) / (osmw*n0g(igsp))
             elseif (recylb(iy,igsp,jx) < -1. .and.
     .               recylb(iy,igsp,jx) > -2.) then  #fix density nglfix
               yldot(iv)=nurlxg*(nglfix - ng(ixt,iy,igsp))/n0g(igsp)
             elseif (recylb(iy,igsp,jx) <= -2.) then #zero gradient
               yldot(iv1) = nurlxn*(ni(ixt1,iy,ifld)-ni(ixt,iy,ifld))/
     .                                                         n0g(ifld)
             endif
             if (is1D_gbx.eq.1) yldot(iv) = nurlxg*(ng(ixt1,iy,igsp) -
     .                                      ng(ixt,iy,igsp))/n0g(igsp)
c   Special coding for Maxim
             if (isnglf==1) yldot(iv) = nurlxg*(nglfix - 
     .                                      ng(ixt,iy,igsp))/n0g(igsp)
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
                    if (zi(ihyd).gt.0.) hflux = hflux + fnix(ixt,iy,ihyd)
                 enddo
                 zflux = 0.
                 do iimp = 1, nzsp(jz) # loop limited by numb. species
                    zflux = zflux + fnix(ixt,iy,nzsp_rt+iimp)
                 enddo
                 sputflxlb(iy,igsp,jx) = 0.
                 zflux_chm = 0.
                 if (isph_sput(igsp) .ge. 1) then  # use fits for phys sput
                   do ifld = ipsputt_s, ipsputt_e 
                     eng_sput = ( 0.5*mi(ifld)*up(ixt,iy,ifld)**2 + 
     .                            ti(ixt,iy) +  zi(ifld)*
     .                            kappal(iy,jx)*te(ixt,iy) )/ev
                     if(zi(ifld)>0.) sputflxlb(iy,igsp,jx) = 
     .                                          sputflxlb(iy,igsp,jx) +
     .                                           fnix(ixt,iy,ifld)*
     .                                           fphysylb(igsp,jx)*
     .                                        yld96(matp,matt,eng_sput)
                   enddo
                   if (isph_sput(igsp) .ge. 2) then  # add chem sput from ions
                     do ifld = 1,1  #(ipsputt_s, ipsputt_e) place holder for imp/imp sputt
                       eng_sput = ( 0.5*mi(ifld)*up(ixt,iy,ifld)**2 + 
     .                            ti(ixt,iy) +  zi(ifld)*
     .                            kappal(iy,jx)*te(ixt,iy) )/ev
                       flx_incid = abs(fnix(ixt,iy,ifld))/sx(ixt,iy)
                       call sputchem (isch_sput(igsp),eng_sput,tvplatlb(iy,jx),
     .                             flx_incid, yld_chm)
                       sputflxlb(iy,igsp,jx) = sputflxlb(iy,igsp,jx) +
     .                       fchemylb(igsp,jx)*fnix(ixt,iy,ifld)*yld_chm
                     enddo
                   endif
                   if (isph_sput(igsp) .eq. 3) then #add chem sput from h neuts
                     do igsp2 = 1, 1+ishymol  #hydrogen neut fluxes only
                       t0p = max(tg(ixt1,iy,igsp2),tgmin*ev)
                       flx_incid = ng(ixt,iy,igsp2)*.25*
     .                                        sqrt(8*t0p/(pi*mg(igsp2)))
                       call sputchem (isch_sput(igsp),t0p/ev,tvplatlb(iy,jx),
     .                                               flx_incid, yld_chm)
                       zflux_chm = zflux_chm - flx_incid*
     .                              fchemylb(igsp,jx)*yld_chm*sx(ixt,iy)
                     enddo
                   endif
                 endif
                 osmw = onesided_maxwellian(
     .                  cdifg(igsp)*tg(ixt1,iy,igsp), 1.0, mg(igsp),
     .                  isoldalbarea*sx(ixt1,iy) + (1-isoldalbarea)*sxnp(ixt1,iy),
     .                  tgmin*ev
     .           )
                 if ( (sputtlb(iy,igsp,jx) .ge. 0.)
     .                   .or. (abs(sputflxlb(iy,igsp,jx)).gt. 0.)
     .           ) then
                    zflux = - sputtlb(iy,igsp,jx) * hflux 
     .                      - sputflxlb(iy,igsp,jx)
     .                      - recylb(iy,igsp,jx) * zflux
     .                      - (1-alblb(iy,igsp,jx))*osmw*ng(ixt1,iy,igsp)
     .                      - zflux_chm 
     .                      + fngxslb(iy,igsp,jx)
     .                      + fngxlb_use(iy,igsp,jx)
                    yldot(iv) = -nurlxg * (fngx(ixt,iy,igsp) 
     .                  - zflux) / (n0(igsp) * vpnorm * sx(ixt,iy))
                 elseif (sputtlb(iy,igsp,jx).ge.-9.9) then # neg. sputtlb ==> albedo
                    yldot(iv) = -nurlxg*( fngx(ixt,iy,igsp)
     .                  - (1+sputtlb(iy,igsp,jx))*osmw*ng(ixt,iy,igsp)
     .                  ) / (osmw*n0g(igsp))
                 




                 else                                # sputtlb < -9.9 ==> fix dens
                    yldot(iv) = -nurlxg*(ng(ixt,iy,igsp)-ngplatlb(igsp,jx))/
     .                                                        n0g(igsp)
                 endif  # end if-test on sputtlb
              endif  # end if-test on isngon
           enddo  # end do-loop in igsp
        endif  # end if-test on isimpon and nzspt

c ... Neutral temperature - test if tg eqn is on, then set BC
	do igsp = 1, ngsp
          if (istgonxy(ixt,iy,igsp) == 1) then
            iv = idxtg(ixt,iy,igsp)
            if (istglb(igsp) == 0) then  #set tg=tgwall
              yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(ixt,iy,igsp))/(temp0*ev)
            elseif (istglb(igsp) == 1)  #extrap. tg from interior cells
              tbound = tg(ixt1,iy,igsp) - gyf(ixt1,iy)*
     .                   (tg(ixt2,iy,igsp)-tg(ixt1,iy,igsp))/gxf(ixt,iy)
              tbound = max(tbound,0.5*temin*ev)
              yldot(iv) = nurlxg*(tbound - tg(ixt,iy,igsp))/(temp0*ev)
            elseif (istglb(igsp) == 2)  #placeholder for gradient BC
         call xerrab("**INPUT ERROR: istglb=2 grad opt not implemented")
            elseif (istglb(igsp) == 3)  #Maxwell thermal flux to wall
              osmw =  onesided_maxwellian( 
     .              cdifg(igsp)*tg(ixt1,iy,igsp), ng(ixt1,iy,igsp), 
     .              mg(igsp), sx(ixt,iy), tgmin*ev
     .        )
              yldot(iv) =  -nurlxg*( fegx(ixt,iy,igsp) 
     .                  + 2*cgengmpl*max(
     .                      cdifg(igsp)*tg(ixt1,iy,igsp), tgmin*ev
     .                  )*osmw ) / (sx(ixt,iy)*vpnorm*ennorm)
            elseif (istglb(igsp) == 4) 
	      if (isupgon(igsp)==1) then
                if (recylb(iy,igsp,jx) .gt. 0.) then
                  fng_alb=(1-alblb(iy,igsp,jx))*onesided_maxwellian(
     .                      tg(ixt1,iy,igsp), ng(ixt1,iy,igsp), 
     .                      mg(igsp), sx(ixt,iy), tgmin*ev
     .            )
                  yldot(iv) = -nurlxg*( fegx(ixt,iy,igsp)
     .                                 +cfalbedo*fng_alb
     .                                 * max(tg(ixt1,iy,igsp),tgmin*ev)
     .                                 +recylb(iy,igsp,jx)*(1.-cfdiss)
     .                                 *fnix(ixt,iy,1)
     .                                 *recyce*cfalbedo
     .                                 *( kappal(iy,jx)*zi(1)*te(ixt,iy)
     .                                   +ti(ixt,iy) ) )/
     .                                 (vpnorm*ennorm*sx(ixt,iy))
                elseif( (recylb(iy,igsp,jx) <=  0.)
     .              .and. (recylb(iy,igsp,jx) >= -1.) 
     .          ) then  # recylb is albedo
                  osmw = onesided_maxwellian(
     .                  tg(ixt1,iy,igsp), ng(ixt1,iy,igsp), 
     .                  mg(igsp), sx(ixt,iy), tgmin*ev
     .            )
                  yldot(iv) = -nurlxg*( fegx(ixt,iy,igsp)
     .                      + cfalbedo*(1-recylb(iy,igsp,jx))*osmw
     .                      * max(tg(ixt1,iy,igsp),tgmin*ev) 
     .                  ) / (vpnorm*ennorm*sx(ixt,iy))
                elseif (recylb(iy,igsp,jx) < -1.) then  #..half Maxwellian
                  yldot(iv) = -nurlxg*( 
     .                  fegx(ixt,iy,igsp) + cfalbedo*fnix(ixt,iy,iigsp)
     .                  * max(tg(ixt1,iy,igsp),tgmin*ev)
     .              ) / (vpnorm*ennorm*sx(ixt,iy))
                endif
              endif
            elseif (istglb(igsp) == 5) then  #set tg=ti*cftgtipltl
              yldot(iv) = nurlxg*(ti(ixt,iy)*cftgtipltl(igsp) -
     .                                    tg(ixt,iy,igsp))/(temp0*ev)
            else
              call xerrab("**INPUT ERROR: istglb set to unknown option")
            endif
          endif
        enddo

      enddo # end big do-loop on iy

c ...   Make electric field at plate uniform over last 2 cells --
        do iy = j2, j5
           if (isphionxy(ixt1,iy)==1 .and. isexunif==1) then
              iv = idxphi(ixt1,iy)  # next formula assumes uniform mesh
              yldot(iv)=-nurlxp*(3*phi(ixt1,iy)-2*phi(ixt,iy)-
     .                                            phi(ixt2,iy))/temp0
           endif
        enddo

      endif # end if-test on i3 and isfixlb

      endif # end if-test on xcnearlb and openbox
      enddo # end do-loop over nxpt mesh regions
c************************************************************************
c     end standard divertor plate conditions for left boundaries
c************************************************************************

        RETURN
        END SUBROUTINE left_boundary