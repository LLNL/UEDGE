        SUBROUTINE north_boundary(neq, yl, yldot)
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


c ====================================================================
c ======================== The iy=ny+1 boundary ======================
c ====================================================================









c...  now do the iy = ny+1 boundary
c...  if extrapolation b.c.on outer wall, isextrw=1, otherwise isextrw=0
      if (j7 .ge. (ny+1)-isextrnw .or. j7 .ge. (ny+1)-isextrtw) then
      


      do 280 ifld = 1 , nisp
        do 278 ix = i4+1-ixmnbcl, i8-1+ixmxbcl
          if (isnionxy(ix,ny+1,ifld)==1) then
            iv1 = idxn(ix,ny+1,ifld)


c ...  First check if this is neutral species (zi=0) & do BC
c ---  Different boundary conditions if neutral momentum equation use;
c ---  typically hydrogen only, so DIVIMP chem sputt not used here
            if (isupgon(1) .eq. 1 .and. zi(ifld) .eq. 0.0) then






























                if (iscplo(ix) .eq. 1) call wsmodo(1)
                fng_chem = 0.
                do ii = 1, ngsp
                   fng_chem = fng_chem + chemsputo(1,ii)*onesided_maxwellian(
     .                  tg(ix,ny+1,1), ng(ix,ny+1,ii), mg(ii), sy(ix,ny), tgmin*ev
     .             )
                enddo
                osmw = onesided_maxwellian(
     .                  tg(ix,ny+1,1), 1.0, mi(ifld), sy(ix,ny), tgmin*ev
     .          ) 
                yldot(iv1) = nurlxg*( fniy(ix,ny,ifld) 
     .                      - (1-albedoo(ix,1))*osmw
     .                      * harmave(ni(ix,ny,ifld), ni(ix,ny+1,ifld))
     .                      + fng_chem ) / (osmw * n0(ifld))
c...   Caution: the wall source models assume gas species 1 only is inertial
                if(matwallo(ix) .gt. 0) then
                  if (recycwot(ix,1) .gt. 0.) then
                    fniy_recy = recycwot(ix,1)*fac2sp*fniy(ix,ny,1)
                    if(isrefluxclip==1) 
     .                  fniy_recy=max(fniy_recy,0.)
                    yldot(iv1) = nurlxg*( fniy(ix,ny,ifld) 
     .                      + fniy_recy
     .                      + fngyo_use(ix,1)
     .                      + fngyso(ix,1)
     .                      - (1-albedoo(ix,1))*osmw
     .                      * harmave(ni(ix,ny,ifld), ni(ix,ny+1,ifld))
     .                      + fng_chem 
     .                  ) / (osmw * n0(ifld))
                  elseif (recycwot(ix,1) < -1) then
                    yldot(iv1) = nurlxg*(
     .                      ngbackg(1) - ni(ix,ny+1,ifld)
     .                  ) / n0(ifld)
                  elseif (recycwot(ix,1) .le. 0.) then  # treat recycwot as albedo
                    yldot(iv1) = nurlxg*( fniy(ix,ny,ifld) 
     .                      - (1+recycwot(ix,1))*ni(ix,ny+1,ifld)*osmw
     .                  ) / (osmw * n0(ifld))
                  endif 
                endif
                if(fngyso(ix,1)+fngyo_use(ix,1).ne.0. .and. matwallo(ix)==0.)
     .                        yldot(iv1) = nurlxg*( fniy(ix,ny,ifld) 
     .                              + fngyo_use(ix,1) 
     .                              + fngyso(ix,1)  
     .                          ) / (osmw * n0(ifld))



            else

























































              if(isnwconoix(ix,ifld) .eq. 0) then
                 yldot(iv1) = nurlxn * ( (1-ifluxni)*
     .                         (niy0(ix,ny,ifld) - niy1(ix,ny,ifld))
     .            + ifluxni*( fniy(ix,ny,ifld)/(sy(ix,ny)*vpnorm) 
     .           - 0.001*ni(ix,ny,ifld)*vy(ix,ny,ifld)/vpnorm ) )
     .                                                      /n0(ifld)
c...  the last term going as 0.001 is to prevent very small densities
c...  This next if is a recalc. if constant ni wanted - could be done better
              elseif (isnwconoix(ix,ifld) .eq. 1) then
                 yldot(iv1) =nurlxn*(nwallo(ix) - ni(ix,ny+1,ifld))/
     .                                                      n0(ifld)
              elseif (isnwconoix(ix,ifld) .eq. 2) then  #extrapolation
                 nbound = ni(ix,ny,ifld) + gyf(ix,ny-1)*
     .                      (ni(ix,ny,ifld)-ni(ix,ny-1,ifld))/gyf(ix,ny)
                 nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                         ni(ix,ny,ifld)-1) ) ) + 0.2*ni(ix,ny,ifld)
ccc                 nbound = max(nbound, 0.3*ni(ix,ny,ifld))
                 yldot(iv1) = nurlxn*(nbound - ni(ix,ny+1,ifld))/n0(ifld)
              elseif (isnwconoix(ix,ifld) .eq. 3) then   #spec. gradient 
                  yldot(iv1) = -nurlxn*( niy1(ix,ny,ifld) -
     .              niy0(ix,ny,ifld)*(2*gyf(ix,ny)*lyniix(2,ix,ifld)-1)/
     .                               (2*gyf(ix,ny)*lyniix(2,ix,ifld)+1) -
     .                         nwomin(ifld) ) / n0(ifld)

              endif            # endif for if (isnwconoix.eq.1 .and. ..
            endif              # endif for neut. ni; i.e., if (isupgon .eq. 1
          endif                # endif for isnionxy
  278   continue               # end loop over ix for ni BC








c...  Reset density at corners
          if (ixmnbcl*iymxbcl.eq.1) then
            do jx = 1, nxpt
              if(isnionxy(ixlb(jx),ny+1,ifld)==1) then
             yldot(idxn(ixlb(jx),ny+1,ifld)) = nurlxn * 
     .             ( harmave(ni(ixlb(jx),ny,ifld),ni(ixlb(jx)+1,ny+1,ifld))
     .                        - ni(ixlb(jx),ny+1,ifld) ) / n0(ifld)
              endif
            enddo
          endif
          if (ixmxbcl*iymxbcl.eq.1) then
            do jx = 1, nxpt
              if(isnionxy(ixrb(jx)+1,ny+1,ifld)==1) then
             yldot(idxn(ixrb(jx)+1,ny+1,ifld)) = nurlxn * 
     .               ( harmave(ni(ixrb(jx)+1,ny,ifld),ni(ixrb(jx),ny+1,ifld))
     .                         - ni(ixrb(jx)+1,ny+1,ifld) ) / n0(ifld)
              endif
            enddo
          endif
  280 continue             # end loop over ifld for ni BC

c...  Do the parallel velocity BC along iy = ny+1

      do ifld = 1, nusp
         do ix = i4+1-ixmnbcl, i8-1+ixmxbcl
            if (isuponxy(ix,ny+1,ifld)==1) then
               iv2 = idxu(ix,ny+1,ifld)






















































               if (isupwoix(ix,ifld)==1) then  #zero parallel momentum flux
                  yldot(iv2) = nurlxu * fmiy(ix,ny,ifld) / 
     .                                 (vpnorm*sy(ix,ny)*fnorm(ifld))
               elseif (isupwoix(ix,ifld)==2) then  #  d(up)/dy = 0
                  yldot(iv2) = nurlxu * nm(ix,ny,ifld) / fnorm(ifld) *
     .                           (up(ix,ny,ifld) - up(ix,ny+1,ifld))
               elseif (isupwoix(ix,ifld)==3) then  # d(up)/dy = up/lyup
                  yldot(iv2) = -nurlxu * nm(ix,ny,ifld) / fnorm(ifld) *
     .                             ( up(ix,ny+1,ifld) - up(ix,ny,ifld)*
     .                                (2*gyf(ix,ny)*lyup(2)-1)/
     .                                      (2*gyf(ix,ny)*lyup(2)+1) )
               else
                  yldot(iv2) = nurlxu * nm(ix,ny,ifld) / fnorm(ifld) *
     .                                      (0. - up(ix,ny+1,ifld))
               endif
            endif 
         enddo
      enddo

c ... Do boundary conditions for impurities along iy = ny+1.
c     Force corner cells (0,ny+1) and (nx+1,ny+1) to relax to an average
c     of the adjacent cells.
      do ifld = 1, nzspt
        iimp = nhsp + ifld
        if (isimpon .ge. 3 .and. isimpon .le. 7) then
          do ix = i2, i5  # i2 and i5 limits omit the ix=0 and nx+1 corners
            if (isnionxy(ix,ny+1,iimp)==1) then
               iv = idxn(ix,ny+1,iimp)  
               






































               if (isnwconoix(ix,iimp) .eq. 0) then       #fix flux
                  yldot(iv) = nurlxn *
     .                        (fniy(ix,ny,iimp) + fnzyso(ix,ifld)) /
     .                            (sy(ix,ny) * n0(iimp) * vpnorm)
               elseif (isnwconoix(ix,iimp) .eq. 1) then   #set to nwallo
                  yldot(iv) =nurlxn*(nwallo(ix) - ni(ix,ny+1,iimp))/
     .                                                      n0(iimp)
               elseif (isnwconoix(ix,iimp) .eq. 2) then   #extrapolation
                  nbound = ni(ix,ny,iimp) + gyf(ix,ny-1)*
     .                      (ni(ix,ny,iimp)-ni(ix,ny-1,iimp))/gyf(ix,ny)
                  nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                         ni(ix,ny,iimp)-1) ) ) + 0.2*ni(ix,ny,iimp)
                  yldot(iv) = nurlxn*(nbound - ni(ix,ny+1,iimp))/n0(iimp)
               elseif (isnwconoix(ix,iimp) .eq. 3) then   #spec. gradient 
                  yldot(iv) = -nurlxn*( niy1(ix,ny,iimp) -
     .              niy0(ix,ny,iimp)*(2*gyf(ix,ny)*lyniix(2,ix,iimp)-1)/
     .                               (2*gyf(ix,ny)*lyniix(2,ix,iimp)+1) -
     .                         nwomin(iimp) ) / n0(iimp)
               endif
            endif       # end if for isnionxy
          


          enddo
c...   set the corner values here
          if (ixmnbcl*iymxbcl.eq.1) then
             do jx = 1, nxpt
               if (isnionxy(ixlb(jx),ny+1,iimp)==1) then
                  iv = idxn(ixlb(jx),ny+1,iimp)
                  yldot(iv) = nurlxn *
     .              ( harmave(ni(ixlb(jx),ny,iimp),ni(ixlb(jx)+1,ny+1,iimp))
     .                           - ni(ixlb(jx),ny+1,iimp) ) / n0(iimp)
               endif
             enddo
          endif
          if (ixmxbcl*iymxbcl.eq.1) then
             do jx = 1, nxpt
               if (isnionxy(ixrb(jx)+1,ny+1,iimp)==1) then
                  iv = idxn(ixrb(jx)+1,ny+1,iimp)
                  yldot(iv) = nurlxn *
     .                 ( harmave(ni(ixrb(jx)+1,ny,iimp),ni(ixrb(jx),ny+1,iimp))
     .                           - ni(ixrb(jx)+1,ny+1,iimp) ) / n0(iimp)
               endif
             enddo
          endif
        endif     # end for if (isimpon .ge. 3 ...
      enddo        # end of impurities along iy = ny+1










































ccc  Now do Te, Ti
ccc  - - - - - - - - - - - - - -
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl  # loops over te,ti,ng,phi
         ix1 = ixm1(ix,ny+1)
         if(isteonxy(ix,ny+1) .eq. 1) then
           iv1 = idxte(ix,ny+1)
          





















           if (istewcix(ix) .eq. 0) then       # fix flux to zero
             yldot(iv1) = nurlxe*(feey(ix,ny)/(n0(1)*vpnorm*sy(ix,ny)))
     .                                                      / (temp0*ev)
           elseif (istewcix(ix) .eq. 1) then   # set to tewallo
             yldot(iv1) = nurlxe*(tewallo(ix)*ev-te(ix,ny+1))/(temp0*ev)
           elseif (istewcix(ix) .eq. 2) then   # extrapolation
             tbound = te(ix,ny) + gyf(ix,ny-1)*(te(ix,ny)-
     .                           te(ix,ny-1))/gyf(ix,ny)
             tbound = max(tbound, tbmin*ev)
             yldot(iv1) = nurlxe*(tbound - te(ix,ny+1)) / (temp0*ev) 
           elseif (istewcix(ix) .eq. 3) then   #spec. gradient 
             yldot(iv1) = nurlxe*( (te(ix,ny) - te(ix,ny+1)) -
     .                         0.5*(te(ix,ny) + te(ix,ny+1))/
     .                           (gyf(ix,ny)*lyteix(2,ix)) )/(temp0*ev)
           elseif (istewcix(ix) .eq. 4) then   #heat flux ~bcee*te*elec_flux
             yldot(iv1) = nurlxe*( feey(ix,ny) - bceew*ne(ix,ny)*
     .                             vey(ix,ny)*sy(ix,ny)*te(ix,ny+1) )/
     .                                      (vpnorm*ennorm*sy(ix,ny))
            endif
         endif

         if(istionxy(ix,ny+1) .eq. 1) then
           iv2 = idxti(ix,ny+1)
           





























           if (istiwcix(ix) .eq. 0) then       # fix flux to zero
             yldot(iv2) = nurlxi*(feiy(ix,ny)/(n0(1)*vpnorm*sy(ix,ny)))
     .                                                      / (temp0*ev)
           elseif (istiwcix(ix) .eq. 1) then   # set to tiwallo
             yldot(iv2) = nurlxi*(tiwallo(ix)*ev-ti(ix,ny+1))/(temp0*ev)
           elseif (istiwcix(ix) .eq. 2) then   # extrapolation
             tbound = ti(ix,ny) + gyf(ix,ny-1)*(ti(ix,ny)-
     .                           ti(ix,ny-1))/gyf(ix,ny)
             tbound = max(tbound, tbmin*ev)
             yldot(iv2) = nurlxi*(tbound - ti(ix,ny+1)) / (temp0*ev)
           elseif (istiwcix(ix) .eq. 3) then   #spec. gradient 
             yldot(iv2) = nurlxi*( (ti(ix,ny) - ti(ix,ny+1)) -
     .                         0.5*(ti(ix,ny) + ti(ix,ny+1))/
     .                            (gyf(ix,ny)*lytiix(2,ix)) )/(temp0*ev)
           elseif (istiwcix(ix) .eq. 4) then   #heat flux ~bcei*ti*elec_flux
             t0 = max(tg(ix,ny+1,1), temin*ev)
             fngyw=0.25*sqrt(8*t0/(pi*mg(1)))*ng(ix,ny,1)*sy(ix,ny)
             yldot(iv2) = nurlxi*( feiy(ix,ny) - bceiw*ne(ix,ny)*
     .                             vey(ix,ny)*sy(ix,ny)*ti(ix,ny+1) -
     .              cftiexclg*( bcenw*fniy(ix,ny,iigsp)*tg(ix,ny+1,1) +
     .                              cgengw*2.*tg(ix,ny+1,1)*fngyw ) ) /
     .                                       (vpnorm*ennorm*sy(ix,ny))
           endif
         endif
       enddo   # ix-loop for Te,Ti



ccc  Now do the diffusive neutral density (ng) and Tg equations
ccc  - - - - - - - - - - - - - -          
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl  # ix-loop for ng & Tg
        nzsp_rt = nhsp
        do igsp = 1, ngsp
          jz = max(igsp - nhgsp, 1)     #identify impurity index
      if (jz > 1) nzsp_rt = nzsp_rt + nzsp(jz-1) #prev index for fniy
          if(isngonxy(ix,ny+1,igsp) .eq. 1) then #skip do-loop if isngon=0
            if (iscplo(ix) .eq. 1) call wsmodo(igsp)
            iv = idxg(ix,ny+1,igsp)
            if (igsp .gt. nhgsp) then
              zflux = 0.
              do iimp = 1, nzsp(jz)
                    zflux = zflux + max(fniy(ix,ny,nzsp_rt+iimp), 0.)
              enddo
            endif




































c ... prepare chemical sputtering info
            fng_chem = 0.
            do igsp2 = 1, 1+ishymol  #include hydrogen neut sputt
                flx_incid = onesided_maxwellian(
     .                      tg(ix,ny+1,igsp2), ng(ix,ny+1,igsp2),
     .                      mg(igsp2), 1.0, tgmin*ev
     .          )
                if (isch_sput(igsp).ne.0) then
                    call sputchem(
     .                      isch_sput(igsp), 
     .                      max(tg(ix,ny+1,igsp2)/ev, tgmin), 
     .                      tvwallo(ix),
     .                      flx_incid, 
     .                      yld_carbo(ix)
     .              )
                    fng_chem = fng_chem + fchemygwo(igsp2)*yld_carbo(ix)*
     .                                            flx_incid*sy(ix,ny)
                else
                    fng_chem = fng_chem + chemsputo(igsp,igsp2)*
     .                                            flx_incid*sy(ix,ny) 
                endif
            enddo
c ... add ion sputtering to gas BC
            sputflxw(ix,igsp) = 0.
            if (isi_sputw(igsp).ge.1 .and. igsp>1) then  # incl phys sput from ions
              do ifld = ipsputt_s, ipsputt_e 
                eng_sput = ( 2*ti(ix,ny+1) + zi(ifld)*
     .                                          3*te(ix,ny+1) )/ev
                if(zi(ifld) > 0) sputflxw(ix,igsp) = sputflxw(ix,igsp) +
     .                            max(fniy(ix,ny,ifld),0.)*fphysyiwo(igsp)*
     .                                    yld96(matp,matt,eng_sput)
              enddo
              if (isi_sputw(igsp) .ge. 2) then  # add chem sput from ions
                do ifld = 1,1  #(ipsputt_s, ipsputt_e) place holder for imp/imp sputt
                  eng_sput = ( 2*ti(ix,ny+1) + zi(ifld)*
     .                                           3*te(ix,ny+1) )/ev
                  flx_incid = max(fniy(ix,ny,ifld)/sy(ix,ny), 0.)
                  call sputchem (isch_sput(igsp), eng_sput, tvwallo(ix),
     .                           flx_incid, yld_chm)
                  sputflxw(ix,igsp) = sputflxw(ix,igsp) + fchemyiwo(igsp)*
     .                                   max(fniy(ix,ny,ifld),0.)*yld_chm
                enddo
              endif
            endif
            osmw = onesided_maxwellian(
     .              cdifg(igsp)*tg(ix,ny+1,igsp), 1.0, mg(igsp),
     .              sy(ix,ny), tgmin*ev
     .      )
            yldot(iv) = nurlxg*( fngy(ix,ny,igsp) 
     .                  - (1-albedoo(ix,igsp))*osmw
     .                      * harmave(ng(ix,ny,igsp), ng(ix,ny+1,igsp))
     .                  + fng_chem
     .                  + sputflxw(ix,igsp) 
     .              ) / (osmw *n0g(igsp))
            
            if(matwallo(ix) .gt. 0) then
              if (recycwot(ix,igsp) .gt. 0.) then
ccc
ccc   MER 01 Apr 2002: need to correct fniy below when drifts are included
ccc
                fniy_recy = fac2sp*fniy(ix,ny,1)
                if(isrefluxclip==1) 
     .                  fniy_recy=max(fniy_recy,0.)
                if (igsp .gt. nhgsp) 
     .                  fniy_recy = zflux
                if (ishymol.eq.1 .and. igsp.eq.2) then # 2 atoms per molecule
                  if (isupgon(1) .eq. 1) then
                        fniy_recy = 0.5*( fniy(ix,ny,1) + fniy(ix,ny,2) )
                  else
                        fniy_recy = 0.5*( fniy(ix,ny,1) + fngy(ix,ny,1) )
                  endif
                  if(isrefluxclip==1)
     .                   fniy_recy=max(fniy_recy,0.)
                endif
                yldot(iv) = nurlxg*( fngy(ix,ny,igsp) 
     .                  + fniy_recy*recycwot(ix,igsp) 
     .                  + fngyso(ix,igsp) 
     .                  + fngyo_use(ix,igsp) 
     .                  - (1-albedoo(ix,igsp))*osmw
     .                      * harmave(ng(ix,ny,igsp), ng(ix,ny+1,igsp))
     .                  + fng_chem 
     .                  + sputflxw(ix,igsp)
     .              ) / (osmw*n0g(igsp))
              

              elseif (recycwot(ix,igsp) < -1) then
                yldot(iv)=nurlxg*(ngbackg(igsp)-ng(ix,ny+1,igsp))/
     .                                                        n0g(igsp)
             elseif (recycwot(ix,igsp) .le. 0.) then # treat recycw as albedo
                nharmave = 2.*(ng(ix,ny,igsp)*ng(ix,ny+1,igsp)) /
     .                        (ng(ix,ny,igsp)+ng(ix,ny+1,igsp))              
                yldot(iv) = nurlxg*( fngy(ix,ny,igsp) -
     .                  - (1+recycwot(ix,igsp))*osmw
     .                      * harmave(ng(ix,ny,igsp), ng(ix,ny+1,igsp))
     .              ) / (osmw * n0g(igsp))
              endif 
            endif



            if( (fngyso(ix,igsp)+fngyo_use(ix,igsp).ne.0.)
     .          .and. (matwallo(ix).eq.0) 
     .      ) then
                    yldot(iv) = nurlxg*(fngy(ix,ny,igsp)
     .                      + fngyo_use(ix,igsp) 
     .                      + fngyso(ix,igsp) 
     .                  ) / (osmw * n0g(igsp))
            endif
          endif     # end if-test on isngonxy
 







c... BC for neutral gas temperature/energy at iy=ny+1
          if (istgonxy(ix,ny+1,igsp) == 1) then
            iv = idxtg(ix,ny+1,igsp)
            



















            if (istgwcix(ix,igsp) == 0) then    # fixed Tg
              yldot(iv) = nurlxg*(tgwall(igsp)*ev-tg(ix,ny+1,igsp))/
     .                                                    (temp0*ev)
            elseif (istgwcix(ix,igsp) == 1)    # extrapolation
              tbound = tg(ix,ny,igsp) + gyf(ix,ny)*
     .                     (tg(ix,ny,igsp)-tg(ix,ny-1,igsp))/gyf(ix,ny)
              tbound = max(tbound, 0.25*tbmin*ev)  #tbmin=.1 eV
              yldot(iv) = nurlxi *(tbound - tg(ix,ny+1,igsp))/(temp0*ev)
            elseif (istgwcix(ix,igsp) == 2)    # specified gradient
              yldot(iv) = nurlxi*( (tg(ix,ny,igsp) - tg(ix,ny+1,igsp)) -
     .                         0.5*(tg(ix,ny,igsp) + tg(ix,ny+1,igsp))/
     .                         (gyf(ix,ny)*lytg(2,igsp)) )/(temp0*ev)
            elseif (istgwcix(ix,igsp) == 3)  #Maxwell thermal flux to wall
              osmw = onesided_maxwellian(
     .              cdifg(igsp)*tg(ix,ny,igsp), ng(ix,ny,igsp),
     .              mg(igsp), sy(ix,ny), tgmin*ev
     .        )
              yldot(iv) =  nurlxg*( 
     .              fegy(ix,ny,igsp) - 2*cgengmw*osmw*max(tg(ix,ny,igsp),tgmin*ev)
     .          ) / (sy(ix,ny)*vpnorm*ennorm)
            elseif (istgwcix(ix,igsp) == 4) 
              fng_chem = 0.
              osmw = onesided_maxwellian(
     .              tg(ix,ny+1,igsp), 
     .              harmave(ng(ix,ny,igsp),ng(ix,ny+1,igsp)),
     .              mg(igsp), sy(ix,ny), tgmin*ev
     .        )
              yldot(iv) = nurlxg*(fegy(ix,ny,igsp) 
     .                  - cfalbedo*max(tg(ix,ny+1,igsp),tgmin*ev)*(1-albedoo(ix,1))*osmw
     .                  + 2.*fng_chem*max(tg(ix,ny+1,igsp), tgmin*ev)
     .              ) / (vpnorm*ennorm*sy(ix,ny))
              fniy_recy = 0.
              if (matwallo(ix) .gt. 0) then
                if (recycwot(ix,igsp) .gt. 0.) then
                  fniy_recy = recycwot(ix,igsp)*fac2sp*fniy(ix,ny,1)
                  if (isrefluxclip==1) fniy_recy=max(fniy_recy,0.)
                  yldot(iv)=nurlxg*(fegy(ix,ny,igsp) - cfalbedo*fng_alb
     .                      *max(tg(ix,ny+1,igsp), tgmin*ev)
     .                      + 2.*fng_chem*max(tg(ix,ny+1,igsp),tgmin*ev)
     .                                             + fniy_recy*(1.-cfdiss)
     .                                              *cfalbedo*recycwe
     .                                              *ti(ix,ny) )
     .                                     /(vpnorm*ennorm*sy(ix,ny))
                endif
              endif
              elseif (istgwc(igsp) == 5) then    # set tg=ti*cftgtiwc
              yldot(iv) = nurlxg*(ti(ix,ny+1)*cftgtiwc(igsp) -
     .                                   tg(ix,ny+1,igsp))/(temp0*ev)
            elseif (istgwcix(ix,igsp) > 5)
               call xerrab("***Input error: invalid istgwc ***")
            endif
          endif

        enddo  # igsp loop over gas species
      enddo  # ix-loop for ng and Tg



ccc  Now do the potential
ccc  - - - - - - - - - - - -
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl  # ix-loop for phi
         if(isphionxy(ix,ny+1) .eq. 1) then
            iv3 = idxphi(ix,ny+1)
	        iv2 = idxphi(ix,ny)
            if(iphibcwoix(ix) == 0) then
              yldot(iv3) = nurlxp*(phi(ix,ny) - phi(ix,ny+1))/temp0
            elseif(iphibcwoix(ix) == 1) then
              yldot(iv3) = nurlxp*
     .               (phintewo*te(ix,ny+1)/ev - phi(ix,ny+1))/temp0
	    elseif(iphibcwoix(ix) == 3) then
              yldot(iv3) = nurlxp*( (phi(ix,ny) - phi(ix,ny+1)) -
     .                         0.5*(phi(ix,ny) + phi(ix,ny+1))/
     .                            (gyf(ix,ny)*lyphiix(2,ix)) )/temp0
	    elseif(iphibcwoix(ix) == 4) then  #fixed prof = phiwo
              yldot(iv3) = nurlxp*(phiwo(ix) - phi(ix,ny+1))/temp0
ccc              yldot(iv2) = nurlxp*(phiwo(ix) - phi(ix,ny))/temp0
            endif
         endif
      enddo  # ix-loop for phi





















c...  Now do the corners: the ion density corners are done elsewhere
       if ((xcnearlb.or.openbox) .and. ixmnbcl*iymxbcl.eq.1) then
         do ifld = 1, nusp
           do jx = 1, nxpt
             if(isuponxy(ixlb(jx),ny+1,ifld)==1) then
               yldot(idxu(ixlb(jx),ny+1,ifld)) = - nurlxu*
     .                      ( up(ixlb(jx),ny+1,ifld) -
     .             0.5*(up(ixlb(jx),ny,ifld)+up(ixlb(jx)+1,ny+1,ifld)) )
     .                                         /vpnorm
             endif
           enddo
         enddo
         do jx = 1, nxpt
            if(isteonxy(ixlb(jx),ny+1)==1) then
               yldot(idxte(ixlb(jx),ny+1)) = nurlxe*
     .              ( 0.5*(te(ixlb(jx)+1,ny+1)+te(ixlb(jx),ny))
     .                               - te(ixlb(jx),ny+1) ) / (temp0*ev)
            endif
         enddo
         do jx = 1, nxpt
            if(istionxy(ixlb(jx),ny+1)==1) then
               yldot(idxti(ixlb(jx),ny+1)) = nurlxi*
     .              ( 0.5*(ti(ixlb(jx)+1,ny+1)+ti(ixlb(jx),ny))
     .                               - ti(ixlb(jx),ny+1) ) / (temp0*ev)
            endif
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
             if(isngonxy(ixlb(jx),ny+1,igsp)==1) then
               yldot(idxg(ixlb(jx),ny+1,igsp)) = nurlxg*
     .            ( ng(ixlb(jx)+1,ny+1,igsp) - ng(ixlb(jx),ny+1,igsp) )
     .                                                      / n0g(igsp)
             endif
             if(istgonxy(ixlb(jx),ny+1,igsp)==1) then 
               yldot(idxtg(ixlb(jx),ny+1,igsp)) = nurlxg*
     .            ( 0.5*(tg(ixlb(jx)+1,ny+1,igsp)+tg(ixlb(jx),ny,igsp))
     .                                       - tg(ixlb(jx),ny+1,igsp) )
     .                                                      / (temp0*ev)
             endif
           enddo
         enddo
       endif

       if ((xcnearrb.or.openbox) .and. ixmxbcl*iymxbcl.eq.1) then
         do ifld = 1, nusp
           do jx = 1, nxpt
             if(isuponxy(ixrb(jx),ny+1,ifld)==1) then
               yldot(idxu(ixrb(jx),ny+1,ifld)) = - nurlxu*
     .                      ( up(ixrb(jx),ny+1,ifld) -
     .             0.5*(up(ixrb(jx)-1,ny+1,ifld)+up(ixrb(jx),ny,ifld)) )
     .                                         /vpnorm
             endif
             if(isuponxy(ixrb(jx)+1,ny+1,ifld)==1) then
               yldot(idxu(ixrb(jx)+1,ny+1,ifld)) = - nurlxu*
     .                   ( up(ixrb(jx)+1,ny+1,ifld)-up(ixrb(jx),ny+1,ifld) )
     .                                            /vpnorm
             endif
           enddo
         enddo
         do jx = 1, nxpt
            if (isteonxy(ixrb(jx)+1,ny+1)==1) then
               yldot(idxte(ixrb(jx)+1,ny+1)) = nurlxe*
     .              ( 0.5*(te(ixrb(jx),ny+1)+te(ixrb(jx)+1,ny))
     .                               - te(ixrb(jx)+1,ny+1) ) / (temp0*ev)
            endif
         enddo
         do jx = 1, nxpt
            if (istionxy(ixrb(jx)+1,ny+1)==1) then
               yldot(idxti(ixrb(jx)+1,ny+1)) = nurlxi*
     .              ( 0.5*(ti(ixrb(jx),ny+1)+ti(ixrb(jx)+1,ny))
     .                               - ti(ixrb(jx)+1,ny+1) ) / (temp0*ev)
            endif
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
             if(isngonxy(ixrb(jx)+1,ny+1,igsp)==1) then
               yldot(idxg(ixrb(jx)+1,ny+1,igsp)) = nurlxg*
     .            ( ng(ixrb(jx),ny+1,igsp) - ng(ixrb(jx)+1,ny+1,igsp) )
     .                                                      / n0g(igsp)
             endif
             if(istgonxy(ixrb(jx)+1,ny+1,igsp)==1) then 
               yldot(idxtg(ixrb(jx)+1,ny+1,igsp)) = nurlxg*
     .            ( 0.5*(tg(ixrb(jx),ny+1,igsp)+tg(ixrb(jx)+1,ny,igsp)) 
     .                               - tg(ixrb(jx)+1,ny+1,igsp) )
     .                                                      / (temp0*ev)
             endif
           enddo
         enddo
       endif  # end of if-test on xcnearrb
      endif   # end of iy = ny+1 boundary conditions except if isnewpot=1

c...  If isnewpot=1, phi boundary involves two eqns at iy=ny and ny+1
CCC   NOTE: special coding for testing alternate phi B.C. (requires 
CCC   isnewpot*isphion=1000, so one can generally ignore this if section) 
      if (isnewpot*isphion.eq.1000 .and. j7.ge.ny-2) then 
        do ix = i4+1-ixmnbcl, i8-1+ixmxbcl
          if(isphionxy(ix,ny+1)*isphionxy(ix,ny)==1) then
            iv  = idxphi(ix,ny+1)
            iv1 = idxphi(ix,ny)
            ix3 = ixm1(ix,ny)
            ix4 = ixm1(ix,ny-1)
            yldot(iv) = nurlxp*( (phi(ix,ny) - phi(ix,ny+1)) -
     .                         0.5*(phi(ix,ny) + phi(ix,ny+1))/
     .                             (gyf(ix,ny)*lyphiix(2,ix)) )/temp0
            yldot(iv1) = -nurlxp*( phi(ix,ny) - phi(ix,ny+1) +
     .                           gyf(ix,ny-1)*(phi(ix,ny)-phi(ix,ny-1))/
     .                                                      gyf(ix,ny) )
          endif    # end of if-test ons isphionxy
        enddo 
      endif

        

        RETURN
        END SUBROUTINE north_boundary