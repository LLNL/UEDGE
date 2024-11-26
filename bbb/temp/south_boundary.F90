        SUBROUTINE south_boundary(neq, yl, yldot)
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

          real yld96

c ====================================================================
c ======================== The iy=0 boundary =========================
c ====================================================================

c...  do the iy = 0 boundary
c...  if extrapolation b.c. on p.f. region, isextrpf=1, otherwise isextrpf=0
      ixc1 = max(0, ixpt1(1)+1)       # 1st  core cell;used for core flux BC
      if ((geometry(1:9)=="snowflake" .and. geometry(10:11).ne."15")
     &                     .or. geometry=="dnXtarget") then
        ix_fl_bc = min(ixpt2(1), nx)
      else
        ix_fl_bc = min(ixpt2(nxpt), nx) # last core cell;use for flux BC
      endif 

c  Note: j3 is local range index for iy passed from pandf in oderhs.m
      if (j3 .le. isextrnpf .or. j3 .le. isextrtpf .or. 
     .                                         j3 .le. isextrngc) then
              # this is a very large if-loop; check for 'j3 index' to find end

      do 274 ifld = 1, nisp
         do 272 ix = i4+1-ixmnbcl, i8-1+ixmxbcl # ix-loop over ion dens
           if (isnionxy(ix,0,ifld)==1) then
             iv1 = idxn(ix,0,ifld)

c ...  First check if this is a neutral species (zi=0) & set BC for it
c ---  Different boundary conditions if neutral momentum equation use;
c ---  typically hydrogen only, so DIVIMP chem sputt not used here.
c ---  Below, isupon(1) & ngbackg(1) used, so implies hydrogen
             if (isupgon(1) .eq. 1 .and. zi(ifld) .eq. 0.0) then
               if (isixcore(ix)==1) then   # ix is part of the core boundary:
                 if (isngcore(1) .eq. 0) then
                   osmw = onesided_maxwellian(
     .                  tg(ix,0,1), harmave(ni(ix,0,ifld),ni(ix,1,ifld)), 
     .                  mi(ifld), sx(ix,0), tgmin*ev)
                   yldot(iv1) = -nurlxg*(fniy(ix,0,ifld) 
     .                  + (1-albedoc(1))*osmw) / (vpnorm*sy(ix,0)*n0(ifld))
                 elseif (isngcore(1) .eq. 1) then
                   yldot(iv1) = nurlxn*(ngcore(1)-ni(ix,0,ifld))/
     .                                                     n0(ifld)
                 elseif (isngcore(1) .eq. 2) then
                   call xerrab("*** isngcore=2 option unvailable  ***")
                   lengg = sqrt(ti(ix,0)/
     .                              (mi(1)*nuix(ix,0,1)*nuiz(ix,0,1)))
                   yldot(iv1) = -nurlxn*( ni(ix,0,ifld) -
     .                           max( ni(ix,1,ifld) -
     .                        0.5*(ni(ix,1,ifld) + ni(ix,0,ifld))/
     .                        (gyf(ix,0)*lengg), 0.1*ngbackg(1) ) )
     .                                                       /n0(ifld)
                 elseif (isngcore(1) .eq. 3) then
                   nbound =  ng(ix,1,1) - gyf(ix,1)*
     .                          (ng(ix,2,1)-ng(ix,1,1))/gyf(ix,0)
                   nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                           ng(ix,1,1)-1) ) ) + 0.2*ng(ix,1,1)
                   yldot(iv1) = nurlxn *(nbound - ng(ix,0,1))/n0(ifld)
                 else  # old condition before 2/21/99
                   yldot(iv1) = nurlxn*(ni(ix,1,ifld)-ni(ix,0,ifld))/
     .                                                    n0(ifld)
                 endif
               else   # ix is not part of the core boundary
                  if (iscpli(ix) .eq. 1) call wsmodi(1)
                  fng_chem = 0.
                  do ii = 1, ngsp  #chem sputt of hydrogen - strange = 0
                     fng_chem= fng_chem + chemsputi(1,ii)*onesided_maxwellian(
     .                      tg(ix,0,1), ng(ix,1,ii), mi(1), sy(ix,0), tgmin*ev)
                  enddo
                  osmw = onesided_maxwellian(
     .                      tg(ix,0,1), 1.0, mi(1), sy(ix,0), tgmin*ev
     .            )
                  yldot(iv1) = -nurlxg*( fniy(ix,0,ifld) 
     .                      + (1-albedoi(ix,1))*osmw
     .                      * harmave(ni(ix,0,1), ni(ix,1,ifld)) -fng_chem
     .                  )/ (osmw*n0(ifld))
c...   Caution: the wall source models assume gas species 1 only is inertial
                  if(matwalli(ix) .gt. 0) then
                    if (recycwit(ix,1,1) .gt. 0.) then  
                      fniy_recy = recycwit(ix,1,1)*fac2sp*fniy(ix,0,1)
                      if (isrefluxclip==1) fniy_recy=min(fniy_recy,0.)
                        osmw = onesided_maxwellian(
     .                          tg(ix,0,1), 1.0, mi(1), sy(ix,0), tgmin*ev
     .                  )
                        yldot(iv1)=-nurlxg*( fniy(ix,0,ifld) + fniy_recy -
     .                      fngyi_use(ix,1) - fngysi(ix,1) 
     .                      + (1-albedoi(ix,1)) * osmw 
     .                      * harmave(ni(ix,0,ifld), ni(ix,1,ifld)) - fng_chem  
     .                      ) / (osmw*n0(ifld))
                    elseif (recycwit(ix,1,1) < -1) then
                      yldot(iv1)=nurlxg*(ngbackg(1)-ni(ix,0,ifld))/n0(ifld)
                    elseif (recycwit(ix,1,1) .le. 0.) then  # treat recycwit as albedo
                      osmw = onesided_maxwellian( 
     .                          tg(ix,0,1), 1.0, mi(1), sy(ix,0), tgmin*ev
     .                )
                      yldot(iv1) = -nurlxg*( fniy(ix,0,ifld) +
     .                   (1+recycwit(ix,1,1))*osmw
     .                   *harmave(ni(ix,0,ifld), ni(ix,1,ifld)) 
     .                   ) / (osmw * n0(ifld))
                    endif                  
                  endif
                  if(fngysi(ix,1)+fngyi_use(ix,1) .ne. 0. 
     .                                           .and. matwalli(ix)==0.)
     .               yldot(iv1)= -nurlxg*(fniy(ix,0,ifld) - fngysi(ix,1) -
     .                            fngyi_use(ix,1) ) / (osmw * n0(ifld))
               endif   # end if-test for core and p.f. boundaries



             else  # ifld is NOT inertial neutrals but still hydrogen ion
c
               if (isixcore(ix)==1) then   # ix is part of the core boundary:
                  if (isnicore(ifld) .eq. 1) then # set uniform, fixed density
                     yldot(iv1) = nurlxn*(ncore(ifld)-ni(ix,0,ifld))/
     .                                                          n0(ifld)
                  elseif (isnicore(ifld) .eq. 0) then # set flux to curcore
                     yldot(iv1) = - nurlxn* ( qe*
     .                       (fniy(ix,0,ifld)-fniycbo(ix,ifld))/sy(ix,0)
     .                            - curcore(ifld)*gyf(ix,0)/sygytotc ) /
     .                                              (qe*vpnorm*n0(ifld))
                  elseif (isnicore(ifld) .eq. 2) then # set flux & ni over range
                     yldot(iv1) = - nurlxn* ( qe*
     .                       (fniy(ix,0,ifld)-fniycbo(ix,ifld))/sy(ix,0)
     .                            - curcore(ifld)*gyf(ix,0)/sygytotc ) /
     .                                              (qe*vpnorm*n0(ifld))
                     if (ix.ge.ixfixnc .and. ix.le.ixfixnc+incixc) then
                        yldot(iv1) = nurlxn*(ncore(ifld)-ni(ix,0,ifld))/
     .                                                           n0(ifld)
                     endif                       
           elseif (isnicore(ifld) .eq. 3) then # const ni; flux set to
                                                      # curcore-recycc*fngy
                     yldot(iv1) = -nurlxn*( ni(ix,0,ifld) -
     .                                ni(ixp1(ix,0),0,ifld) ) / n0(ifld)




                     if (ix.eq.ix_fl_bc) then
                        ii=max(0,ixpt1(1)+1)
                        fniytotc=fniy(ii,0,ifld)-fniycbo(ii,ifld)
                        fngytotc = fngy(ii,0,1)  # needs generalization
                        do
                          ii=ixp1(ii,0)
                          fniytotc = fniytotc + fniy(ii,0,ifld) -
     .                                          fniycbo(ii,ifld)
                          fngytotc = fngytotc + fngy(ii,0,1) #generalize
                          if (ii==ix_fl_bc) break
                        enddo
            if (ifld > 1) fngytotc = 0.  # crude fix for imp.
                        yldot(iv1)= -nurlxn*(qe*fniytotc +
     .                                       qe*recycc(1)*fngytotc -
     .                               curcore(ifld))/(qe*vpnorm*n0(ifld))
                     endif
                  elseif (isnicore(ifld)==5) then # set dni/dy=-ni/lynicore
                                                  # at midpl; the pol const ni
                     if (ix .ne. ixmp) then  #  flux func for ni
                       yldot(iv1) = -nurlxn*( ni(ix,0,ifld) -
     .                                 ni(ixp1(ix,0),0,ifld) ) / n0(ifld)
                     else  # close by setting midp gradient 1/lynicore
                       yldot(iv1)=-nurlxn*( ni(ix,0,ifld) - ni(ix,1,ifld)*
     .                                 (2*gyf(ixmp,0)*lynicore(ifld)+1)/
     .                                 (2*gyf(ixmp,0)*lynicore(ifld)-1)- 
     .                                     ncoremin(ifld) ) / n0(ifld)
                     endif
                  else
                    call xerrab ('** isnicore value not valid option **')
                  endif
               elseif (isnwconiix(ix,ifld) .eq. 0) then
               # ix not on core boundary; set zero gradient (or flux)
                  yldot(iv1) = nurlxn * ( (1-ifluxni)*
     .                      (niy1(ix,0,ifld) - niy0(ix,0,ifld))
     .                 - ifluxni*( fniy(ix,0,ifld)/(sy(ix,0)*vpnorm)    
     .                 - 0.001*ni(ix,1,ifld)*vy(ix,0,ifld)/vpnorm 
     .                 ))/n0(ifld)
c...  the last term going as 0.001 is to prevent very small densities
               elseif (isnwconiix(ix,ifld) .eq. 1) then
               # ix is not part of the core boundary; set fixed density
                  yldot(iv1) = nurlxn*(nwalli(ix) - ni(ix,0,ifld))/
     .                                                        n0(ifld)
               elseif (isnwconiix(ix,ifld) .eq. 2) then
               # ix is not part of the core boundary; set extrapolation bc
                  nbound =  ni(ix,1,ifld) - gyf(ix,1)*
     .                         (ni(ix,2,ifld)-ni(ix,1,ifld))/gyf(ix,0)
                  nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                          ni(ix,1,ifld)-1) ) ) + 0.2*ni(ix,1,ifld)
                  yldot(iv1) = nurlxn *(nbound - ni(ix,0,ifld))/n0(ifld)
               elseif (isnwconiix(ix,ifld) .eq. 3) then  #spec. gradient
                    yldot(iv1) =-nurlxn*( niy0(ix,0,ifld) -
     .                niy1(ix,0,ifld)*(2*gyf(ix,0)*lyniix(1,ix,ifld)-1)/
     .                                (2*gyf(ix,0)*lyniix(1,ix,ifld)+1) -
     .                                  nwimin(ifld) ) / n0(ifld)
               endif   # end if-test for core and p.f. boundaries

             endif           # endif for if (isupgon.eq.1 ...
           endif             # end for if-test on isnionxy
  272    continue  # ix-loop for ion dens







c...  Reset density at corners
         do jx = 1, nxpt
           if(isfixlb(jx).ne.2 .and. ixmnbcl*iymnbcl.eq.1 .and.
     .                       isnionxy(ixlb(jx),0,ifld)==1) then
         yldot(idxn(ixlb(jx),0,ifld)) = nurlxn * 
     .            ( harmave(ni(ixlb(jx),1,ifld),ni(ixlb(jx)+1,0,ifld))
     .                        - ni(ixlb(jx),0,ifld) ) / n0(ifld)
           endif
           if (isfixrb(jx).ne.2 .and. ixmxbcl*iymnbcl.eq.1 .and.
     .                      isnionxy(ixrb(jx)+1,0,ifld)==1) then
          yldot(idxn(ixrb(jx)+1,0,ifld)) = nurlxn * 
     .           ( harmave(ni(ixrb(jx)+1,1,ifld),ni(ixrb(jx),0,ifld))
     .                       - ni(ixrb(jx)+1,0,ifld) ) / n0(ifld)
           endif
         enddo
      
  274 continue  # ifld loop over ion species





c...  Now do the parallel velocity BC at iy = 0 for full mom eqn species
      do ifld = 1, nusp
         do ix = i4+1-ixmnbcl, i8-1+ixmxbcl  #large ix-loop for up
            if (isuponxy(ix,0,ifld) .eq. 1) then
               ix2 = ixp1(ix,0)
               iv2 = idxu(ix,0,ifld)
               if (isixcore(ix)==1) then   # ix is part of the core boundary
                  upbdry = up(ix,0,ifld)
                  upbdry1 = up(ix,1,ifld)
                  upbdry2 = up(ix,2,ifld)
	          if (isupcore(ifld) == 0) then
		    yldot(iv2) = nurlxu * (upcore(ifld) - upbdry)/vpnorm
		  elseif (isupcore(ifld) ==1) then # d(up)/dy=0
                    yldot(iv2) = nurlxu * (upbdry1 - upbdry) / vpnorm
		  elseif(isupcore(ifld) ==2) then  # d^2(up)/dy^2 = 0
		    yldot(iv2) = nurlxu * ( (upbdry1-upbdry )*gy(ix,1) -
     .                                    (upbdry2-upbdry1)*gy(ix,2) )/
     .                                                 (gy(ix,1)*vpnorm)
		  elseif(isupcore(ifld)==3) then # set fmiy=0 on core bdry
                     yldot(iv2) = -nurlxu * fmiy(ix,0,ifld) / 
     .                                    (vpnorm*sy(ix,0)*fnorm(ifld))
		  elseif (isupcore(ifld)==4 .or. isupcore(ifld)==5) then 
                      # n*uz/rm=const & for =4- radial flux=lzflux; 
                      # or, for =5- ave uz=utorave 
                    yldot(iv2)=-nurlxe*( uz(ix,0,ifld)-uz(ixp1(ix,0),0,ifld)*
     .                           rm(ix,0,2)*ni(ixp1(ix,0),0,ifld)/
     .                          (rm(ixp1(ix,0),0,2)*ni(ix,0,ifld)) )/vpnorm
                    if (ix.eq.ix_fl_bc) then # not all Jac elems included
                      # sum core parallel momentum flux:
                      ii = max(0, ixpt1(1)+1)
                      fmiytotc = rm(ii,0,0)*fmiy(ii,0,ifld)
                      sytotc = sy(ii,0)
                      uztotc = uz(ii,0,ifld)/gxf(ii,0)
                      xtotc = dx(ii,0)
                      do    # loop over ii as changed by ii=ixp1 statement
                        ii = ixp1(ii,0)
                        if (isudsym==1 .and. (ii==nxc .or. ii==nxc+1)) then
                          continue  # skip m.p. guard cells
                        else
                          fmiytotc=fmiytotc+ rm(ii,0,0)*fmiy(ii,0,ifld)
                          sytotc = sytotc + sy(ii,0)
                          uztotc = uztotc + uz(ii,0,ifld)/gxf(ii,0)
                          xtotc = xtotc + dx(ii,0)
                        endif
                        if (ii==ix_fl_bc) break
                      enddo
                      uztotc = uztotc/xtotc   # divide by total core distance
                      if (isupcore(ifld)==4) then
                         yldot(iv2)= -nurlxe*( fmiytotc/sytotc -
     .                                                 lzflux(ifld) )/
     .                                (vpnorm*fnorm(ifld)*rm(ixmp,0,0))
                      elseif (isupcore(ifld)==5) then
                         yldot(iv2)= -nurlxe*(uztotc - utorave(ifld))/vpnorm
                      endif
                    endif
                  else   # isupcore does not correspond to cases above
                    call xerrab ("*** Illegal setting of isupcore ***")
                  endif

               elseif (isupwiix(ix,ifld)==1) then   # isixcore if-test; PF bdry
                  yldot(iv2) = -nurlxu * fmiy(ix,0,ifld) / 
     .                                    (vpnorm*sy(ix,0)*fnorm(ifld))
               elseif (isupwiix(ix,ifld)==2) then   # PF bdry
                  yldot(iv2) = nurlxu * nm(ix,0,ifld) / fnorm(ifld) *
     .                               (up(ix,1,ifld) - up(ix,0,ifld))
               elseif (isupwiix(ix,ifld)==3) then   # PF bdry
                  yldot(iv2) = -nurlxu * nm(ix,0,ifld) / fnorm(ifld) *
     .                         ( up(ix,0,ifld) - up(ix,1,ifld)*
     .                            (2*gyf(ix,0)*lyup(1)-1)/
     .                                     (2*gyf(ix,0)*lyup(1)+1) )
               else                          # PF bdry
                  yldot(iv2) = nurlxu * nm(ix,0,ifld) / fnorm(ifld) *
     .                                      (0. - up(ix,0,ifld))
               endif   # end if-test for core and p.f. boundaries
            endif      # endif for if (isuponxy .eq. 1)
	 enddo         # ix-loop for up eqns
      enddo

c ... Do boundary conditions for impurities along iy = 0.
c     Force corner cells (0,0) and (nx+1,0) to relax to an average of
c     the adjacent cells.
      do ifld = 1, nzspt
         iimp = nhsp + ifld
         if (isimpon .ge. 3 .and. isimpon .le. 7) then 
           do ix = i2, i5   # ix-loop for ion dens
             if (isnionxy(ix,0,iimp)==1) then
               iv = idxn(ix,0,iimp)
               if (isixcore(ix)==1) then   # ix is part of the core boundary:
                  if (isnicore(iimp) .eq. 1) then # set uniform, fixed density
                     yldot(iv) = nurlxn*(ncore(iimp)-ni(ix,0,iimp))/
     .                                                         n0(iimp)
                  elseif (isnicore(iimp) .eq. 0) then # set flux to curcore
                     yldot(iv) = -nurlxn * ( qe*
     .                           (fniy(ix,0,iimp)-fniycbo(ix,iimp))/sy(ix,0)
     .                             - curcore(iimp)*gyf(ix,0)/sygytotc ) /
     .                                               (qe*vpnorm*n0(iimp))
                  elseif (isnicore(iimp) .eq. 2) then # set flux & ni over range
                     yldot(iv) = -nurlxn * ( qe*
     .                           (fniy(ix,0,iimp)-fniycbo(ix,iimp))/sy(ix,0)
     .                             - curcore(iimp)*gyf(ix,0)/sygytotc ) /
     .                                               (qe*vpnorm*n0(iimp))
                     if (ix.ge.ixfixnc .and. ix.le.ixfixnc+incixc) then
                        yldot(iv) = nurlxn*(ncore(iimp)-ni(ix,0,iimp))/
     .                                                         n0(iimp)
                     endif
                  elseif (isnicore(iimp) .eq. 3) then # set integr. flux, const ni 
                     yldot(iv1) = nurlxn*( ni(ix,0,iimp) -
     .                                 ni(ixp1(ix,0),0,iimp) ) / n0(iimp)
                     if (ix.eq.ix_fl_bc) then
                        ii=max(0,ixpt1(1)+1)
                        fniytotc=fniy(ii,0,iimp)-fniycbo(ii,iimp)
                        do
                          ii=ixp1(ii,0)
                          fniytotc=fniytotc+fniy(ii,0,iimp)-fniycbo(ii,iimp)
                          if (ii==ix_fl_bc) break
                        enddo
                        yldot(iv1)= -nurlxn*(qe*fniytotc-
     .                               curcore(iimp))/(qe*vpnorm*n0(iimp))
                     endif
                  elseif (isnicore(iimp) .eq. 4) then  # use bndry impur source
                     yldot(iv) = -nurlxn *
     .                           (fniy(ix,0,iimp) - fnzysi(ix,ifld)) /
     .                           (sy(ix,0) * n0(iimp) * vpnorm)
                  elseif (isnicore(iimp) == 5) then # set dni/dy=0
                     yldot(iv1) = nurlxn*( ni(ix,1,iimp) -
     .                                       ni(ix,0,iimp) ) / n0(ifld)
                  endif
               else   # ix is not part of the core boundary:
                  if (isnwconiix(ix,iimp) .eq. 0) then
                     yldot(iv) = -nurlxn *
     .                           (fniy(ix,0,iimp) - fnzysi(ix,ifld)) /
     .                           (sy(ix,0) * n0(iimp) * vpnorm)
                  elseif (isnwconiix(ix,iimp) .eq. 1) then
                     yldot(iv) = nurlxn*(nwalli(ix) - ni(ix,0,iimp))/
     .                                                        n0(iimp)
                  elseif (isnwconiix(ix,iimp) .eq. 2) then  #extrapolation bc
                     nbound =  ni(ix,1,iimp) - gyf(ix,1)*
     .                         (ni(ix,2,iimp)-ni(ix,1,iimp))/gyf(ix,0)
                     nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                          ni(ix,1,iimp)-1) ) ) + 0.2*ni(ix,1,iimp)
                     yldot(iv) = nurlxn *(nbound - ni(ix,0,iimp))/n0(iimp)
                  elseif (isnwconiix(ix,iimp) .eq. 3) then  #spec. gradient
                     yldot(iv) = -nurlxn*( niy0(ix,0,iimp) -
     .                 niy1(ix,0,iimp)*(2*gyf(ix,0)*lyniix(1,ix,iimp)-1)/
     .                                 (2*gyf(ix,0)*lyniix(1,ix,iimp)+1) -
     .                                     nwimin(iimp) ) / n0(iimp)
                  endif
               endif   # end if-test for core and p.f. boundaries
             endif     # end if-test for isnionxy
       enddo       # end for do ix = i2, i5 loop for ion dens
c...   set the corner values here
            if (ixmnbcl*iymnbcl.eq.1) then
               do jx = 1, nxpt
          if (isnionxy(ixlb(jx),0,iimp)==1) then
                    iv = idxn(ixlb(jx),0,iimp)
                    yldot(iv) = nurlxn *
     .                 ( harmave(ni(ixlb(jx),1,iimp),ni(ixlb(jx)+1,0,iimp))
     .                           - ni(ixlb(jx),0,iimp) ) / n0(iimp)
                  endif
               enddo
            endif
            if (ixmxbcl*iymnbcl.eq.1) then
               do jx = 1, nxpt
             if (isnionxy(ixrb(jx)+1,0,iimp)==1) then
                   iv = idxn(ixrb(jx)+1,0,iimp)
                   yldot(iv) = nurlxn *
     .                 ( harmave(ni(ixrb(jx)+1,1,iimp),ni(ixrb(jx),0,iimp))
     .                           - ni(ixrb(jx)+1,0,iimp) ) / n0(iimp)
                 endif
               enddo
            endif

         endif
      enddo          # end of impurities along iy = 0

CCC  Do Te, Ti & later diffusive neutrals, phi along iy=0 boundary
ccc  -------------------------------------------------------------------
c
c    if flux energy condition, precompute feeytotc & feiytotc for use
      if (iflcore==1) then  #sum e-i core pwr fee,iytotc; unnecess if te,i frozen
        feeytotc = 0.
        feiytotc = 0.
        do ix = i4+1-ixmnbcl, i8-1+ixmxbcl
          if (iymnbcl==1 .and. isixcore(ix)==1) then  # domain part of core bdry
            ii = max(0, ixpt1(1)+1)
            feeytotc = feey(ii,0)-feeycbo(ii)
            feiytotc = feiy(ii,0)-feiycbo(ii)
            do    # loop over ii as changed by ii=ixp1 statement
              ii = ixp1(ii,0)
              if (isudsym==1 .and.(ii==nxc .or. ii==nxc+iexclnxc1)) then
                 continue  # skip m.p. guard cells
              else
                 feeytotc = feeytotc+feey(ii,0)-feeycbo(ii)
                 feiytotc = feiytotc+feiy(ii,0)-feiycbo(ii)
              endif
              if (ii==ix_fl_bc) break
            enddo
          endif
        enddo
      endif
c ... Sum for serial case
      ispwrbcl = 0
      if (ndomain == 1) then
        sfeeytotc = feeytotc
        sfeiytotc = feiytotc
        ispwrbcl = 1
      else
        if (ispwrbc(mype+1) == 1) ispwrbcl = 1
      endif

c     Sum for parallel case, sum/distribute to all processors
c_mpi      call mpi_allreduce(feeytotc, sfeeytotc, 1,
c_mpi     .      MPI_DOUBLE_PRECISION, MPI_SUM, uedgeComm, ierr)
c_mpi      call mpi_allreduce(feiytotc, sfeiytotc, 1,
c_mpi     .      MPI_DOUBLE_PRECISION, MPI_SUM, uedgeComm, ierr)

c ... Set Te and Ti BCs
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl #long ix-loop for Te, Ti
         if (isteonxy(ix,0) .eq. 1) then
           iv1 = idxte(ix,0)
            if (isixcore(ix)==1) then   # ix is part of the core boundary:
              yldot(iv1)=nurlxe*(tcoree*ev-te(ix,0))*1.5*ne(ix,0)/ennorm
              if (iflcore .eq. 1) then
                 yldot(iv1) = -nurlxe*(te(ix,0)-te(ixp1(ix,0),0))*n0(1)/ennorm
                 if (ix.eq.ix_fl_bc) then # not all Jac elems included
                    # sum core power:
                    ii = max(0, ixpt1(1)+1)
                    feeytotc = feey(ii,0)-feeycbo(ii)
                    do    # loop over ii as changed by ii=ixp1 statement
                       ii = ixp1(ii,0)
                       if (isudsym==1 .and. (ii==nxc .or. ii==nxc+iexclnxc1)) then
                          continue  # skip m.p. guard cells
                       else
                          feeytotc = feeytotc+feey(ii,0)-feeycbo(ii)
                       endif
                       if (ii==ix_fl_bc) break
                    enddo
                    yldot(iv1)= -nurlxe*(feeytotc-pcoree)/(vpnorm*ennorm)
                 endif
              elseif (iflcore .eq. -1) then
                 yldot(iv1) = -nurlxe*(te(ix,0)-te(ix,1))*n0(1)/ennorm
              endif
            else   # ix is not part of the core boundary; various PF cases:
               if (istepfcix(ix) .eq. 0) then          # zero electron energy flux
                 yldot(iv1) = -nurlxe*(feey(ix,0)/(n0(1)*vpnorm*sy(ix,0)))
     .                                                       / (temp0*ev) 
               elseif (istepfcix(ix) .eq. 1) then      # fixed Te
                 yldot(iv1) =nurlxe*(tewalli(ix)*ev-te(ix,0))/(temp0*ev)
               elseif (istepfcix(ix) .eq. 2) then       # extrapolation
                 tbound = te(ix,1) - gyf(ix,1)*
     .                                     (te(ix,2)-te(ix,1))/gyf(ix,0)
                 tbound = max(tbound, tbmin*ev)
                 yldot(iv1) = nurlxe *(tbound - te(ix,0))/(temp0*ev)
               elseif (istepfcix(ix) .eq. 3) then         # specified gradient
                 yldot(iv1) = nurlxe*( (te(ix,1) - te(ix,0)) -
     .                             0.5*(te(ix,1) + te(ix,0))/
     .                              (gyf(ix,0)*lyteix(1,ix)) )/(temp0*ev)
               elseif (istepfcix(ix) .eq. 4) then      
                 yldot(iv1) = -nurlxe*(feey(ix,0) - bceew*ne(ix,0)*
     .                               vey(ix,0)*sy(ix,0)*te(ix,0)) /
     .                              (vpnorm*ennorm*sy(ix,0))
               endif
            endif   # end if-test for core and p.f. boundaries
         endif   # end if-test on isteonxy

         if (istionxy(ix,0) .eq. 1) then
            iv2 = idxti(ix,0)
            if (isixcore(ix)==1) then   # ix is part of the core boundary:
              yldot(iv2)=nurlxi*(tcorei*ev-ti(ix,0))*1.5*ne(ix,0)/ennorm
              if (iflcore .eq. 1) then
                 yldot(iv2) = -nurlxi*(ti(ix,0)-ti(ixp1(ix,0),0))*
     .                                                     n0(1) /ennorm
                 if (ix.eq.ix_fl_bc) then # not all Jac elems included
                    # sum core power:
                    ii = max(0, ixpt1(1)+1)
                    feiytotc = feiy(ii,0)-feiycbo(ii)
                    do
                       ii = ixp1(ii,0)
                       if (isudsym==1 .and. (ii==nxc .or. ii==nxc+iexclnxc1)) then
                          continue  # skip m.p. guard cells
                       else
                          feiytotc = feiytotc+feiy(ii,0)-feiycbo(ii)
                       endif
                       if (ii==ix_fl_bc) break
                    enddo
                    yldot(iv2)= -nurlxi*(feiytotc-pcorei)/(vpnorm*ennorm)
                 endif
              elseif (iflcore .eq. -1) then
                 yldot(iv2) = -nurlxe*(ti(ix,0)-ti(ix,1))*n0(1)/ennorm
              endif
ccc            elseif (matwalli(ix) .gt. 0) then
ccc               bcen = matwalli(ix)
ccc               yldot(iv2) = -nurlxi*(feiy(ix,0) -
ccc     .             (bcei*fniy(ix,0,1) + bcen*fniy(ix,0,iigsp))*ti(ix,0)) /
ccc     .             (vpnorm*ennorm*sy(ix,0))
            else                                 # various PF cases
               if (istipfcix(ix) .eq. 0) then          # zero ion energy flux
                 yldot(iv2) = -nurlxi*(feiy(ix,0)/(n0(1)*vpnorm*sy(ix,0)))
     .                                                       / (temp0*ev) 
               elseif (istipfcix(ix) .eq. 1) then      # fixed Ti
                 yldot(iv2) =nurlxi*(tiwalli(ix)*ev-ti(ix,0))/(temp0*ev)
               elseif (istipfcix(ix) .eq. 2) then      # extrapolation
                 tbound = ti(ix,1) - gyf(ix,1)*
     .                                     (ti(ix,2)-ti(ix,1))/gyf(ix,0)
                 tbound = max(tbound, tbmin*ev)
                 yldot(iv2) = nurlxi *(tbound - ti(ix,0))/(temp0*ev)
               elseif (istipfcix(ix) .eq. 3) then      # specified gradient
                 yldot(iv2) = nurlxi*( (ti(ix,1) - ti(ix,0)) -
     .                           0.5*(ti(ix,1) + ti(ix,0))/

     .                             (gyf(ix,0)*lytiix(1,ix)) )/(temp0*ev)
               elseif (istipfcix(ix) .eq. 4) then      # sheath-like condition
                 t0 = max(tg(ix,0,1), temin*ev)
                 fngyw=0.25*sqrt(8*t0/(pi*mg(1)))*ng(ix,1,1)*sy(ix,0)
                 yldot(iv2) = -nurlxi*( feiy(ix,0) - bceiw*fniy(ix,0,1)*
     .                                                        ti(ix,0) -
     .                     cftiexclg*( bcenw*fniy(ix,0,iigsp)*tg(ix,0,1) -
     .                                cgengw*2.*tg(ix,0,iigsp)*fngyw ) ) /
     .                                          (vpnorm*ennorm*sy(ix,0))
               endif
            endif   # end if-test for core and p.f. boundaries
         endif   # end if-test on istionxy
      enddo  # long ix-loop for Te and Ti BCs

ccc  Diffusive neutrals on iy=0 boundary
ccc  - - - - - - - - - - - - -
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl #long ix-loop for diff neut dens
         nzsp_rt = nhsp
         do 275 igsp = 1, ngsp
       jz = max(igsp - nhgsp, 1)  # identify impurity index
       if (jz > 1) nzsp_rt = nzsp_rt + nzsp(jz-1) #prev index for fniy
           if (isngonxy(ix,0,igsp) .eq. 1) then  # ends just before 275 continue statem.
            if (iscpli(ix) .eq. 1) call wsmodi(igsp)
            iv = idxg(ix,0,igsp)   
c ... prepare impurity ion flux for possible recycling
            zflux = 0.
            if (igsp .gt. nhgsp) then
              do iimp = 1, nzsp(jz)
        zflux = zflux + min(fniy(ix,0,nzsp_rt+iimp), 0.)
              enddo
        endif
c ... Set core BC
            if (isixcore(ix)==1) then   # ix is part of the core boundary:
               if (isngcore(igsp) .eq. 0) then
                 osmw = onesided_maxwellian( 
     .                  cdifg(igsp)*tg(ix,0,igsp), 1.0, mg(igsp), 
     .                  sy(ix,0), tgmin*ev
     .           )
                 yldot(iv) = -nurlxg*(fngy(ix,0,igsp) 
     .              + (1-albedoc(igsp))*osmw
     .              * harmave(ng(ix,0,igsp), ng(ix,1,igsp))
     .              ) / (osmw*n0g(igsp))
               elseif (isngcore(igsp).eq.1) then
                 yldot(iv) = nurlxg*(ngcore(igsp)-ng(ix,0,igsp)) /
     .                                                       n0g(igsp)
               elseif (isngcore(igsp) .eq. 2) then
                 lengg = sqrt(tg(ix,0,igsp)/
     .                  (mg(igsp)*(nuix(ix,0,igsp)*nuiz(ix,0,igsp))) )
                 yldot(iv) = nurlxn*(
     .                         (ng(ix,1,igsp) - ng(ix,0,igsp)) -
     .                        0.5*(ng(ix,1,igsp) + ng(ix,0,igsp))/
     .                                 (gyf(ix,0)*lengg) )/n0g(igsp)
               elseif (isngcore(igsp) .eq. 3) then
                 nbound =  ng(ix,1,igsp) - gyf(ix,1)*
     .                        (ng(ix,2,igsp)-ng(ix,1,igsp))/gyf(ix,0)
                 nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                         ng(ix,1,igsp)-1) ) ) + 0.2*ng(ix,1,igsp)
                 yldot(iv) = nurlxn *(nbound - ng(ix,0,igsp))/n0g(igsp)
               elseif (isngcore(igsp) .eq. 4) then
                 yldot(iv) = nurlxn *(ng(ix,1,igsp) - ng(ix,0,igsp))/
     .                                                       n0g(igsp)
               endif

###  End of core bounary conditions 

            else  # cell is not part of core boundary
### wall BCs. calc diff neut sputt fluxes
               fng_chem = 0.
               do igsp2 = 1, 1+ishymol  #include hydrogen neut sputt
                 t0 = max(tg(ix,0,igsp2),tgmin*ev)
                 flx_incid=ng(ix,1,igsp2)*0.25*sqrt(8*t0/(pi*mg(igsp2)))
		 if (isch_sput(igsp).ne.0) then
                   call sputchem (isch_sput(igsp), t0/ev, tvwalli(ix),
     .                                       flx_incid, yld_carbi(ix))
                   fng_chem = fng_chem + fchemygwi(igsp2)*yld_carbi(ix)*
     .                                                flx_incid*sy(ix,0)
                 else
                   fng_chem = fng_chem + chemsputi(igsp,igsp2)* 
     .                                                flx_incid*sy(ix,0)
                 endif
               enddo






c ... Include gas BC from sputtering by ions
               sputflxpf(ix,igsp) = 0.
               if (isi_sputpf(igsp).ge.1 .and. igsp>1) then # incl phys sput from ions
                 do ifld = ipsputt_s, ipsputt_e 
                   eng_sput = ( 2*ti(ix,0) + zi(ifld)*
     .                                           3*te(ix,0) )/ev
                   if (zi(ifld) > 0.) sputflxpf(ix,igsp) = sputflxpf(ix,igsp) +
     .                                min(fniy(ix,0,ifld),0.)*fphysyiwi(igsp)*
     .                                      yld96(matp,matt,eng_sput)
                 enddo
                  if (isi_sputpf(igsp) .ge. 2) then # add chem sput from ions
                    do ifld = 1,1  # (ipsputt_s, ipsputt_e) place holder for imp/imp sputt
                      eng_sput = ( 2*ti(ix,0) + zi(ifld)*
     .                                             3*te(ix,0) )/ev
              flx_incid = -min(fniy(ix,0,ifld)/sy(ix,0),0.)
                      call sputchem (isch_sput(igsp), eng_sput, tvwalli(ix),
     .                                            flx_incid, yld_chm)
                      sputflxpf(ix,igsp) = sputflxpf(ix,igsp)+fchemyiwi(igsp)*
     .                                        min(fniy(ix,0,ifld),0.)*yld_chm
                    enddo
                  endif
               endif
               osmw = onesided_maxwellian( 
     .                  cdifg(igsp)*tg(ix,0,igsp), 1.0, mg(igsp), 
     .                  sy(ix,0), tgmin*ev
     .         )
               yldot(iv) = -nurlxg*( fngy(ix,0,igsp) 
     .              + (1-albedoi(ix,igsp))*osmw
     .              * harmave(ng(ix,0,igsp), ng(ix,1,igsp))
     .              - fng_chem + sputflxpf(ix,igsp) ) / (osmw*n0g(igsp)) 
              



               if(matwalli(ix) .gt. 0) then
                 if (recycwit(ix,igsp,1) .gt. 0.) then
                  


                   fniy_recy = fac2sp*fniy(ix,0,1)
                   if(isrefluxclip==1) 
     .                  fniy_recy=min(fniy_recy,0.)
                   if (igsp .gt. nhgsp) 
     .                  fniy_recy = zflux
                   if (ishymol.eq.1 .and. igsp.eq.2) then # 2 atoms per molecule
                     if (isupgon(1) .eq. 1) then
                        fniy_recy = 0.5*( fniy(ix,0,1) + fniy(ix,0,2) )
                     else
                       fniy_recy = 0.5*( fniy(ix,0,1) + fngy(ix,0,1) )
                     endif
                     if(isrefluxclip==1) fniy_recy=min(fniy_recy,0.)
                   endif
                   osmw = onesided_maxwellian( 
     .                  cdifg(igsp)*tg(ix,0,igsp), 1.0, mg(igsp), 
     .                  sy(ix,0), tgmin*ev
     .             ) 
                   yldot(iv) = -nurlxg*( fngy(ix,0,igsp) 
     .                  + fniy_recy*recycwit(ix,igsp,1) 
     .                  - fngyi_use(ix,igsp) - fngysi(ix,igsp) 
     .                  + (1-albedoi(ix,igsp))*osmw
     .                  * harmave(ng(ix,0,igsp), ng(ix, 1, igsp))
     .                  - fng_chem + sputflxpf(ix,igsp)
     .                  ) / (osmw*n0g(igsp))
                
                elseif (recycwit(ix,igsp,1) < -1) then
                   yldot(iv)=nurlxg*(ngbackg(igsp)-ng(ix,0,igsp))/
     .                                                          n0g(igsp)
                 elseif (recycwit(ix,igsp,1) .le. 0.) then # treat recycwit as albedo
                   osmw = onesided_maxwellian(
     .                  cdifg(igsp)*tg(ix,0,igsp), 1.0, mg(igsp), 
     .                  sy(ix,0), tgmin*ev
     .             )
                   yldot(iv) = -nurlxg*( fngy(ix,0,igsp)
     .                  + (1+recycwit(ix,igsp,1))*osmw
     .                  * harmave(ng(ix,0,igsp), ng(ix,1,igsp))
     .                  ) / (osmw*n0g(igsp))
                 endif 
               endif

               if(  (fngysi(ix,igsp)+fngyi_use(ix,igsp) .ne.0.) 
     .              .and. (matwalli(ix).eq.0.)
     .         ) then
                    osmw = onesided_maxwellian( 
     .                  cdifg(igsp)*tg(ix,0,igsp), 1.0, mg(igsp), 
     .                  sy(ix,0), tgmin*ev 
     .              )
                    yldot(iv) = -nurlxg*( fngy(ix,0,igsp) 
     .                  - fngysi(ix,igsp) ) / (osmw*n0g(igsp))
               endif 
            endif   # end if-test for core and p.f. boundaries
           endif    # End for if(isngon(igsp).eq.1) after do 275
 275     continue   # igsp loop over species
      enddo  # end of long ix-loop for diff neuts

c... BC for neutral gas temperature/energy at iy=0
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl # ix-loop for diff neut temp 
        do igsp = 1, ngsp
          if (istgonxy(ix,0,igsp) == 1) then
            iv = idxtg(ix,0,igsp)
            if (isixcore(ix)==1) then
              if (istgcore(igsp) == 0) then  #set to ti
                yldot(iv)=nurlxg*(ti(ix,0)*cftgticore(igsp) 
     .                            -tg(ix,0,igsp))/(temp0*ev)
              elseif (istgcore(igsp) == 1) then   # set to tgcore
                yldot(iv)=nurlxg*(tgcore(igsp)*ev-tg(ix,0,igsp))/
     .                                                      (temp0*ev)
              elseif (istgcore(igsp) == 2) then 
                  osmw = onesided_maxwellian(
     .                  tg(ix,0,igsp), harmave(ng(ix,0,igsp), ng(ix,1,igsp)),
     .                  mg(igsp), sy(ix,0), tgmin*ev
     .            )
                  yldot(iv)=-nurlxg*(fegy(ix,0,igsp) 
     .              + cfalbedo*max(tg(ix,0,igsp),tgmin*ev)*(1-albedoc(igsp))*osmw
     .              ) / (vpnorm*ennorm*sy(ix,0))
              else # all others, set to zero y-gradient
                yldot(iv)=nurlxg*(tg(ix,1,igsp)-tg(ix,0,igsp))/
     .                                                      (temp0*ev)
              endif
            else   # PF wall
              if (istgpfcix(ix,igsp) == 0) then   # fixed Tg
                yldot(iv) = nurlxg*(tgwall(igsp)*ev-tg(ix,0,igsp))/
     .                                                      (temp0*ev)
              elseif (istgpfcix(ix,igsp) == 1)    # extrapolation
                tbound = tg(ix,1,igsp) - gyf(ix,1)*
     .                      (tg(ix,2,igsp)-tg(ix,1,igsp))/gyf(ix,0)
                tbound = max(tbound, 0.25*tbmin*ev)  #tbmin=.1 eV
                yldot(iv) = nurlxi *(tbound - tg(ix,0,igsp))/(temp0*ev)
              elseif (istgpfcix(ix,igsp) == 2)    # specified gradient
                yldot(iv) = nurlxi*( (tg(ix,1,igsp) - tg(ix,0,igsp)) -
     .                           0.5*(tg(ix,1,igsp) + tg(ix,0,igsp))/
     .                           (gyf(ix,0)*lytg(1,igsp)) )/(temp0*ev)
              elseif (istgpfcix(ix,igsp) == 3)  #Maxwell thermal flux to wall
                osmw = onesided_maxwellian(
     .                  cdifg(igsp)*tg(ix,1,igsp), ng(ix,1,igsp), 
     .                  mg(igsp), sy(ix,0), tgmin*ev
     .          )
                yldot(iv) =  -nurlxg*( fegy(ix,0,igsp) 
     .              + 2*cgengmw*osmw) / (sy(ix,0)*vpnorm*ennorm)
              elseif (istgpfcix(ix,igsp) == 4) 
                fng_chem = 0.  #..double check
                osmw = onesided_maxwellian(
     .              tg(ix,0,igsp), harmave(ng(ix,0,igsp), ng(ix,1,igsp)),
     .              mg(igsp), sy(ix,0), tgmin*ev
     .          )
                yldot(iv)=-nurlxg*(fegy(ix,0,igsp) 
     .              + cfalbedo*max(tg(ix,0,igsp),tgmin*ev)*(1-albedoi(ix,1))*osmw
     .              - 2.*fng_chem*max(tg(ix,0,igsp), tgmin*ev)) / (vpnorm*ennorm*sy(ix,0))
                fniy_recy = 0.
                if (matwalli(ix) .gt. 0) then
                  if (recycwit(ix,igsp,1) .gt. 0) then
                    fniy_recy = recycwit(ix,igsp,1)*fac2sp*fniy(ix,0,1)
                    if (isrefluxclip==1) 
     .                  fniy_recy=min(fniy_recy,0.)
                    osmw = onesided_maxwellian(
     .                  tg(ix,0,igsp), harmave(ng(ix,0,igsp), ng(ix,1,igsp)),
     .                  mg(igsp), sy(ix,0), tgmin*ev
     .              )
                    yldot(iv)=-nurlxg*(fegy(ix,0,igsp) 
     .                  + cfalbedo*max(tg(ix,0,igsp),tgmin*ev)*(1-albedoi(ix,1))*osmw
     .                  - 2.*fng_chem*max(tg(ix,0,igsp),tgmin*ev)
     .                  + fniy_recy*(1.-cfdiss)*cfalbedo*recycwe*ti(ix,0)
     .                  ) / (vpnorm*ennorm*sy(ix,0))
                  endif
                endif
              elseif (istgpfc(igsp) == 5) then   # tg=ti*cftgtipfc
                yldot(iv) = nurlxg*(ti(ix,0)*cftgtipfc(igsp)
     .                                        -tg(ix,0,igsp))/
     .                                                      (temp0*ev)
              elseif (istgpfc(igsp) > 5)
                 call xerrab("***Input error: invalid istgpfc ***")
              endif

            endif           
          endif
        enddo  
      enddo  # end ix-loop for neutral temps

ccc  Now do potential along iy=0 boundary if isnewpot=-1 or -2; others next page
ccc  - - - - - - - - - - - - - - 
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl #long ix-loop for potential
         if (isphionxy(ix,0) .eq. 1) then
            iv3 = idxphi(ix,0)
            yldot(iv3) = nurlxp*( (phi(ix,1) - phi(ix,0)) -
     .                             0.5*(phi(ix,1) + phi(ix,0))/
     .                              (gyf(ix,0)*lyphiix(1,ix)) )/temp0
            if (isnewpot==-1) then
               if (isixcore(ix)==1) then   # ix is part of the core boundary:
                  ix4 = ixm1(ix,0)
                  ut0 = (rr(ix,0)+rr(ix,1))/(btot(ix,0)+btot(ix,1))*
     .                  ( ey(ix,0) - gpiy(ix,0,1)/
     .                         (qe*zi(1)*niy0(ix,0,1)) ) +
     .                        0.25*(rbfbt(ix,0)+rbfbt(ix,1))*
     .                              (up(ix,0,1)+up(ix4,0,1))
                  yldot(iv3) = -nurlxp*ut0/vpnorm
               endif   # end if-test for core and p.f. boundaries
            endif
            if (isnewpot==-2) then  # phi is constant on core boundary
               if (isixcore(ix)==1) then   # ix is part of the core boundary:
                  if (ix .ne. ix_fl_bc) then
                     yldot(iv3) = -nurlxp*(phi(ix,0)-phi(ixp1(ix,0),0))/temp0
                  elseif (ix==ix_fl_bc) then # total core current = icoreelec
                     fqytotc = 0.
                     ii=ixc1
                     fqytotc=fqyao(ii,0)
                     do
                        ii=ixp1(ii,0)
                        fqytotc = fqytotc + fqyao(ii,0)
                        if (ii==ix_fl_bc) break
                     enddo
                     yldot(iv3) = -nurlxp*(fqytotc - icoreelec)/
     .                                     (qe*n0(1)*vpnorm*sy(ixc1,0))
                  endif   # end if-test for integral BC on core
               endif   # end if-test for core and p.f. boundaries
            endif   # end if-test for isnewpot==-2
         endif

      enddo # long ix loop for potential
ccc  -------------------------------------------------------------------

c...  Now do the corners: the ion density corners are done elsewhere
       if ((xcnearlb.or.openbox) .and. ixmnbcl*iymnbcl.eq.1) then  
         do ifld = 1, nusp
            do jx = 1, nxpt
	     if(isuponxy(ixlb(jx),0,ifld)==1) then
               yldot(idxu(ixlb(jx),0,ifld)) = - nurlxu*
     .                      ( up(ixlb(jx),0,ifld) -
     .             0.5*(up(ixlb(jx),1,ifld)+up(ixlb(jx)+1,0,ifld)) )
     .                                         /vpnorm
             endif
           enddo
         enddo
         do jx = 1, nxpt
	     if (isteonxy(ixlb(jx),0)==1) then
               yldot(idxte(ixlb(jx),0)) = nurlxe*
     .             ( 0.5*(te(ixlb(jx)+1,0)+te(ixlb(jx),1))
     .                               - te(ixlb(jx),0) ) / (temp0*ev)
             endif
         enddo
         do jx = 1, nxpt
	    if (istionxy(ixlb(jx),0)==1) then
               yldot(idxti(ixlb(jx),0)) = nurlxi*
     .             ( 0.5*(ti(ixlb(jx)+1,0)+ti(ixlb(jx),1))
     .                               - ti(ixlb(jx),0) ) / (temp0*ev)
            endif
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
	     if (isngonxy(ixlb(jx),0,igsp)==1) then
               yldot(idxg(ixlb(jx),0,igsp)) = nurlxg*
     .           (ng(ixlb(jx)+1,0,igsp) - ng(ixlb(jx),0,igsp))/n0g(igsp)
             endif
           enddo
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
	     if(istgonxy(ixlb(jx),0,igsp)==1) then
               yldot(idxtg(ixlb(jx),0,igsp)) = nurlxg*
     .          (tg(ixlb(jx)+1,0,igsp) - tg(ixlb(jx),0,igsp))/(temp0*ev)
             endif
           enddo
         enddo
       endif
       if ((xcnearrb.or.openbox) .and. ixmxbcl*iymnbcl.eq.1) then
         do ifld = 1, nusp
           do jx = 1, nxpt
	     if(isuponxy(ixrb(jx),0,ifld)==1) then
               yldot(idxu(ixrb(jx),0,ifld)) = - nurlxu*
     .                      ( up(ixrb(jx),0,ifld) -
     .             0.5*(up(ixrb(jx)-1,0,ifld)+up(ixrb(jx),1,ifld)) )
     .                                         /vpnorm
               yldot(idxu(ixrb(jx)+1,0,ifld)) = - nurlxu*
     .                   ( up(ixrb(jx)+1,0,ifld)-up(ixrb(jx),0,ifld) )
     .                                            /vpnorm
             endif
           enddo
         enddo
         do jx = 1, nxpt
	    if(isteonxy(ixrb(jx)+1,0)==1) then
               yldot(idxte(ixrb(jx)+1,0)) = nurlxe*
     .             ( 0.5*(te(ixrb(jx)+1,1)+te(ixrb(jx),0))
     .                               - te(ixrb(jx)+1,0) ) / (temp0*ev)
            endif
         enddo
         do jx = 1, nxpt
	    if(istionxy(ixrb(jx)+1,0)==1) then
               yldot(idxti(ixrb(jx)+1,0)) = nurlxi*
     .             ( 0.5*(ti(ixrb(jx)+1,1)+ti(ixrb(jx),0))
     .                               - ti(ixrb(jx)+1,0) ) / (temp0*ev)
            endif
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
	     if(isngonxy(ixrb(jx),0,igsp)==1) then
               yldot(idxg(ixrb(jx)+1,0,igsp)) = nurlxg*
     .           (ng(ixrb(jx),0,igsp) - ng(ixrb(jx)+1,0,igsp))/n0g(igsp)
             endif
           enddo
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
	     if(istgonxy(ixrb(jx),0,igsp)==1) then
               yldot(idxtg(ixrb(jx)+1,0,igsp)) = nurlxg*
     .          (tg(ixrb(jx),0,igsp) - tg(ixrb(jx)+1,0,igsp))/(temp0*ev)
             endif
           enddo
         enddo
       endif

      endif           # end of iy = 0 boundary conditions except for isnewpot=1





c  ####################################################################
c  ###  Now do the potential for isnewpot=1 case along iy=0 boundary
c  ####################################################################
c...  If isnewpot=1, phi boundary involves two eqns at iy=0 and 1
c...  Note: j3 is local range index for iy passed from pandf in oderhs.m
      if (isnewpot*isphion.eq.1 .and. j3.le.3) then
        do ix = min(i4+1-ixmnbcl,ixpt1(1)+1), max(i8-1+ixmxbcl,ixpt2(nxpt))
          if(isphionxy(ix,0)*isphionxy(ix,1)==1) then
            iv  = idxphi(ix,0)
            iv1 = idxphi(ix,1)

c  ################ do BC on core bdry, then on PF bdry #################
            if (isixcore(ix)==1 .and. isphicore0==0) then # ix is core bdry
               ix3 = ixm1(ix,1)
               ix4 = ixm1(ix,0)
               r_major = 0.5*(rm(ix,0,0) + rm(ix,1,0)) + r0slab
c         Set two B.C. (int. cur. and int. ut0 + pol. uniform phi(,0), te(,1)
c         Set phi poloidally constant at iy=0 on core bdry
                 yldot(iv) = -nurlxp*(phi(ix,0)-phi(ixp1(ix,0),0))/temp0

c         Next set phi at iy=1 depending on iphibcc flag

c  ######################################################################
c  ### First core phi BC at iy=0 & 1 over full ixcore range; then ix=ixmp
c  ######################################################################
cc     First check that iphibcc=1,2, or 3; othewise abort with message
               if (iphibcc < 1 .or. iphibcc > 3) then
               call xerrab("**INPUT ERROR: only iphibcc=1,2,3 available")
               endif

                 if (iphibcc==1) then
                   yldot(iv1) = -nurlxp*( (ey(ix,1)-ey(ix,0))*gy(ix,1) -
     .                                  (ey(ix,2)-ey(ix,1))*gy(ix,2) )/
     .                                                  (gy(ix,1)*temp0)
                 elseif (iphibcc==2) then
            yldot(iv1)=-nurlxp*(te(ix,1)-te(ixp1(ix,1),1))/
     .                                                          (ev*temp0)
                 elseif (iphibcc==3) then
                   yldot(iv1)= -nurlxp*(phi(ix,1)-phi(ixp1(ix,1),1))/temp0
                 else  # no longer allowed; need to reconsider
                   yldot(iv1)= -nurlxp*( phi(ix,1)-phi(ixp1(ix,1),1) +
     .                                      dphi_iy1(ix) ) / temp0
                   if (isteon.eq.0) yldot(iv1)=nurlxp*
     .                                     (tes(ix,1)-te(ix,1))/temp0
                 endif

c  ####################################################################
c  ### Reset phi BC at iy=0,1 at ix=ixmp only using iphibcc & isutcore
c  ####################################################################
               if (ix == ixmp) then  #redefine eqn for midplane ix
                 fqytotc = 0.
                 uztotc = 0.
             uztotc1 = 0.
                 uztotc2 = 0.
c NOTE: total poloidal distance around the core:
                 ii=ixc1
                 xtotc = dx(ii,0)
                 do
                    ii=ixp1(ii,1)
                    xtotc = xtotc + dx(ii,0)
                    if (ii==ix_fl_bc) break
                 enddo
                 ii=ixc1
c NOTE: bndry surface flux should exclude fqyb and include fqyd
                 fqytotc = fqya(ii,1)+ cfqyn*fqyn(ii,1) + cfqym*fqym(ii,1) 
     .                     + cfqybbo*fqyb(ii,1) + cfqydbo*fqyd(ii,1)
                 uztotc  = uz(ii,0,1) / ( gxf(ii,0)*xtotc )
                 uztotc1 = uz(ii,1,1) / ( gxf(ii,1)*xtotc )
                 uztotc2 = uz(ii,2,1) / ( gxf(ii,2)*xtotc )
                 do
                    ii=ixp1(ii,1)
                    fqytotc = fqytotc +
     .                        fqya(ii,1) + cfqyn*fqyn(ii,1) + cfqym*fqym(ii,1) 
     .                           + cfqybbo*fqyb(ii,1) + cfqydbo*fqyd(ii,1)
                    uztotc  = uztotc  + uz(ii,0,1) / ( gxf(ii,0)*xtotc )
                    uztotc1 = uztotc1 + uz(ii,1,1) / ( gxf(ii,1)*xtotc )
                    uztotc2 = uztotc2 + uz(ii,2,1) / ( gxf(ii,2)*xtotc )
                    if (ii==ix_fl_bc) break
                 enddo
c ... Midplane phi(ixmp,0) BC is always to set rad current to icoreelec
                 yldot(iv)=-nurlxp*(fqytotc - icoreelec)/
     .                                     (qe*n0(1)*vpnorm*sy(ixc1,0))

c ... Midplane phi(ixmp,1) BC set by iphibcc if =1,2,3; if iphibbc otherwise,
c ... isutcore sets phi(ixmp,1) BC
         if (iphibcc.eq.1) then  # set d^2(ey)/dy^2=0
                   yldot(iv1) = -nurlxp*( (ey(ix,1)-ey(ix,0))*gy(ix,1) -
     .                                   (ey(ix,2)-ey(ix,1))*gy(ix,2) )/
     .                                                  (gy(ix,1)*temp0)
                 elseif (iphibcc==2 .or. iphibcc==3) then
                   yldot(iv1)=-nurlxp*(ey(ix,0)-eycore)/(gyf(ix,0)*temp0)
                 else  # any iphibcc > 3 uses isutcore BC at ix=ixmp
                   if (isutcore .eq. 0) then # flux-ave R*n*m*uz = lzcore
                     yldot(iv1)=-nurlxp*(uztotc -
     .                        lzcore(1)/(mi(1)*n0(1)*r_major)) / vpnorm
                   elseif (isutcore == 1) then  # zero radial deriv. of ave uz
                     yldot(iv1) = nurlxp*(uztotc-uztotc1)/vpnorm
           else     # second deriv. of ey at the midplane
                     yldot(iv1) = -nurlxp*( (ey(ix,1)-ey(ix,0))*gy(ix,1) -
     .                                   (ey(ix,2)-ey(ix,1))*gy(ix,2) )/
     .                                                  (gy(ix,1)*temp0)
                   endif
                 endif
               endif  #end outer midplane cells ix=ixmp potential BCs
c  ####################################################################
c  ... Completed phi BC on core bdry
c  ####################################################################
c  ... Special case rarely used
               if (cfvycf .gt. 1e-20) then  # preserves PET97 cases
                 yldot(iv) = -nurlxp*( ey(ix,0) - gpiy(ix,0,1)/
     .                     (qe*zi(1)*niy0(ix,0,1)) )/(btot(ix,0)*vpnorm)
                 yldot(iv1) =nurlxp*(fqy(ix,1)-(fqx(ix,1)-fqx(ix3,1))) /
     .                              (rrv(ix,0)*sy(ix,0)*vpnorm*ev*n0(1))
               endif

c  ###################################################################
c  ### Finally set private-flux (PF) phi BC at iy=0 & 1
c  ###################################################################
            else   # ix is not part of the core boundary:
               if(iphibcwiix(ix) == 0) then
                 yldot(iv) = nurlxp*( phi(ix,1) - phi(ix,0) )/temp0
               elseif(iphibcwiix(ix) == 1) then
                 yldot(iv) = nurlxp*
     .                     (phintewi*te(ix,0)/ev - phi(ix,0))/temp0
           elseif(iphibcwiix(ix) == 3) then
                 yldot(iv) = nurlxp*( (phi(ix,1) - phi(ix,0)) -
     .                             0.5*(phi(ix,1) + phi(ix,0))/
     .                               (gyf(ix,0)*lyphiix(1,ix)) )/temp0
               elseif(iphibcwiix(ix) == 4) then
                 yldot(iv) = nurlxp*
     .                     (phiwi(ix) - phi(ix,0))/temp0
               endif
            endif   # end if-test for core and p.f. boundaries
          endif     # end if-test on isphionxy for ix=0,1
        enddo       # end loop over ix
      endif         # large if testing isnewpot*isphion=1, j3

        RETURN
        END SUBROUTINE south_boundary
