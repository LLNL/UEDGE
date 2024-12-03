		SUBROUTINE symmetric_boundary(neq, ix, ixd yl, yldot)
		IMPLICIT NONE
        INTEGER, INTENT(IN) :: neq, ix, ixd
        LOGICAL, INTENT(IN) :: isleftbound
        REAL, INTENT(OUT) :: yl(neq), yldot(neq) 


c 			LEFT
            ix = 0
            ixd = 1
            ix2 = ix
            isfix = isfixlb
            cssging = -1
            yyb = yylb
            alb = alblb
            recyb = recylb
            fngx_use = fngxslb
            fngxs = fngxslb
c 			RIGTH
            ix = ny + 1
            ixd = ny
            ix2 = ixd
            isfix = isfixrb
            cssign = 1
            yyb = yyrb
           	alb = albrb
           	recyb = recyrb
           	fngx_use = fngxrb_use
           	fngxs = fngxsrb



c********************************************************************
c...  First check if ix=0 has fixed boundary values, no potential
c...  isfixlb=1 sets all profiles; isfixlb=2 sets reflection boundary 
c...  conditions
c********************************************************************
c 		LEFT
		if (i3 .le. 0 .and. isfixlb(1) .ne. 0) then
c 		RIGHT
      	if (i6 .ge. nx+1 .and. isfixrb(1) .gt. 0) then  #begin symmetry BC at nx+1


      do iy = j2, j5
         do ifld = 1 , nisp
            if(isnionxy(ix,iy,ifld) .eq. 1) then
               iv1 = idxn(ix,iy,ifld)
               yldot(iv1) = nurlxn *
     .                (nib(ifld)*nibprof(iy)-ni(ix,iy,ifld))/n0(ifld)
               if(isfix(1).eq.2) 
     .				yldot(iv1) = nurlxn * (1/n0(ifld))
     .                           * (ni(ixd,iy,ifld) - ni(ix,iy,ifld))
            endif
         enddo



         do ifld = 1, nusp
            if(isuponxy(ix,iy,ifld) .eq. 1) then
               iv2 = idxu(ix,iy,ifld)
       		if isleftbound then
               yldot(iv2) = nurlxu *
     .           (upb(ifld)*upbprof(iy) - up(ix,iy,ifld))/vpnorm
			else
               yldot(iv2) = nurlxu*
     .			 (up(ixd,iy,ifld) - up(ix,iy,ifld))/vpnorm
			endif
            if(isfix(1).eq.2) yldot(iv2) = nurlxu *
     .                           (0. - up(ix,iy,ifld))/vpnorm

            if(isfix(1).eq.2 .and. yyb(iy,1).gt.rlimiter) then
                cs = sqrt( (te(ix,iy)+ti(ix,iy))/mi(ifld) )
                	yldot(iv2) = nurlxu*
     .                          (cssign*cs -up(ix,iy,ifld))/vpnorm

               endif
            endif
         enddo



		if(isteonxy(ix,iy) .eq. 1) then
         	iv1 = idxte(ix,iy)
           	yldot(iv1) = nurlxe * ne(ix,iy) *
     .                     (teb*ev*tebprof(iy) - te(ix,iy))/ennorm
           	if(isfix(1).eq.2) yldot(iv1) = nurlxe * ne(ix,iy) *
     .                               (te(ixd,iy) - te(ix,iy))/ennorm
           	if(isfix(1).eq.2 .and. yyb(iy,1).gt.rlimiter) then
            	yldot(iv1) = - nurlxe*(feex(ix2,iy)/sx(ix2,iy) - bcee*
     .                     ne(ix,iy)*vex(ix2,iy)*te(ix,iy))/(vpnorm*ennorm)
         	endif
        endif


        if(istionxy(ix,iy) .eq. 1) then
        	iv2 = idxti(ix,iy)
            yldot(iv2) = nurlxi * ne(ix,iy) *
     .                     (tib*ev*tibprof(iy) - ti(ix,iy))/ennorm

            if(isfix(1).eq.2) yldot(iv2) = nurlxi * ne(ix,iy) *
     .                     (ti(ixd,iy) - ti(ix,iy))/ennorm
        	if(isfix(1).eq.2 .and. yyb(iy,1).gt.rlimiter) then
            	yldot(iv2) = -nurlxi*
     .          ( feix(ix2,iy) - bcei*ti(ix,iy)*fac2sp*fnix(ix2,iy,1) ) / 
     .                                       (vpnorm*ennorm*sx(ix2,iy))
            endif
         endif




        do igsp = 1, ngsp
        	if(isngonxy(ix,iy,igsp) .eq. 1) then
            	iv = idxg(ix,iy,igsp)
               	yldot(iv) = nurlxg * (ngbackg(igsp) - 
     .                                         ng(ix,iy,igsp))/n0g(igsp)
               	if(isfix(1).eq.2) yldot(iv) = nurlxg * 
     .                        (ng(ixd,iy,igsp) - ng(ix,iy,igsp))/n0g(igsp)
               	if(isfix(1).eq.2 .and. yyb(iy,1).gt.rlimiter) then
                  	flux_inc = fac2sp*fnix(ix2,iy,1)
                  	if (ishymol.eq.1 .and. igsp.eq.2) then
c 					Previously, lb used alb 3rd index 1, whereas 
c 					lb used nxpt. Now uses 1 for both. Min Tg was
c 					previously set by temin, now tgmin				                  	
                      flxa= ismolcrm*(1-alb(iy,1,nxpt))
     .               	* onesided_maxwellian(
     .                  	engbsr*tg(ixd,iy,1), ng(ixd,iy,1),
     .                  	mg(1), sx(ix2, iy), engbsr*tgmin*ev
     .              	)

                      if (isupgon(1) .eq. 1) then  # two atoms for one molecule
                      	flux_inc = 0.5*( fnix(ix2,iy,1) + fnix(ix2,iy,2) -cssign*flxa) 
                      else
                      	flux_inc = 0.5*( fnix(ix2,iy,1) + fngx(ix2,iy,1) -cssign*flxa) 
                      endif                  	
 					endif
c 					Previously lb had now engsbr-factor, now included in 
c 					both. Previously, lb used igsp index 1 whereas rb
c 					used igsp: now both use igsp
                  	osmw = onesided_maxwellian(
     .                  engbsr*tg(ixd,iy,1), ng(ixd,iy,igsp), mg(igsp),
     .                  isoldalbarea*sx(ix2,iy) + (1-isoldalbarea)*sxnp(ix2,iy),
     .                  engbsr*tgmin*ev
     .            	)       					
                  	yldot(iv) = -nurlxg * ( fngx(ix2,iy,igsp)
     .                      + cssign*fngx_use(iy,igsp,1)
     .                      - cssign*fngxs(iy,igsp,1) 
     .                      + recyb(iy,igsp,1)*flux_inc
     .                      - cssign*(1-alb(iy,igsp,1))*osmw
     .              	) / (vpnorm*n0g(igsp)*sx(ix2,iy))
            endif


            if (
     .			(isleftbound)
     .			.and.
     .			(is1D_gbx.eq.1)
     .     	) yldot(iv) = nurlxg*(ng(ixd,iy,igsp) -
     .                                    ng(ix,iy,igsp))/n0g(igsp)


           	if (istgonxy(ix,iy,igsp) == 1) then
            	iv = idxtg(ix,iy,igsp)
             	yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(ix,iy,igsp))/(temp0*ev)
             	if(isfix(1)==2) then #just above applies if isfixlb=1
               		yldot(iv)=nurlxg*
     .					(tg(ixd,iy,igsp)-tg(ix,iy,igsp))/(temp0*ev)
             	endif
             	if(isfix(1)==2 .and. yyb(iy,1) > rlimiter) then
               		yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .					tg(ix,iy,igsp))/(temp0*ev)
               	endif
            endif
        enddo


        if (isphionxy(ix,iy) .eq. 1) then
        	iv = idxphi(ix,iy)
            yldot(iv) = nurlxp*(phi(ixd,iy) - phi(ix,iy))/temp0
        endif
    enddo

c 		COMMON        
c 		LEFT
		ixpt = ixpt2
		iysptrx = iysptrx2
c 		RIGHT
		ixpt = ixpt1
		iysptrx = iysptrx1

 

c...  If isfixlb=2, check if i2,i5 range for yldot loop in pandf includes ixpt(1)
c...  Then overwrite pandf value is up(ixpt(1)) --> 0 eqn.
      if (isfix(1) .eq. 2) then
         if (i2.le.ixpt(1) .and. i5.ge.ixpt(1) .and. j2.le.iysptrx(1)) then  
            do ifld = 1, nusp
            	do iy = 0+1-iymnbcl, iysptrx(1)
					if(isuponxy(ixpt(1),iy,ifld)==1) then
						iv = idxu(ixpt(1),iy,ifld)
						if isleftbound then
							yldot(iv) = nurlxu*(
	 .							0.-up(ixpt(1),iy,ifld))/vpnorm	
						else
	 .						if (isupcore(ifld).eq.0)
	 .							yldot(iv) = nurlxu*(
	 .								0.-up(ixpt(1),iy,ifld))/vpnorm



               
                 if(isuponxy(ixpt(1),iy,ifld)==1) then
                     iv = idxu(ixpt(1),iy,ifld)
                     yldot(iv) = nurlxu*(0.-up(ixpt(1),iy,ifld))/vpnorm
                 endif
               enddo
            enddo
         endif 
      endif            # end of isfixlb=2, check for ix=ixpt2



             if (isupcore(ifld).eq.0) then
               do iy = 0+1-iymnbcl, iysptrx(1)
                 if(isuponxy(ixpt(1),iy,ifld)==1) then
                   iv = idxu(ixpt(1),iy,ifld)
                   yldot(iv) = nurlxu*(0.-up(ixpt(1),iy,ifld))/vpnorm
                 endif
               enddo
             endif
           enddo
         endif 
      endif                             #end symm. BC at ix=nx+1 for isfixrb=2












    END SUBROUTINE symmetric_boundary

	SUBROUTINE left_symmetric_boundary(neq, yl, yldot)
	IMPLICIT NONE
    INTEGER, INTENT(IN) :: neq
    REAL, INTENT(OUT) :: yl(neq), yldot(neq)   


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


      END SUBROUTINE left_ymmetric_boundary



	SUBROUTINE right_symmetric_boundary(neq, yl, yldot)
	IMPLICIT NONE
    INTEGER, INTENT(IN) :: neq
    REAL, INTENT(OUT) :: yl(neq), yldot(neq)   


c********************************************************************
c...  First, check if isfixrb=2 for using symmetry BC at ix = nx+1
c*******************************************************************
      if (i6 .ge. nx+1 .and. isfixrb(1) .gt. 0) then  #begin symmetry BC at nx+1
       do iy = j2, j5
c...  First do the ion density
         do ifld = 1, nisp
            if(isnionxy(ix,iy,ifld) .eq. 1) then
               iv1 = idxn(ix,iy,ifld)
               yldot(iv1) = nurlxn *
     .              (nib(ifld)*nibprof(iy)-ni(ix,iy,ifld))/n0(ifld)
               if(isfixrb(1).eq.2) yldot(iv1) = nurlxn * (1/n0(ifld)) *
     .              (ni(ixd,iy,ifld) - ni(ix,iy,ifld))
            endif
         enddo

c...  Now do the parallel velocity
         do ifld = 1, nusp
            if(isuponxy(ix,iy,ifld) .eq. 1) then
               iv2 = idxu(ix,iy,ifld)
               yldot(iv2) = nurlxu*(up(ixd,iy,ifld) - up(ix,iy,ifld))/
     .                                                         vpnorm
               if(isfixrb(1).eq.2) yldot(iv2) = nurlxu*
     .                            (0. - up(ix,iy,ifld))/vpnorm
               if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
                  cs = sqrt( (te(ix,iy)+ti(ix,iy))/mi(ifld) )
                  yldot(iv2) = nurlxu*(cs -up(ix,iy,ifld))/vpnorm
               endif
            endif
         enddo

c...  now do the gas and temperatures
         if(isteonxy(ix,iy) .eq. 1) then
           iv1 = idxte(ix,iy)
           yldot(iv1) = nurlxe * ne(ix,iy) *
     .                     (teb*ev*tebprof(iy) - te(ix,iy))/ennorm
           if(isfixrb(1).eq.2) yldot(iv1) = nurlxe * ne(ix,iy) *
     .                               (te(ixd,iy) - te(ix,iy))/ennorm
           if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
              yldot(iv1) = - nurlxe*(feex(ixd,iy)/sx(ixd,iy) - bcee*
     .                           ne(ix,iy)*vex(ixd,iy)*te(ix,iy))/
     .                                 (vpnorm*ennorm)
           endif
         endif
         if(istionxy(ix,iy) .eq. 1) then
            iv2 = idxti(ix,iy)
            yldot(iv2) = nurlxi * ne(ix,iy) *
     .                     (tib*ev*tibprof(iy) - ti(ix,iy))/ennorm
            if(isfixrb(1).eq.2) yldot(iv2) = nurlxi * ne(ix,iy) *
     .                               (ti(ixd,iy) - ti(ix,iy))/ennorm
            if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
               yldot(iv2) = -nurlxi*
     .          ( feix(ixd,iy) - bcei*ti(ix,iy)*fac2sp*fnix(ixd,iy,1) ) /
     .                                         (vpnorm*ennorm*sx(ixd,iy))
            endif
         endif
         do igsp = 1, nhgsp # not valid for ngsp > nhgsp; only on hydrog. gas
            if (isngonxy(ix,iy,igsp) .eq. 1) then
               iv = idxg(ix,iy,igsp)
               yldot(iv) = nurlxg * (ngbackg(igsp) - 
     .                                     ng(ix,iy,igsp)) / n0g(igsp)
               if(isfixrb(1).eq.2) yldot(iv) = nurlxg * 
     .                     (ng(ixd,iy,igsp) - ng(ix,iy,igsp))/n0g(igsp)
               if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
                  flux_inc = fac2sp*fnix(ixd,iy,1)
                  if (ishymol.eq.1 .and. igsp.eq.2) then
                    flxa= ismolcrm*(1-albrb(iy,1,nxpt))* onesided_maxwellian(
     .                  engbsr*tg(ixd,iy,1), ng(ixd,iy,1), mg(1), 
     .                  sx(ixd, iy), engbsr*tgmin*ev
     .              )
                    if (isupgon(1) .eq. 1) then  # two atoms for one molecule
                      flux_inc = 0.5*( fnix(ixd,iy,1) + fnix(ixd,iy,2) -flxa) 
                    else
                      flux_inc = 0.5*( fnix(ixd,iy,1) + fngx(ixd,iy,1) -flxa) 
                    endif
                  endif
                  osmw = onesided_maxwellian(
     .                  engbsr*tg(ixd,iy,1), ng(ixd,iy,igsp), mg(igsp),
     .                  isoldalbarea*sx(ixd,iy) + (1-isoldalbarea)*sxnp(ixd,iy),
     .                  engbsr*tgmin*ev
     .            )      
                  yldot(iv) = -nurlxg * ( fngx(ixd,iy,igsp)
     .                      + fngxrb_use(iy,igsp,1)
     .                      - fngxsrb(iy,igsp,1) 
     .                      + recyrb(iy,igsp,1)*flux_inc
     .                      - (1-albrb(iy,igsp,nxpt))*osmw
     .                  ) / (vpnorm*n0g(igsp)*sx(ixd,iy))
               endif
            endif
         enddo      # end of igsp loop over gas


c ... Neutral temperature - test if tg eqn is on, then set BC
	 do igsp = 1, ngsp
           if (istgonxy(ix,iy,igsp) == 1) then
             iv = idxtg(ix,iy,igsp)
             yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(ix,iy,igsp))/(temp0*ev)
             if(isfixrb(1)==2) then #just above applies if isfixrb=1
               yldot(iv)=nurlxg*(tg(ixd,iy,igsp)-tg(ix,iy,igsp))/
     .                                                      (temp0*ev)
             endif
             if(isfixrb(1)==2 .and. yyrb(iy,1) > rlimiter) then
               yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(ix,iy,igsp))/(temp0*ev)
             endif
           endif
         enddo

         if (isphionxy(ix,iy) .eq. 1) then
            iv = idxphi(ix,iy)
            yldot(iv) = nurlxp*(phi(ixd,iy) - phi(ix,iy))/temp0
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




    END SUBROUTINE right_symmetric_boundary