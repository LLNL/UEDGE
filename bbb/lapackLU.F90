
subroutine CopyMatrixBanded(wpcopy,wp,lowd,neq)
implicit none
double precision wp(1:lowd,1:neq)
integer lowd,neq
double precision wpcopy(1:lowd,1:neq)
wpcopy(1:lowd,1:neq)=wp(1:lowd,1:neq)
end subroutine CopyMatrixBanded

subroutine CompareMatrixLU(wpcopy,wp,lowd,neq,lbw,ubw)
implicit none
double precision wp(1:lowd,1:neq)
integer lowd,neq,i,j,lbw,ubw,upsilon
double precision wpcopy(1:lowd,1:neq)

do i=1,lowd
! Multipiers stored in output matrix from dgbco_u and have opposite signs between
if (i>lbw+ubw+1) then
upsilon=-1
else
upsilon=1
endif
do j=1,neq
if (abs(wpcopy(i,j)-wp(i,j)).gt.1e-14) then
                    if (max(abs(wp(i,j)),abs(wpcopy(i,j)))>0) then
                    if (abs(wpcopy(i,j)-upsilon*wp(i,j))/max(abs(wp(i,j)),abs(wpcopy(i,j)))>1e-14) then
                    if (i<=lbw+ubw+1) then
                    write(*,*) '>>>>',i,j,wpcopy(i,j),wp(i,j),abs(wpcopy(i,j)-upsilon*wp(i,j))/max(abs(wpcopy(i,j)),abs(wp(i,j)))
                    call xerrab('Mismatch in LU decomposition between dgbco_u and ')
                    endif
                    endif
                    else
                    write(*,*) '>>>>',i,j
                    endif
                    !call xerrab('diff in rhsnk')
                endif
enddo
enddo

end subroutine CompareMatrixLU

subroutine lapacksgbco(wp,lowd,neq,lbw,ubw,ipvt,rcond,rwk1)
      use LapackLU,only:CheckBandedPrecond,BandedPrecondVerbose
      integer lowd,neq,lbw,ubw,ipvt(neq)
      double precision wp(lowd,neq),rwk1(neq),wpcopy(lowd,neq)
      double precision rcond
      integer ierr,t_start
      real,external::tock


      if (CheckBandedPrecond.gt.0) then
            wpcopy(1:lowd,1:neq)=wp(1:lowd,1:neq)
      endif

      if (BandedPrecondVerbose.gt.0) call tick(t_start)
      call DGBTRF( neq, neq, lbw, ubw, wp, lowd, ipvt, ierr )
      if (BandedPrecondVerbose.gt.0)   write(*,*) 'dgbtrf time:',tock(t_start),'ierr',ierr

      if (CheckBandedPrecond.gt.0) then
         call tick(t_start)
         call dgbco_u (wpcopy, lowd, neq, lbw, ubw, ipvt, rcond, rwk1)
         write(*,*) 'sbgco time:',tock(t_start)
         write(*,*) 'rcond=',rcond
         write(*,*) 'lbw','ubw','lowd',lbw,ubw,lowd
         call CompareMatrixLU(wpcopy,wp,lowd,neq,lbw,ubw)
      endif
    rcond=0
end subroutine lapacksgbco
