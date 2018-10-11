# newuoa.m,v 1.6 2006/08/12 18:27:06 bulmer Exp

      subroutine xnewuoa(n,npt,x,rhobeg,rhoend,maxfun,iprint)

      implicit real (a-h,o-z), integer (i-n)
      Use(UOA)
      dimension x(*)

    # This subroutine seeks the least value of a function of many variables,
    # by a trust region method that forms quadratic models by interpolation.
    # There can be some freedom in the interpolation conditions, which is
    # taken up by minimizing the Frobenius norm of the change to the second
    # derivative of the quadratic model, beginning with a zero matrix. The
    # arguments of the subroutine are as follows.
    #
    # N must be set to the number of variables and must be at least two.
    # NPT is the number of interpolation conditions. Its value must be in the
    #   interval [N+2,(N+1)(N+2)/2].
    # Initial values of the variables must be set in X(1),X(2),...,X(N). They
    #   will be changed to the values that give the least calculated F.
    # RHOBEG and RHOEND must be set to the initial and final values of a trust
    #   region radius, so both must be positive with RHOEND<=RHOBEG. Typically
    #   RHOBEG should be about one tenth of the greatest expected change to a
    #   variable, and RHOEND should indicate the accuracy that is required in
    #   the final values of the variables.
    # MAXFUN must be set to an upper bound on the number of calls of CALFUN.
    # The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
    #   amount of printing. Specifically, there is no output if IPRINT=0 and
    #   there is output only at the return if IPRINT=1. Otherwise, each new
    #   value of RHO is printed, with the best vector of variables so far and
    #   the corresponding value of the objective function. Further, each new
    #   value of F with its variables are output if IPRINT=3.
    #
    # SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
    # the value of the objective function for the variables X(1),X(2),...,X(N).
    #
    # Partition the working space array, so that different parts of it can be
    # treated separately by the subroutine that performs the main calculation.

      np=n+1
      nptm=npt-np

    # Allocate workspace
      n_uoa=n
      npt_uoa=npt
      m_uoa=(npt+13)*(npt+n)+3*n*(n+3)/2
      call gchange("UOA",0)

    # Provide a partition of W for subroutine NEWUOB.
    # The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
    # W plus the space that is needed by the last array of NEWUOB.
      ndim=npt+n
      ixb=1
      ixo=ixb+n
      ixn=ixo+n
      ixp=ixn+n
      ifv=ixp+n*npt
      igq=ifv+npt
      ihq=igq+n
      ipq=ihq+(n*np)/2
      ibmat=ipq+npt
      izmat=ibmat+ndim*n
      id=izmat+npt*nptm
      ivl=id+n
      iw=ivl+ndim

      call newuob(n, npt, x, rhobeg, rhoend, iprint, maxfun, 
     .            w_uoa(ixb), w_uoa(ixo), w_uoa(ixn), w_uoa(ixp), 
     .            w_uoa(ifv), w_uoa(igq), w_uoa(ihq), w_uoa(ipq), 
     .            w_uoa(ibmat), w_uoa(izmat), ndim, w_uoa(id), 
     .            w_uoa(ivl), w_uoa(iw))

      return

      end # newuoa

c----------------------------------------------------------------------
# newuob.m,v 1.4 2006/08/12 18:27:06 bulmer Exp

      subroutine newuob(n,npt,x,rhobeg,rhoend,iprint,maxfun,xbase, 
     .             xopt,xnew,xpt,fval,gq,hq,pq,bmat,zmat,ndim,d,vlag,w)

      implicit real (a-h,o-z), integer (i-n)
      Use(UOA)
      dimension x(*),xbase(*),xopt(*),xnew(*),xpt(npt,*),fval(*)
      dimension gq(*),hq(*),pq(*),bmat(ndim,*),zmat(npt,*),d(*),vlag(*),w(*)

    # The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
    #   to the corresponding arguments in SUBROUTINE NEWUOA.
    # XBASE will hold a shift of origin that should reduce the contributions
    #   from rounding errors to values of the model and Lagrange functions.
    # XOPT will be set to the displacement from XBASE of the vector of
    #   variables that provides the least calculated F so far.
    # XNEW will be set to the displacement from XBASE of the vector of
    #   variables for the current calculation of F.
    # XPT will contain the interpolation point coordinates relative to XBASE.
    # FVAL will hold the values of F at the interpolation points.
    # GQ will hold the gradient of the quadratic model at XBASE.
    # HQ will hold the explicit second derivatives of the quadratic model.
    # PQ will contain the parameters of the implicit second derivatives of
    #   the quadratic model.
    # BMAT will hold the last N columns of H.
    # ZMAT will hold the factorization of the leading NPT by NPT submatrix of
    #   H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
    #   the elements of DZ are plus or minus one, as specified by IDZ.
    # NDIM is the first dimension of BMAT and has the value NPT+N.
    # D is reserved for trial steps from XOPT.
    # VLAG will contain the values of the Lagrange functions at a new point X.
    #   They are part of a product that requires VLAG to be of length NDIM.
    # The array W will be used for working space. Its length must be at least
    #   10*NDIM = 10*(NPT+N).

      half=0.5
      one=1.0
      tenth=0.1
      zero=0.0
      np=n+1
      nh=(n*np)/2
      nptm=npt-np
      nftest=max0(maxfun,1)

    # Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.

      do  j=1,n
            xbase(j)=x(j)
          do k=1,npt
              xpt(k,j)=zero
          enddo
          do i=1,ndim
              bmat(i,j)=zero
          enddo
      enddo
      do ih=1,nh
          hq(ih)=zero
      enddo
      do k=1,npt
          pq(k)=zero
          do j=1,nptm
              zmat(k,j)=zero
          enddo
      enddo

    # Begin the initialization procedure. NF becomes one more than the number
    # of function values so far. The coordinates of the displacement of the
    # next initial interpolation point from XBASE are set in XPT(NF,.).

      rhosq=rhobeg*rhobeg
      recip=one/rhosq
      reciq=sqrt(half)/rhosq
      nf=0
 50   nfm=nf
      nfmm=nf-n
      nf=nf+1
      if (nfm <= 2*n) then
        if (nfm >= 1 .and. nfm <= n) then
            xpt(nf,nfm)=rhobeg
        elseif (nfm > n) then
            xpt(nf,nfmm)=-rhobeg
        endif
      else
        itemp=(nfmm-1)/n
        jpt=nfm-itemp*n-n
        ipt=jpt+itemp
        if (ipt > n) then
            itemp=jpt
            jpt=ipt-n
            ipt=itemp
        endif
        xipt=rhobeg
        if (fval(ipt+np) < fval(ipt+1)) xipt=-xipt
        xjpt=rhobeg
        if (fval(jpt+np) < fval(jpt+1)) xjpt=-xjpt
        xpt(nf,ipt)=xipt
        xpt(nf,jpt)=xjpt
      endif

    # Calculate the next value of F, label 70 being reached immediately
    # after this calculation. The least function value so far and its index
    # are required.

      do j=1,n
        x(j)=xpt(nf,j)+xbase(j)
      enddo
      goto 310
 70   fval(nf)=f
      if (nf == 1) then
        fbeg=f
        fopt=f
        kopt=1
      elseif (f < fopt) then
        fopt=f
        kopt=nf
      endif

    # Set the nonzero initial elements of BMAT and the quadratic model in
    # the cases when NF is at most 2*N+1.

      if (nfm <= 2*n) then
        if (nfm >= 1 .and. nfm <= n) then
            gq(nfm)=(f-fbeg)/rhobeg
            if (npt < nf+n) then
                bmat(1,nfm)=-one/rhobeg
                bmat(nf,nfm)=one/rhobeg
                bmat(npt+nfm,nfm)=-half*rhosq
            endif
          elseif (nfm > n) then
            bmat(nf-n,nfmm)=half/rhobeg
            bmat(nf,nfmm)=-half/rhobeg
            zmat(1,nfmm)=-reciq-reciq
            zmat(nf-n,nfmm)=reciq
            zmat(nf,nfmm)=reciq
            ih=(nfmm*(nfmm+1))/2
            temp=(fbeg-f)/rhobeg
            hq(ih)=(gq(nfmm)-temp)/rhobeg
            gq(nfmm)=half*(gq(nfmm)+temp)
        endif

      else
        # Set the off-diagonal second derivatives of the Lagrange functions and
        # the initial quadratic model.
            ih=(ipt*(ipt-1))/2+jpt
          if (xipt < zero) ipt=ipt+n
          if (xjpt < zero) jpt=jpt+n
          zmat(1,nfmm)=recip
          zmat(nf,nfmm)=recip
            zmat(ipt+1,nfmm)=-recip
          zmat(jpt+1,nfmm)=-recip
          hq(ih)=(fbeg-fval(ipt+1)-fval(jpt+1)+f)/(xipt*xjpt)
      endif
      if (nf < npt) goto 50

    # Begin the iterative procedure, because the initial model is complete.

      rho=rhobeg
      rho_uoa=rho
      delta=rho
      idz=1
      diffa=zero
      diffb=zero
      itest=0
      xoptsq=zero
      do i=1,n
        xopt(i)=xpt(kopt,i)
        xoptsq=xoptsq+xopt(i)**2
      enddo
 90   nfsav=nf

    # Generate the next trust region step and test its length. Set KNEW
    # to -1 if the purpose of the next F will be to improve the model.

100   knew=0
      call trsapp(n,npt,xopt,xpt,gq,hq,pq,delta,d,w,w(np), 
     .            w(np+n),w(np+2*n),crvmin)
      dsq=zero
      do i=1,n
        dsq=dsq+d(i)**2
      enddo
      dnorm=min(delta,sqrt(dsq))
      if (dnorm < half*rho) then
        knew=-1
        delta=tenth*delta
        ratio=-1.0
        if (delta <= 1.5*rho) delta=rho
        if (nf <= nfsav+2) goto 460
        temp=0.125*crvmin*rho*rho
        if (temp <= max(diffa,diffb,diffc)) goto 460
        goto 490
      endif

    # Shift XBASE if XOPT may be too far from XBASE. First make the changes
    # to BMAT that do not depend on ZMAT.

120   if (dsq <= 1.0e-3*xoptsq) then
        tempq=0.25*xoptsq
        do k=1,npt
            sum=zero
            do i=1,n
                sum=sum+xpt(k,i)*xopt(i)
            enddo
            temp=pq(k)*sum
            sum=sum-half*xoptsq
            w(npt+k)=sum
            do i=1,n
                gq(i)=gq(i)+temp*xpt(k,i)
                xpt(k,i)=xpt(k,i)-half*xopt(i)
                vlag(i)=bmat(k,i)
                w(i)=sum*xpt(k,i)+tempq*xopt(i)
                ip=npt+i
                do j=1,i
                    bmat(ip,j)=bmat(ip,j)+vlag(i)*w(j)+w(i)*vlag(j)
                enddo
            enddo
        enddo

        # Then the revisions of BMAT that depend on ZMAT are calculated.

        do k=1,nptm
            sumz=zero
            do i=1,npt
                sumz=sumz+zmat(i,k)
                w(i)=w(npt+i)*zmat(i,k)
            enddo
            do j=1,n
                sum=tempq*sumz*xopt(j)
                do i=1,npt
                    sum=sum+w(i)*xpt(i,j)
                enddo
                vlag(j)=sum
                if (k < idz) sum=-sum
                do i=1,npt
                    bmat(i,j)=bmat(i,j)+sum*zmat(i,k)
                enddo
            enddo
            do i=1,n
                ip=i+npt
                temp=vlag(i)
                if (k < idz) temp=-temp
                do j=1,i
                    bmat(ip,j)=bmat(ip,j)+temp*vlag(j)
                enddo
            enddo
        enddo

        # The following instructions complete the shift of XBASE, including
        # the changes to the parameters of the quadratic model.

          ih=0
          do j=1,n
            w(j)=zero
            do k=1,npt
                w(j)=w(j)+pq(k)*xpt(k,j)
                xpt(k,j)=xpt(k,j)-half*xopt(j)
            enddo
            do i=1,j
                ih=ih+1
                if (i < j) gq(j)=gq(j)+hq(ih)*xopt(i)
                gq(i)=gq(i)+hq(ih)*xopt(j)
                hq(ih)=hq(ih)+w(i)*xopt(j)+xopt(i)*w(j)
                bmat(npt+i,j)=bmat(npt+j,i)
            enddo
          enddo
        do j=1,n
            xbase(j)=xbase(j)+xopt(j)
            xopt(j)=zero
        enddo
        xoptsq=zero
      endif

    # Pick the model step if KNEW is positive. A different choice of D
    # may be made later, if the choice of D by BIGLAG causes substantial
    # cancellation in DENOM.

      if (knew > 0) then
          call biglag(n,npt,xopt,xpt,bmat,zmat,idz,ndim,knew,dstep, 
     .                d,alpha,vlag,vlag(npt+1),w,w(np),w(np+n))
      endif

    # Calculate VLAG and BETA for the current choice of D. The first NPT
    # components of W_check will be held in W.

      do k=1,npt
        suma=zero
        sumb=zero
        sum=zero
        do j=1,n
            suma=suma+xpt(k,j)*d(j)
            sumb=sumb+xpt(k,j)*xopt(j)
            sum=sum+bmat(k,j)*d(j)
        enddo
        w(k)=suma*(half*suma+sumb)
        vlag(k)=sum
      enddo
      beta=zero
      do k=1,nptm
        sum=zero
        do i=1,npt
            sum=sum+zmat(i,k)*w(i)
        enddo
        if (k < idz) then
            beta=beta+sum*sum
            sum=-sum
        else
            beta=beta-sum*sum
        endif
        do i=1,npt
            vlag(i)=vlag(i)+sum*zmat(i,k)
        enddo
      enddo
      bsum=zero
      dx=zero
      do j=1,n
        sum=zero
        do i=1,npt
            sum=sum+w(i)*bmat(i,j)
        enddo
        bsum=bsum+sum*d(j)
        jp=npt+j
        do k=1,n
            sum=sum+bmat(jp,k)*d(k)
        enddo
        vlag(jp)=sum
        bsum=bsum+sum*d(j)
        dx=dx+d(j)*xopt(j)
      enddo
      beta=dx*dx+dsq*(xoptsq+dx+dx+half*dsq)+beta-bsum
      vlag(kopt)=vlag(kopt)+one

    # If KNEW is positive and if the cancellation in DENOM is unacceptable,
    # then BIGDEN calculates an alternative model step, XNEW being used for
    # working space.

      if (knew > 0) then
        temp=one+alpha*beta/vlag(knew)**2
        if (abs(temp) <= 0.8) then
            call bigden(n,npt,xopt,xpt,bmat,zmat,idz,ndim,kopt, 
     .                  knew,d,w,vlag,beta,xnew,w(ndim+1),w(6*ndim+1))
        endif
      endif

    # Calculate the next value of the objective function.

290   do i=1,n
        xnew(i)=xopt(i)+d(i)
        x(i)=xbase(i)+xnew(i)
      enddo
      nf=nf+1
310   if (nf > nftest) then
        nf=nf-1
        if (iprint > 0) then
            call remark("newuoa: CALFUN has been called MAXFUN times.")
        endif
        goto 530
      endif
      call calfun(n,x,f)
      if (iprint == 3) then
        print 330, nf,f,(x(i),i=1,n)
330     format (/4x,'Function number',i6,'    F =',1pd18.10, 
     .               '    The corresponding X is:'/(2x,5d15.6))
      endif
      if (nf <= npt) goto 70
      if (knew == -1) goto 530

    # Use the quadratic model to predict the change in F due to the step D,
    # and set DIFF to the error of this prediction.

      vquad=zero
      ih=0
      do j=1,n
        vquad=vquad+d(j)*gq(j)
        do i=1,j
            ih=ih+1
            temp=d(i)*xnew(j)+d(j)*xopt(i)
            if (i == j) temp=half*temp
            vquad=vquad+temp*hq(ih)
        enddo
      enddo
      do k=1,npt
        vquad=vquad+pq(k)*w(k)
      enddo
      diff=f-fopt-vquad
      diffc=diffb
      diffb=diffa
      diffa=abs(diff)
      if (dnorm > rho) nfsav=nf

    # Update FOPT and XOPT if the new F is the least value of the objective
    # function so far. The branch when KNEW is positive occurs if D is not
    # a trust region step.

      fsave=fopt
      if (f < fopt) then
        fopt=f
        xoptsq=zero
        do i=1,n
            xopt(i)=xnew(i)
            xoptsq=xoptsq+xopt(i)**2
        enddo
      endif
      ksave=knew
      if (knew > 0) goto 410

    # Pick the next value of DELTA after a trust region step.

      if (vquad >= zero) then
        if (iprint > 0) then
            call remark("newuoa: a trust region step has failed to reduce Q.")
        endif
        goto 530
      endif
      ratio=(f-fsave)/vquad
      if (ratio <= tenth) then
        delta=half*dnorm
      elseif (ratio <= 0.7) then
        delta=max(half*delta,dnorm)
      else
        delta=max(half*delta,dnorm+dnorm)
      endif
      if (delta <= 1.5*rho) delta=rho

    # Set KNEW to the index of the next interpolation point to be deleted.

      rhosq=max(tenth*delta,rho)**2
      ktemp=0
      detrat=zero
      if (f >= fsave) then
        ktemp=kopt
        detrat=one
      endif
      do k=1,npt
        hdiag=zero
        do j=1,nptm
            temp=one
            if (j < idz) temp=-one
            hdiag=hdiag+temp*zmat(k,j)**2
        enddo
        temp=abs(beta*hdiag+vlag(k)**2)
        distsq=zero
        do j=1,n
            distsq=distsq+(xpt(k,j)-xopt(j))**2
        enddo
        if (distsq > rhosq) temp=temp*(distsq/rhosq)**3
        if (temp > detrat .and. k .ne. ktemp) then
            detrat=temp
            knew=k
        endif
      enddo
      if (knew == 0) goto 460

    # Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
    # can be moved. Begin the updating of the quadratic model, starting
    # with the explicit second derivative term.

410   call update(n,npt,bmat,zmat,idz,ndim,vlag,beta,knew,w)
      fval(knew)=f
      ih=0
      do i=1,n
          temp=pq(knew)*xpt(knew,i)
        do j=1,i
            ih=ih+1
            hq(ih)=hq(ih)+temp*xpt(knew,j)
        enddo
      enddo
      pq(knew)=zero

    # Update the other second derivative parameters, and then the gradient
    # vector of the model. Also include the new interpolation point.

      do j=1,nptm
        temp=diff*zmat(knew,j)
        if (j < idz) temp=-temp
        do k=1,npt
            pq(k)=pq(k)+temp*zmat(k,j)
        enddo
      enddo
      gqsq=zero
      do i=1,n
        gq(i)=gq(i)+diff*bmat(knew,i)
        gqsq=gqsq+gq(i)**2
        xpt(knew,i)=xnew(i)
      enddo

    # If a trust region step makes a small change to the objective function,
    # then calculate the gradient of the least Frobenius norm interpolant at
    # XBASE, and store it in W, using VLAG for a vector of right hand sides.

      if (ksave == 0 .and. delta == rho) then
        if (abs(ratio) > 1.0e-2) then
            itest=0
        else
            do k=1,npt
                vlag(k)=fval(k)-fval(kopt)
            enddo
            gisq=zero
            do i=1,n
                sum=zero
                do k=1,npt
                    sum=sum+bmat(k,i)*vlag(k)
                enddo
                gisq=gisq+sum*sum
                w(i)=sum
            enddo

            # Test whether to replace the new quadratic model by the
            # least Frobenius norm interpolant, making the replacement
            # if the test is satisfied.

            itest=itest+1
            if (gqsq < 1.0e+2*gisq) itest=0
            if (itest >= 3) then
                do i=1,n
                    gq(i)=w(i)
                enddo
                do ih=1,nh
                    hq(ih)=zero
                enddo
                do j=1,nptm
                    w(j)=zero
                    do k=1,npt
                        w(j)=w(j)+vlag(k)*zmat(k,j)
                    enddo
                    if (j < idz) w(j)=-w(j)
                enddo
                do k=1,npt
                    pq(k)=zero
                    do j=1,nptm
                        pq(k)=pq(k)+zmat(k,j)*w(j)
                    enddo
                enddo
                itest=0
            endif
        endif
      endif
      if (f < fsave) kopt=knew

    # If a trust region step has provided a sufficient decrease in F, then
    # branch for another trust region calculation. The case KSAVE>0 occurs
    # when the new function value was calculated by a model step.

      if (f <= fsave+tenth*vquad) goto 100
      if (ksave > 0) goto 100

    # Alternatively, find out if the interpolation points are close enough
    # to the best point so far.

      knew=0
460   distsq=4.0*delta*delta
      do k=1,npt
        sum=zero
        do j=1,n
            sum=sum+(xpt(k,j)-xopt(j))**2
        enddo
        if (sum > distsq) then
            knew=k
            distsq=sum
        endif
      enddo

    # If KNEW is positive, then set DSTEP, and branch back for the next
    # iteration, which will generate a "model step".

      if (knew > 0) then
        dstep=max(min(tenth*sqrt(distsq),half*delta),rho)
        dsq=dstep*dstep
        goto 120
      endif
      if (ratio > zero) goto 100
      if (max(delta,dnorm) > rho) goto 100

    # The calculations with the current value of RHO are complete. Pick the
    # next values of RHO and DELTA.

490   if (rho > rhoend) then
        delta=half*rho
        ratio=rho/rhoend
        if (ratio <= 16.0) then
            rho=rhoend
        elseif (ratio <= 250.0) then
            rho=sqrt(ratio)*rhoend
        else
            rho=tenth*rho
        endif
        rho_uoa=rho
        delta=max(delta,rho)
        if (iprint >= 2) then
            if (iprint >= 3) print 500
500         format (5x)
            print 510, rho,nf
510         format (/4x,'New RHO =',1pd11.4,5x,'Number of', 
     .                  ' function values =',i6)
            print 520, fopt,(xbase(i)+xopt(i),i=1,n)
520         format (4x,'Least value of F =',1pd23.15,9x, 
     .                 'The corresponding X is:'/(2x,5d15.6))
        endif
        goto 90
      endif

    # Return from the calculation, after another Newton-Raphson step, if
    # it is too short to have been tried before.

      if (knew == -1) goto 290
530   if (fopt <= f) then
        do i=1,n
            x(i)=xbase(i)+xopt(i)
        enddo
        f=fopt
      endif
      if (iprint >= 1) then
        print 550, nf
550     format (/4x,'At the return from NEWUOA',5x, 
     .               'Number of function values =',i6)
        print 520, f,(x(i),i=1,n)
      endif
      return

      end # newuob

c----------------------------------------------------------------------
# bigden.m,v 1.3 2006/08/12 18:27:06 bulmer Exp

      subroutine bigden (n,npt,xopt,xpt,bmat,zmat,idz,ndim,kopt, 
     .              knew,d,w,vlag,beta,s,wvec,prod)

      implicit real*8 (a-h,o-z), integer (i-n)

      dimension xopt(*),xpt(npt,*),bmat(ndim,*),zmat(npt,*),d(*)
      dimension w(*),vlag(*),s(*),wvec(ndim,*),prod(ndim,*)
      dimension den(9),denex(9),par(9)

    # N is the number of variables.
    # NPT is the number of interpolation equations.
    # XOPT is the best interpolation point so far.
    # XPT contains the coordinates of the current interpolation points.
    # BMAT provides the last N columns of H.
    # ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
    # NDIM is the first dimension of BMAT and has the value NPT+N.
    # KOPT is the index of the optimal interpolation point.
    # KNEW is the index of the interpolation point that is going to be moved.
    # D will be set to the step from XOPT to the new point, and on entry it
    #   should be the D that was calculated by the last call of BIGLAG. The
    #   length of the initial D provides a trust region bound on the final D.
    # W will be set to Wcheck for the final choice of D.
    # VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
    # BETA will be set to the value that will occur in the updating formula
    #   when the KNEW-th interpolation point is moved to its new position.
    # S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used
    #   for working space.
    #
    # D is calculated in a way that should provide a denominator with a large
    # modulus in the updating formula when the KNEW-th interpolation point is
    # shifted to the new position XOPT+D.

      half=0.5
      one=1.0
      quart=0.25
      two=2.0
      zero=0.0
      twopi=8.0*atan(one)
      nptm=npt-n-1

    # Store the first NPT elements of the KNEW-th column of H in W(N+1)
    # to W(N+NPT).

      do k=1,npt
        w(n+k)=zero
      enddo
      do j=1,nptm
        temp=zmat(knew,j)
        if (j < idz) temp=-temp
        do k=1,npt
            w(n+k)=w(n+k)+temp*zmat(k,j)
        enddo
      enddo
      alpha=w(n+knew)

    # The initial search direction D is taken from the last call of BIGLAG,
    # and the initial S is set below, usually to the direction from X_OPT
    # to X_KNEW, but a different direction to an interpolation point may
    # be chosen, in order to prevent S from being nearly parallel to D.

      dd=zero
      ds=zero
      ss=zero
      xoptsq=zero
      do i=1,n
        dd=dd+d(i)**2
        s(i)=xpt(knew,i)-xopt(i)
        ds=ds+d(i)*s(i)
        ss=ss+s(i)**2
        xoptsq=xoptsq+xopt(i)**2
      enddo
      if (ds*ds > 0.99*dd*ss) then
        ksav=knew
        dtest=ds*ds/ss
        do k=1,npt
            if (k .ne. kopt) then
                dstemp=zero
                sstemp=zero
                do i=1,n
                    diff=xpt(k,i)-xopt(i)
                    dstemp=dstemp+d(i)*diff
                    sstemp=sstemp+diff*diff
                enddo
                if (dstemp*dstemp/sstemp < dtest) then
                    ksav=k
                    dtest=dstemp*dstemp/sstemp
                    ds=dstemp
                    ss=sstemp
                endif
            endif
        enddo
        do i=1,n
            s(i)=xpt(ksav,i)-xopt(i)
        enddo
      endif
      ssden=dd*ss-ds*ds
      iterc=0
      densav=zero

    # Begin the iteration by overwriting S with a vector that has the
    # required length and direction.

70    iterc=iterc+1
      temp=one/sqrt(ssden)
      xoptd=zero
      xopts=zero
      do i=1,n
        s(i)=temp*(dd*s(i)-ds*d(i))
        xoptd=xoptd+xopt(i)*d(i)
        xopts=xopts+xopt(i)*s(i)
      enddo

    # Set the coefficients of the first two terms of BETA.

      tempa=half*xoptd*xoptd
      tempb=half*xopts*xopts
      den(1)=dd*(xoptsq+half*dd)+tempa+tempb
      den(2)=two*xoptd*dd
      den(3)=two*xopts*dd
      den(4)=tempa-tempb
      den(5)=xoptd*xopts
      do i=6,9
        den(i)=zero
      enddo

    # Put the coefficients of Wcheck in WVEC.

      do k=1,npt
        tempa=zero
        tempb=zero
        tempc=zero
        do i=1,n
            tempa=tempa+xpt(k,i)*d(i)
            tempb=tempb+xpt(k,i)*s(i)
            tempc=tempc+xpt(k,i)*xopt(i)
        enddo
        wvec(k,1)=quart*(tempa*tempa+tempb*tempb)
        wvec(k,2)=tempa*tempc
        wvec(k,3)=tempb*tempc
        wvec(k,4)=quart*(tempa*tempa-tempb*tempb)
        wvec(k,5)=half*tempa*tempb
      enddo
      do i=1,n
        ip=i+npt
        wvec(ip,1)=zero
        wvec(ip,2)=d(i)
        wvec(ip,3)=s(i)
        wvec(ip,4)=zero
        wvec(ip,5)=zero
      enddo

    # Put the coefficents of THETA*Wcheck in PROD.

      do jc=1,5
        nw=npt
        if (jc == 2 .or. jc == 3) nw=ndim
        do k=1,npt
            prod(k,jc)=zero
        enddo
        do j=1,nptm
            sum=zero
            do k=1,npt
                sum=sum+zmat(k,j)*wvec(k,jc)
            enddo
            if (j < idz) sum=-sum
            do k=1,npt
                prod(k,jc)=prod(k,jc)+sum*zmat(k,j)
            enddo
        enddo
        if (nw == ndim) then
            do k=1,npt
                sum=zero
                do j=1,n
                    sum=sum+bmat(k,j)*wvec(npt+j,jc)
                enddo
                prod(k,jc)=prod(k,jc)+sum
            enddo
        endif
        do  j=1,n
            sum=zero
            do i=1,nw
                sum=sum+bmat(i,j)*wvec(i,jc)
            enddo
            prod(npt+j,jc)=sum
        enddo
      enddo

    # Include in DEN the part of BETA that depends on THETA.

      do k=1,ndim
        sum=zero
        do i=1,5
            par(i)=half*prod(k,i)*wvec(k,i)
            sum=sum+par(i)
        enddo
        den(1)=den(1)-par(1)-sum
        tempa=prod(k,1)*wvec(k,2)+prod(k,2)*wvec(k,1)
        tempb=prod(k,2)*wvec(k,4)+prod(k,4)*wvec(k,2)
        tempc=prod(k,3)*wvec(k,5)+prod(k,5)*wvec(k,3)
        den(2)=den(2)-tempa-half*(tempb+tempc)
        den(6)=den(6)-half*(tempb-tempc)
        tempa=prod(k,1)*wvec(k,3)+prod(k,3)*wvec(k,1)
        tempb=prod(k,2)*wvec(k,5)+prod(k,5)*wvec(k,2)
        tempc=prod(k,3)*wvec(k,4)+prod(k,4)*wvec(k,3)
        den(3)=den(3)-tempa-half*(tempb-tempc)
        den(7)=den(7)-half*(tempb+tempc)
        tempa=prod(k,1)*wvec(k,4)+prod(k,4)*wvec(k,1)
        den(4)=den(4)-tempa-par(2)+par(3)
        tempa=prod(k,1)*wvec(k,5)+prod(k,5)*wvec(k,1)
        tempb=prod(k,2)*wvec(k,3)+prod(k,3)*wvec(k,2)
        den(5)=den(5)-tempa-half*tempb
        den(8)=den(8)-par(4)+par(5)
        tempa=prod(k,4)*wvec(k,5)+prod(k,5)*wvec(k,4)
        den(9)=den(9)-half*tempa
      enddo

      # Extend DEN so that it holds all the coefficients of DENOM.

      sum=zero
      do i=1,5
        par(i)=half*prod(knew,i)**2
        sum=sum+par(i)
      enddo
      denex(1)=alpha*den(1)+par(1)+sum
      tempa=two*prod(knew,1)*prod(knew,2)
      tempb=prod(knew,2)*prod(knew,4)
      tempc=prod(knew,3)*prod(knew,5)
      denex(2)=alpha*den(2)+tempa+tempb+tempc
      denex(6)=alpha*den(6)+tempb-tempc
      tempa=two*prod(knew,1)*prod(knew,3)
      tempb=prod(knew,2)*prod(knew,5)
      tempc=prod(knew,3)*prod(knew,4)
      denex(3)=alpha*den(3)+tempa+tempb-tempc
      denex(7)=alpha*den(7)+tempb+tempc
      tempa=two*prod(knew,1)*prod(knew,4)
      denex(4)=alpha*den(4)+tempa+par(2)-par(3)
      tempa=two*prod(knew,1)*prod(knew,5)
      denex(5)=alpha*den(5)+tempa+prod(knew,2)*prod(knew,3)
      denex(8)=alpha*den(8)+par(4)-par(5)
      denex(9)=alpha*den(9)+prod(knew,4)*prod(knew,5)

    # Seek the value of the angle that maximizes the modulus of DENOM.

      sum=denex(1)+denex(2)+denex(4)+denex(6)+denex(8)
      denold=sum
      denmax=sum
      isave=0
      iu=49
      temp=twopi/float(iu+1)
      par(1)=one
      do i=1,iu
        angle=float(i)*temp
        par(2)=cos(angle)
        par(3)=sin(angle)
        do j=4,8,2
            par(j)=par(2)*par(j-2)-par(3)*par(j-1)
            par(j+1)=par(2)*par(j-1)+par(3)*par(j-2)
        enddo
        sumold=sum
        sum=zero
        do j=1,9
            sum=sum+denex(j)*par(j)
        enddo
        if (abs(sum) > abs(denmax)) then
            denmax=sum
            isave=i
            tempa=sumold
        elseif (i == isave+1) then
            tempb=sum
        endif
      enddo
      if (isave == 0) tempa=sum
      if (isave == iu) tempb=denold
      step=zero
      if (tempa .ne. tempb) then
        tempa=tempa-denmax
        tempb=tempb-denmax
        step=half*(tempa-tempb)/(tempa+tempb)
      endif
      angle=temp*(float(isave)+step)

    # Calculate the new parameters of the denominator, the new VLAG vector
    # and the new D. Then test for convergence.

      par(2)=cos(angle)
      par(3)=sin(angle)
      do j=4,8,2
        par(j)=par(2)*par(j-2)-par(3)*par(j-1)
        par(j+1)=par(2)*par(j-1)+par(3)*par(j-2)
      enddo
      beta=zero
      denmax=zero
      do j=1,9
        beta=beta+den(j)*par(j)
        denmax=denmax+denex(j)*par(j)
      enddo
      do k=1,ndim
        vlag(k)=zero
        do j=1,5
            vlag(k)=vlag(k)+prod(k,j)*par(j)
        enddo
      enddo
      tau=vlag(knew)
      dd=zero
      tempa=zero
      tempb=zero
      do i=1,n
        d(i)=par(2)*d(i)+par(3)*s(i)
        w(i)=xopt(i)+d(i)
        dd=dd+d(i)**2
        tempa=tempa+d(i)*w(i)
        tempb=tempb+w(i)*w(i)
      enddo
      if (iterc >= n) goto 340
      if (iterc > 1) densav=max(densav,denold)
      if (abs(denmax) <= 1.1*abs(densav)) goto 340
      densav=denmax

    # Set S to half the gradient of the denominator with respect to D.
    # Then branch for the next iteration.

      do i=1,n
        temp=tempa*xopt(i)+tempb*d(i)-vlag(npt+i)
        s(i)=tau*bmat(knew,i)+alpha*temp
      enddo
      do k=1,npt
        sum=zero
        do j=1,n
            sum=sum+xpt(k,j)*w(j)
        enddo
        temp=(tau*w(n+k)-alpha*vlag(k))*sum
        do i=1,n
            s(i)=s(i)+temp*xpt(k,i)
        enddo
      enddo
      ss=zero
      ds=zero
      do i=1,n
        ss=ss+s(i)**2
        ds=ds+d(i)*s(i)
      enddo
      ssden=dd*ss-ds*ds
      if (ssden >= 1.0e-8*dd*ss) goto 70

    # Set the vector W before the RETURN from the subroutine.

340   do k=1,ndim
          w(k)=zero
          do j=1,5
            w(k)=w(k)+wvec(k,j)*par(j)
          enddo
      enddo
      vlag(kopt)=vlag(kopt)+one

      return

      end # bigden

c----------------------------------------------------------------------
# biglag.m,v 1.3 2006/08/12 18:27:06 bulmer Exp

      subroutine biglag(n,npt,xopt,xpt,bmat,zmat,idz,ndim,knew, 
     .             delta,d,alpha,hcol,gc,gd,s,w)

      implicit real (a-h,o-z), integer (i-n)
      dimension xopt(*),xpt(npt,*),bmat(ndim,*),zmat(npt,*),d(*)
      dimension hcol(*),gc(*),gd(*),s(*),w(*)

    # N is the number of variables.
    # NPT is the number of interpolation equations.
    # XOPT is the best interpolation point so far.
    # XPT contains the coordinates of the current interpolation points.
    # BMAT provides the last N columns of H.
    # ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
    # NDIM is the first dimension of BMAT and has the value NPT+N.
    # KNEW is the index of the interpolation point that is going to be moved.
    # DELTA is the current trust region bound.
    # D will be set to the step from XOPT to the new point.
    # ALPHA will be set to the KNEW-th diagonal element of the H matrix.
    # HCOL, GC, GD, S and W will be used for working space.
    #
    # The step D is calculated in a way that attempts to maximize the modulus
    # of LFUNC(XOPT+D), subject to the bound ||D|| <= DELTA, where LFUNC is
    # the KNEW-th Lagrange function.

      half=0.5
      one=1.0
      zero=0.0
      twopi=8.0*atan(one)
      delsq=delta*delta
      nptm=npt-n-1

    # Set the first NPT components of HCOL to the leading elements of the
    # KNEW-th column of H.

      iterc=0
      do k=1,npt
        hcol(k)=zero
      enddo
      do j=1,nptm
        temp=zmat(knew,j)
        if (j < idz) temp=-temp
        do k=1,npt
            hcol(k)=hcol(k)+temp*zmat(k,j)
        enddo
      enddo
      alpha=hcol(knew)

    # Set the unscaled initial direction D. Form the gradient of LFUNC at
    # XOPT, and multiply D by the second derivative matrix of LFUNC.

      dd=zero
      do i=1,n
        d(i)=xpt(knew,i)-xopt(i)
        gc(i)=bmat(knew,i)
        gd(i)=zero
        dd=dd+d(i)**2
      enddo
      do k=1,npt
        temp=zero
        sum=zero
        do j=1,n
            temp=temp+xpt(k,j)*xopt(j)
            sum=sum+xpt(k,j)*d(j)
        enddo
        temp=hcol(k)*temp
        sum=hcol(k)*sum
        do i=1,n
            gc(i)=gc(i)+temp*xpt(k,i)
            gd(i)=gd(i)+sum*xpt(k,i)
        enddo
      enddo

    # Scale D and GD, with a sign change if required. Set S to another
    # vector in the initial two dimensional subspace.

      gg=zero
      sp=zero
      dhd=zero
      do i=1,n
        gg=gg+gc(i)**2
        sp=sp+d(i)*gc(i)
        dhd=dhd+d(i)*gd(i)
      enddo
      scale=delta/dsqrt(dd)
      if (sp*dhd < zero) scale=-scale
      temp=zero
      if (sp*sp > 0.99*dd*gg) temp=one
      tau=scale*(abs(sp)+half*scale*abs(dhd))
      if (gg*delsq < 0.01*tau*tau) temp=one
      do i=1,n
        d(i)=scale*d(i)
        gd(i)=scale*gd(i)
        s(i)=gc(i)+temp*gd(i)
      enddo

    # Begin the iteration by overwriting S with a vector that has the
    # required length and direction, except that termination occurs if
    # the given D and S are nearly parallel.

 80   iterc=iterc+1
      dd=zero
      sp=zero
      ss=zero
      do i=1,n
        dd=dd+d(i)**2
        sp=sp+d(i)*s(i)
        ss=ss+s(i)**2
      enddo
      temp=dd*ss-sp*sp
      if (temp <= 1.0e-8*dd*ss) goto 160
      denom=sqrt(temp)
      do i=1,n
        s(i)=(dd*s(i)-sp*d(i))/denom
        w(i)=zero
      enddo

    # Calculate the coefficients of the objective function on the circle,
    # beginning with the multiplication of S by the second derivative matrix.

      do k=1,npt
        sum=zero
        do j=1,n
            sum=sum+xpt(k,j)*s(j)
        enddo
        sum=hcol(k)*sum
        do i=1,n
            w(i)=w(i)+sum*xpt(k,i)
        enddo
      enddo
      cf1=zero
      cf2=zero
      cf3=zero
      cf4=zero
      cf5=zero
      do i=1,n
        cf1=cf1+s(i)*w(i)
        cf2=cf2+d(i)*gc(i)
        cf3=cf3+s(i)*gc(i)
        cf4=cf4+d(i)*gd(i)
        cf5=cf5+s(i)*gd(i)
      enddo
      cf1=half*cf1
      cf4=half*cf4-cf1

    # Seek the value of the angle that maximizes the modulus of TAU.

      taubeg=cf1+cf2+cf4
      taumax=taubeg
      tauold=taubeg
      isave=0
      iu=49
      temp=twopi/float(iu+1)
      do i=1,iu
        angle=float(i)*temp
        cth=cos(angle)
        sth=sin(angle)
        tau=cf1+(cf2+cf4*cth)*cth+(cf3+cf5*cth)*sth
        if (abs(tau) > abs(taumax)) then
            taumax=tau
            isave=i
            tempa=tauold
        elseif (i == isave+1) then
            tempb=tau
        endif
        tauold=tau
      enddo
      if (isave == 0) tempa=tau
      if (isave == iu) tempb=taubeg
      step=zero
      if (tempa .ne. tempb) then
        tempa=tempa-taumax
        tempb=tempb-taumax
        step=half*(tempa-tempb)/(tempa+tempb)
      endif
      angle=temp*(float(isave)+step)

    # Calculate the new D and GD. Then test for convergence.

      cth=cos(angle)
      sth=sin(angle)
      tau=cf1+(cf2+cf4*cth)*cth+(cf3+cf5*cth)*sth
      do i=1,n
        d(i)=cth*d(i)+sth*s(i)
        gd(i)=cth*gd(i)+sth*w(i)
        s(i)=gc(i)+gd(i)
      enddo

      if (dabs(tau) <= 1.1*abs(taubeg)) goto 160
      if (iterc < n) goto 80
160   return

      end # biglag

c----------------------------------------------------------------------
# trsapp.m,v 1.5 2006/08/12 18:27:06 bulmer Exp

      subroutine trsapp(n,npt,xopt,xpt,gq,hq,pq,delta,step,d,g,hd,hs,crvmin)

      implicit real (a-h,o-z), integer (i-n)
      dimension xopt(*),xpt(npt,*),gq(*),hq(*),pq(*),step(*),d(*),g(*),hd(*),hs(*)

    # N is the number of variables of a quadratic objective function, Q say.
    # The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
    #   in order to define the current quadratic model Q.
    # DELTA is the trust region radius, and has to be positive.
    # STEP will be set to the calculated trial step.
    # The arrays D, G, HD and HS will be used for working space.
    # CRVMIN will be set to the least curvature of H along the conjugate
    # directions that occur, except that it is set to zero if STEP goes
    # all the way to the trust region boundary.
    #
    # The calculation of STEP begins with the truncated conjugate gradient
    # method. If the boundary of the trust region is reached, then further
    # changes to STEP may be made, each one being in the 2D space spanned
    # by the current STEP and the corresponding gradient of Q. Thus STEP
    # should provide a substantial reduction to Q within the trust region.
    #
    # Initialization, which includes setting HD to H times XOPT.

      half=0.5
      zero=0.0
      twopi=8.0*atan(1.0)
      delsq=delta*delta
      iterc=0
      itermax=n
      itersw=itermax
      do i=1,n
          d(i)=xopt(i)
      enddo
      goto 170

    # Prepare for the first line search.

 20   qred=zero
      dd=zero
      do i=1,n
        step(i)=zero
        hs(i)=zero
        g(i)=gq(i)+hd(i)
        d(i)=-g(i)
        dd=dd+d(i)**2
      enddo
      crvmin=zero
      if (dd == zero) goto 160
      ds=zero
      ss=zero
      gg=dd
      ggbeg=gg

    # Calculate the step to the trust region boundary and the product HD.

 40   iterc=iterc+1
      temp=delsq-ss
      arg=ds*ds+dd*temp
      if (arg <= 0.0) then
        call remark('***FATAL in UOA (trsapp: divide-by-zero, reduce rhobeg')
        call kaboom(0)
      endif
      bstep=temp/(ds+sqrt(ds*ds+dd*temp))
      goto 170
 50   dhd=zero
      do j=1,n
        dhd=dhd+d(j)*hd(j)
      enddo

    # Update CRVMIN and set the step-length ALPHA.

      alpha=bstep
      if (dhd > zero) then
        temp=dhd/dd
        if (iterc == 1) crvmin=temp
        crvmin=min(crvmin,temp)
        alpha=min(alpha,gg/dhd)
      endif
      qadd=alpha*(gg-half*alpha*dhd)
      qred=qred+qadd

    #  Update STEP and HS.

      ggsav=gg
      gg=zero
      do i=1,n
        step(i)=step(i)+alpha*d(i)
        hs(i)=hs(i)+alpha*hd(i)
        gg=gg+(g(i)+hs(i))**2
      enddo

    #  Begin another conjugate direction iteration if required.

      if (alpha < bstep) then
        if (qadd <= 0.01*qred) goto 160
        if (gg <= 0.0001*ggbeg) goto 160
        if (iterc == itermax) goto 160
        temp=gg/ggsav
        dd=zero
        ds=zero
        ss=zero
        do i=1,n
            d(i)=temp*d(i)-g(i)-hs(i)
            dd=dd+d(i)**2
            ds=ds+d(i)*step(i)
            ss=ss+step(i)**2
        enddo
        if (ds <= zero) goto 160
        if (ss < delsq) goto 40
      endif
      crvmin=zero
      itersw=iterc

    # Test whether an alternative iteration is required.

 90   if (gg <= 1.0e-4*ggbeg) goto 160
      sg=zero
      shs=zero
      do i=1,n
        sg=sg+step(i)*g(i)
        shs=shs+step(i)*hs(i)
      enddo
      sgk=sg+shs
      angtest=sgk/sqrt(gg*delsq)
      if (angtest <= -0.99) goto 160

    #  Begin the alternative iteration by calculating D and HD and some
    #  scalar products.

      iterc=iterc+1
      temp=sqrt(delsq*gg-sgk*sgk)
      tempa=delsq/temp
      tempb=sgk/temp
      do i=1,n
        d(i)=tempa*(g(i)+hs(i))-tempb*step(i)
      enddo
      goto 170
120   dg=zero
      dhd=zero
      dhs=zero
      do i=1,n
        dg=dg+d(i)*g(i)
        dhd=dhd+hd(i)*d(i)
        dhs=dhs+hd(i)*step(i)
      enddo

    # Seek the value of the angle that minimizes Q.

      cf=half*(shs-dhd)
      qbeg=sg+cf
      qsav=qbeg
      qmin=qbeg
      isave=0
      iu=49
      temp=twopi/dfloat(iu+1)
      do i=1,iu
        angle=dfloat(i)*temp
        cth=cos(angle)
        sth=sin(angle)
        qnew=(sg+cf*cth)*cth+(dg+dhs*cth)*sth
        if (qnew < qmin) then
            qmin=qnew
            isave=i
            tempa=qsav
        elseif (i == isave+1) then
            tempb=qnew
        endif
        qsav=qnew
      enddo
      if (isave == zero) tempa=qnew
      if (isave == iu) tempb=qbeg
      angle=zero
      if (tempa .ne. tempb) then
        tempa=tempa-qmin
        tempb=tempb-qmin
        angle=half*(tempa-tempb)/(tempa+tempb)
      endif
      angle=temp*(dfloat(isave)+angle)
    
    # Calculate the new STEP and HS. Then test for convergence.

      cth=cos(angle)
      sth=sin(angle)
      reduc=qbeg-(sg+cf*cth)*cth-(dg+dhs*cth)*sth
      gg=zero
      do i=1,n
        step(i)=cth*step(i)+sth*d(i)
        hs(i)=cth*hs(i)+sth*hd(i)
        gg=gg+(g(i)+hs(i))**2
      enddo
      qred=qred+reduc
      ratio=reduc/qred
      if (iterc < itermax .and. ratio > 0.01) goto 90
160   return

    # The following instructions act as a subroutine for setting the vector
    # HD to the vector D multiplied by the second derivative matrix of Q.
    # They are called from three different places, which are distinguished
    # by the value of ITERC.

170   do i=1,n
        hd(i)=zero
      enddo
      do k=1,npt
        temp=zero
        do j=1,n
            temp=temp+xpt(k,j)*d(j)
        enddo
        temp=temp*pq(k)
        do i=1,n
            hd(i)=hd(i)+temp*xpt(k,i)
        enddo
      enddo
      ih=0
      do j=1,n
        do i=1,j
            ih=ih+1
            if (i < j) hd(j)=hd(j)+hq(ih)*d(i)
            hd(i)=hd(i)+hq(ih)*d(j)
        enddo
      enddo
      if (iterc == 0) goto 20
      if (iterc <= itersw) goto 50
      goto 120
    
      end # trsapp

c----------------------------------------------------------------------
# update.m,v 1.3 2006/08/12 18:27:06 bulmer Exp

      subroutine update(n,npt,bmat,zmat,idz,ndim,vlag,beta,knew,w)

      implicit real (a-h,o-z), integer (i-n)
      dimension bmat(ndim,*),zmat(npt,*),vlag(*),w(*)

    # The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
    # interpolation point that has index KNEW. On entry, VLAG contains the
    # components of the vector Theta*Wcheck+e_b of the updating formula
    # (6.11), and BETA holds the value of the parameter that has this name.
    # The vector W is used for working space.
    # 
    # Set some constants.

      one=1.0
      zero=0.0
      nptm=npt-n-1

    # Apply the rotations that put zeros in the KNEW-th row of ZMAT.

      jl=1
      do j=2,nptm
        if (j == idz) then
            jl=idz
        elseif (zmat(knew,j) .ne. zero) then
            temp=sqrt(zmat(knew,jl)**2+zmat(knew,j)**2)
            tempa=zmat(knew,jl)/temp
            tempb=zmat(knew,j)/temp
            do i=1,npt
                temp=tempa*zmat(i,jl)+tempb*zmat(i,j)
                zmat(i,j)=tempa*zmat(i,j)-tempb*zmat(i,jl)
                zmat(i,jl)=temp
            enddo
            zmat(knew,j)=zero
        endif
      enddo

    # Put the first NPT components of the KNEW-th column of HLAG into W,
    # and calculate the parameters of the updating formula.

      tempa=zmat(knew,1)
      if (idz >= 2) tempa=-tempa
      if (jl > 1) tempb=zmat(knew,jl)
      do i=1,npt
        w(i)=tempa*zmat(i,1)
        if (jl > 1) w(i)=w(i)+tempb*zmat(i,jl)
      enddo
      alpha=w(knew)
      tau=vlag(knew)
      tausq=tau*tau
      denom=alpha*beta+tausq
      vlag(knew)=vlag(knew)-one
    
    # Complete the updating of ZMAT when there is only one nonzero element
    # in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
    # then the first column of ZMAT will be exchanged with another one later.

      iflag=0
      if (jl == 1) then
        temp=sqrt(abs(denom))
        tempb=tempa/temp
        tempa=tau/temp
        do i=1,npt
            zmat(i,1)=tempa*zmat(i,1)-tempb*vlag(i)
        enddo
        if (idz == 1 .and. temp < zero) idz=2
        if (idz >= 2 .and. temp >= zero) iflag=1
      else
        # Complete the updating of ZMAT in the alternative case.
        ja=1
        if (beta >= zero) ja=jl
        jb=jl+1-ja
        temp=zmat(knew,jb)/denom
        tempa=temp*beta
        tempb=temp*tau
        temp=zmat(knew,ja)
        scala=one/sqrt(abs(beta)*temp*temp+tausq)
        scalb=scala*sqrt(abs(denom))
        do i=1,npt
            zmat(i,ja)=scala*(tau*zmat(i,ja)-temp*vlag(i))
            zmat(i,jb)=scalb*(zmat(i,jb)-tempa*w(i)-tempb*vlag(i))
        enddo
        if (denom <= zero) then
            if (beta < zero) idz=idz+1
            if (beta >= zero) iflag=1
        endif
      endif

    # IDZ is reduced in the following case, and usually the first column
    # of ZMAT is exchanged with a later one.

      if (iflag == 1) then
        idz=idz-1
        do i=1,npt
            temp=zmat(i,1)
            zmat(i,1)=zmat(i,idz)
            zmat(i,idz)=temp
        enddo
      endif

    # Finally, update the matrix BMAT.

      do j=1,n
        jp=npt+j
        w(jp)=bmat(knew,j)
        tempa=(alpha*vlag(jp)-tau*w(jp))/denom
        tempb=(-beta*w(jp)-tau*vlag(jp))/denom
        do i=1,jp
            bmat(i,j)=bmat(i,j)+tempa*vlag(i)+tempb*w(i)
            if (i > npt) bmat(jp,i-npt)=bmat(i,j)
        enddo
      enddo

      return

      end # update

c----------------------------------------------------------------------
# calfun.m,v 1.2 2004/10/10 19:51:43 bulmer Exp


      subroutine calfun(n,x,f)

      implicit real (a-h,o-z), integer (i-n)
      Use(UOA)

    # This routine, which was the user-supplied object function evaluator
    # in the original Fortran NEWUOA, now executes a Basis script function
    # which is a wrapper for the new user-supplied script (of arbitrary
    # name) that evaluates the object function.

      dimension x(*)

      do i=1,n
        x_uoa(i)=x(i)
      enddo
c_tdr      call parsestr("f_uoa=objfcn_newuoa(x_uoa)")
      call remark('***Must provide functionality just above to use')
      call remark('***parsestr does not work for c++ uedge')
      call kaboom(0)
      f=f_uoa

      return

      end # calfun
