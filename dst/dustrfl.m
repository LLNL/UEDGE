      subroutine rflmir(vx,vy,vz,amu,vnx,vny,vnz)
c
c---direction for mirror reflection
c   relatevely surface normal vector
c
      implicit none

      real*8 vx,vy,vz,amu,vnx,vny,vnz
      real*8 a2,s
      data a2/-2.d0/

      s=a2*amu

      vx=vx+s*vnx
      vy=vy+s*vny
      vz=vz+s*vnz

      s=1.0/dsqrt(vx*vx+vy*vy+vz*vz)
      vx=vx*s
      vy=vy*s
      vz=vz*s

      end


c-----------------------------------------------------------------------
      subroutine rfldif(wx,wy,wz,vnx,vny,vnz)

c
c-----direction for diffusive reflection
c     relatevely surface normal vector
c
      implicit none

      Use (std_cns)

      real*8 wx,wy,wz,vnx,vny,vnz
      real*8 seps,a,b,ax,az,amu,bmu,cf,sf

      data seps/1.d-9/
c
      az=1.d0-vnz*vnz
      if(az.lt.seps) go to 21
c
c-----deflection angles from normal
c     direction in local-normal-vector system
c
      call std_rndd7x(a)
      amu = dsqrt(a)
      bmu = dsqrt(1.d0 - a)

      az=dsqrt(az)
      ax=bmu/az
c
c---azimuthal angle
c
      call std_rndd7x(a)
      a=twopi*a
      sf=dsin(a)
      cf=dcos(a)
c
      a=amu-ax*cf*vnz
      b=ax*sf
c
c---laboratory system
c
      wz=amu*vnz+cf*bmu*az
      wx=a*vnx+b*vny
      wy=a*vny-b*vnx

      a=1.d0/dsqrt(wx*wx+wy*wy+wz*wz)
      wz=a*wz
      wx=a*wx
      wy=a*wy
      goto 1
c
c----laboratory & local-normal-vector systems coincide
c
  21  call std_rndd7x(a)
      amu = dsqrt(a)
      bmu = dsqrt(1.d0 - a)

      call std_rndd7x(a)
      a=twopi*a
      sf=dsin(a)
      cf=dcos(a)

      if(vnz.lt.0.d0) then
       wz=-amu
      else
       wz=amu
      endif

      wy=bmu*sf
      wx=bmu*cf

 1    return
      end


c-----------------------------------------------------------------------
      subroutine dirfix(wx,wy,wz, u0,phi0, vnx,vny,vnz)

c
c---direction for diffusive reflection
c   relatevely surface normal vector
c
      implicit none

      Use (std_cns)
      
      real*8 wx,wy,wz, u0,phi0, vnx,vny,vnz
      real*8 seps,a,b,ax,az,amu,bmu,cf,sf, u,phi

      data seps/1.d-9/
c
      if(u0 .lt. 0.d0) then
       u=-u0*u2pi
      else
       u=u0*u2pi
      endif
c
c---deflection angles from normal
c   direction in local-normal-vector system
c
      amu = dcos(u)
      if(amu .lt. 0.d0) amu=-amu
      bmu = 1.d0-amu*amu
      if(bmu .lt. 0.d0) bmu=-bmu
      bmu = dsqrt(bmu)
c
c---azimuthal angle cosines
c
      phi=phi0*u2pi

      sf=dsin(phi)
      cf=dcos(phi)

      az=1.d0-vnz*vnz
      if(az .lt. seps) go to 21

      az=dsqrt(az)
      ax=bmu/az

      a=amu-ax*cf*vnz
      b=ax*sf
c
c---laboratory system
c
      wz=amu*vnz+cf*bmu*az
      wx=a*vnx+b*vny
      wy=a*vny-b*vnx

      a=1.d0/dsqrt(wx*wx+wy*wy+wz*wz)
      wz=a*wz
      wx=a*wx
      wy=a*wy
      goto 1
c
c----laboratory & local-normal-vector systems coincide
c
  21  if(vnz .lt. 0.d0) then
       wz=-amu
      else
       wz=amu
      endif

      wy=bmu*sf
      wx=bmu*cf

 1    return
      end
