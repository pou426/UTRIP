c======================================================================|
      subroutine ccfspt(cmat,src,rkap0,gm,dt,te,ro
     &  ,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
c======================================================================|
c
c NAME  ccfspt
c
c PURPOSE
c    set coefficient matrix for conduction equation
c        * Spitzer type
c
c OUTPUTS
c    cmat(ix,jx,kx,7): [double] coefficient matrix of heat conduction eq
c    src(ix,jx,kx): [double] source vector of heat conduction eq
c
c INPUTS
c    NOTE: ??m(ix,jx,kx) is the variable array defined at grid bounds
c
c    ix,jx,kx: [integer] dimension size
c    te(ix,jx,kx): [double] temperature
c    ro(ix,jx,kx): [double] density
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    rkap0: [double] constant part of heat conduction coefficient
c    gm: [double] polytropic index gamma
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,kx,7),src(ix,jx,kx)
      dimension ro(ix,jx,kx)
      dimension te(ix,jx,kx)
      dimension a0t(ix,jx,kx),axx(ix,jx,kx),ayy(ix,jx,kx),azz(ix,jx,kx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)
c----------------------------------------------------------------------|

      do k=1,kx
      do j=1,jx
      do i=1,ix
         a0t(i,j,k)=ro(i,j,k)/(gm-1.)
      enddo
      enddo
      enddo

      telim=1.e4

      do k=1,kx
      do j=1,jx
      do i=1,ix-1
c        te00=0.5*(te(i,j,k)+te(i+1,j,k))
         te00=sqrt(te(i,j,k)*te(i+1,j,k))
         axx(i,j,k)=rkap0*sqrt(min(te00,telim))**5
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx-1
      do i=1,ix
c        te00=0.5*(te(i,j,k)+te(i,j+1,k))
         te00=sqrt(te(i,j,k)*te(i,j+1,k))
         ayy(i,j,k)=rkap0*sqrt(min(te00,telim))**5
      enddo
      enddo
      enddo

      do k=1,kx-1
      do j=1,jx
      do i=1,ix
c        te00=0.5*(te(i,j,k)+te(i,j+1,k))
         te00=sqrt(te(i,j,k)*te(i,j,k+1))
         azz(i,j,k)=rkap0*sqrt(min(te00,telim))**5
      enddo
      enddo
      enddo

c-----------------------------------------------------------------------
c    coefficient matrices
c    non-uniform section effect included
c    aa0=0 & src=0   for laprace equation
c    aa0=0 & src=!0  for poison equation
c    aa0=ccf(i)*dsh(i) & src = aa0*te    for diffusion equation
c-----------------------------------------------------------------------
       do k=2,kx-1
       do j=2,jx-1
       do i=2,ix-1
          src(i,j,k) =a0t(i,j,k)/dt*te(i,j,k)
       enddo
       enddo
       enddo

      do k=1,kx
      do m=1,7
      do j=1,jx
      do i=1,ix
        cmat(i,j,k,m)=0.0d0
      enddo
      enddo
      enddo
      enddo

       do k=2,kx-1
       do j=2,jx-1
       do i=2,ix-1
          aa0=a0t(i,j,k)/dt

          aaxp=axx(i,j,k)/dxm(i)/dx(i)
          aaxm=axx(i-1,j,k)/dxm(i-1)/dx(i)
          cmat(i,j,k,2)=-aaxm
          cmat(i,j,k,3)=-aaxp

          aayp=ayy(i,j,k)/dym(j)/dy(j)
          aaym=ayy(i,j-1,k)/dym(j-1)/dy(j)
          cmat(i,j,k,4)=-aaym
          cmat(i,j,k,5)=-aayp

          aazp=azz(i,j,k)/dzm(k)/dz(k)
          aazm=azz(i,j,k-1)/dzm(k-1)/dz(k)
          cmat(i,j,k,6)=-aazm
          cmat(i,j,k,7)=-aazp

          cmat(i,j,k,1)=aa0+aaxp+aaxm+aayp+aaym+aazp+aazm
       enddo
       enddo
       enddo

      return
      end
