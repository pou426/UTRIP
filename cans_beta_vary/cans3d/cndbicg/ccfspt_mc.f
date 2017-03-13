c======================================================================|
      subroutine ccfspt_mc(cmat,src,rkap0,gm,dt,te,ro,bx,by,bz
     &              ,x,xh,dx,dxm,ix,dy,dym,jx,dz,dzm,kx)
c======================================================================|
c 
c NAME  ccfspt_mc
c 
c PURPOSE
c    set coefficient matrix for conduction equation
c        * Spitzer type
c        * magnetic field
c        * Cylindrical coordinate
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
c    bx(ix,jx,kx) : [double] magnetic field 
c    by(ix,jx,kx) : [double] magnetic field
c    bz(ix,jx,kx) : [double] magnetic field 
c    x(ix), xm(ix): [double] coordinate
c    margin: [integer] margin, i.e. # of grid points outside the boundary
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,kx,7),src(ix,jx,kx)
      dimension ro(ix,jx,kx)
      dimension te(ix,jx,kx)
      dimension teh(ix,jx,kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension a0t(ix,jx,kx),axx(ix,jx,kx),axy(ix,jx,kx),axz(ix,jx,kx)
     &                       ,ayy(ix,jx,kx),ayz(ix,jx,kx),ayx(ix,jx,kx)
     &                       ,azz(ix,jx,kx),azx(ix,jx,kx),azy(ix,jx,kx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension dz(kx),dzm(kx)
      dimension x(ix),xh(ix)
c----------------------------------------------------------------------|

      do k=1,kx
      do j=1,jx
      do i=1,ix
         a0t(i,j,k)=ro(i,j,k)/(gm-1.)*x(i)
      enddo
      enddo
      enddo

      telim=1.e4
      sft=0.03

      do k=1,kx
      do j=1,jx
      do i=1,ix-1
c        te00=0.5*(te(i,j,k)+te(i+1,j,k))
         te00=sqrt(te(i,j,k)*te(i+1,j,k))
         bx00=0.5*(bx(i,j,k)+bx(i+1,j,k))
         by00=0.5*(by(i,j,k)+by(i+1,j,k))
         bz00=0.5*(bz(i,j,k)+bz(i+1,j,k))
         bf2=(sqrt(bx00**2+by00**2+bz00**2)+sft)**2
         tmp=rkap0*sqrt(min(te00,telim))**5*bx00/bf2*xh(i)
         axx(i,j,k)=tmp*bx00
         axy(i,j,k)=tmp*by00
         axz(i,j,k)=tmp*bz00
      enddo
      enddo
      enddo

      do k=1,kx
      do j=1,jx-1
      do i=1,ix
c        te00=0.5*(te(i,j,k)+te(i,j+1,k))
         te00=sqrt(te(i,j,k)*te(i,j+1,k))
         bx00=0.5*(bx(i,j,k)+bx(i,j+1,k))
         by00=0.5*(by(i,j,k)+by(i,j+1,k))
         bz00=0.5*(bz(i,j,k)+bz(i,j+1,k))
         bf2=(sqrt(bx00**2+by00**2+bz00**2)+sft)**2
         tmp=rkap0*sqrt(min(te00,telim))**5*by00/bf2*x(i)
         ayy(i,j,k)=tmp*by00
         ayz(i,j,k)=tmp*bz00
         ayx(i,j,k)=tmp*bx00
      enddo
      enddo
      enddo

      do k=1,kx-1
      do j=1,jx
      do i=1,ix
c        te00=0.5*(te(i,j,k)+te(i,j,k+1))
         te00=sqrt(te(i,j,k)*te(i,j,k+1))
         bx00=0.5*(bx(i,j,k)+bx(i,j,k+1))
         by00=0.5*(by(i,j,k)+by(i,j,k+1))
         bz00=0.5*(bz(i,j,k)+bz(i,j,k+1))
         bf2=(sqrt(bx00**2+by00**2+bz00**2)+sft)**2
         tmp=rkap0*sqrt(min(te00,telim))**5*bz00/bf2*x(i)
         azz(i,j,k)=tmp*bz00
         azx(i,j,k)=tmp*bx00
         azy(i,j,k)=tmp*by00
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
       do k=1,kx-1
       do j=1,jx-1
       do i=1,ix-1
        teh(i,j,k)= 0.125*
     &  ((te(i,j,k  )+te(i+1,j,k  ))+(te(i,j+1,k  )+te(i+1,j+1,k  ))
     & + (te(i,j,k+1)+te(i+1,j,k+1))+(te(i,j+1,k+1)+te(i+1,j+1,k+1)))
       enddo
       enddo
       enddo

       do k=2,kx-1
       do j=2,jx-1
       do i=2,ix-1
         src(i,j,k) =a0t(i,j,k)/dt*te(i,j,k)
     &   +axy(i  ,j,k)/dx(i)/dy(j)/2*((teh(i  ,j,k-1)-teh(i  ,j-1,k-1))
     &                               +(teh(i  ,j,k  )-teh(i  ,j-1,k  )))
     &   -axy(i-1,j,k)/dx(i)/dy(j)/2*((teh(i-1,j,k-1)-teh(i-1,j-1,k-1))
     &                               +(teh(i-1,j,k  )-teh(i-1,j-1,k  )))
     &   +axz(i  ,j,k)/dx(i)/dz(k)/2*((teh(i  ,j-1,k)-teh(i  ,j-1,k-1))
     &                               +(teh(i  ,j  ,k)-teh(i  ,j  ,k-1)))
     &   -axz(i-1,j,k)/dx(i)/dz(k)/2*((teh(i-1,j-1,k)-teh(i-1,j-1,k-1))
     &                               +(teh(i-1,j  ,k)-teh(i-1,j  ,k-1)))
     &   +ayz(i,j  ,k)/dy(j)/dz(k)/2*((teh(i-1,j  ,k)-teh(i-1,j  ,k-1))
     &                               +(teh(i  ,j  ,k)-teh(i  ,j  ,k-1)))
     &   -ayz(i,j-1,k)/dy(j)/dz(k)/2*((teh(i-1,j-1,k)-teh(i-1,j-1,k-1))
     &                               +(teh(i  ,j-1,k)-teh(i  ,j-1,k-1)))
     &   +ayx(i,j  ,k)/dy(j)/dx(i)/2*((teh(i,j  ,k-1)-teh(i-1,j  ,k-1))
     &                               +(teh(i,j  ,k  )-teh(i-1,j  ,k  )))
     &   -ayx(i,j-1,k)/dy(j)/dx(i)/2*((teh(i,j-1,k-1)-teh(i-1,j-1,k-1))
     &                               +(teh(i,j-1,k  )-teh(i-1,j-1,k  )))
     &   +azx(i,j,k  )/dz(k)/dx(i)/2*((teh(i,j-1,k  )-teh(i-1,j-1,k  ))
     &                               +(teh(i,j  ,k  )-teh(i-1,j  ,k  )))
     &   -azx(i,j,k-1)/dz(k)/dx(i)/2*((teh(i,j-1,k-1)-teh(i-1,j-1,k-1))
     &                               +(teh(i,j  ,k-1)-teh(i-1,j  ,k-1)))
     &   +azy(i,j,k  )/dz(k)/dy(j)/2*((teh(i-1,j,k  )-teh(i-1,j-1,k  ))
     &                               +(teh(i  ,j,k  )-teh(i  ,j-1,k  )))
     &   -azy(i,j,k-1)/dz(k)/dy(j)/2*((teh(i-1,j,k-1)-teh(i-1,j-1,k-1))
     &                               +(teh(i  ,j,k-1)-teh(i  ,j-1,k-1)))
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
