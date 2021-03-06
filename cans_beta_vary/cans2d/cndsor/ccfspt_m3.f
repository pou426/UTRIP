c======================================================================|
      subroutine ccfspt_m3(cmat,src,rkap0,gm,dt,te,ro,bx,by,bz
     &              ,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c 
c NAME  ccfspt_m3
c 
c PURPOSE
c    set coefficient matrix for conduction equation
c        * Spitzer type
c        * magnetic field 3 components
c 
c OUTPUTS
c    cmat(ix,jx,5): [double] coefficient matrix of heat conduction eq
c    src(ix,jx): [double] source vector of heat conduction eq
c 
c INPUTS
c    NOTE: ??m(ix,jx) is the variable array defined at grid bounds
c    ix,jx: [integer] dimension size
c    te(ix,jx): [double] temperature
c    ro(ix,jx): [double] density
c    dx(ix), dxm(ix): [double] grid spacing
c    dt: [double] delta time
c    rkap0: [double] constant part of heat conduction coefficient
c    gm: [double] polytropic index gamma
c    bx(ix,jx) : [double] magnetic field
c    by(ix,jx) : [double] magnetic field
c    bz(ix,jx) : [double] magnetic field
c 
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,5),src(ix,jx)
      dimension ro(ix,jx)
      dimension te(ix,jx)
      dimension teh(ix,jx)
      dimension bx(ix,jx),by(ix,jx),bz(ix,jx)
      dimension a0t(ix,jx),axx(ix,jx),axy(ix,jx)
     &                    ,ayy(ix,jx),ayx(ix,jx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix
         a0t(i,j)=ro(i,j)/(gm-1.)
      enddo
      enddo

      telim=1.e4
      sft=0.03

      do j=1,jx
      do i=1,ix-1
c        te00=0.5*(te(i,j)+te(i+1,j))
         te00=sqrt(te(i,j)*te(i+1,j))
         bx00=0.5*(bx(i,j)+bx(i+1,j))
         by00=0.5*(by(i,j)+by(i+1,j))
         bz00=0.5*(bz(i,j)+bz(i+1,j))
         bf2=(sqrt(bx00**2+by00**2+bz00**2)+sft)**2
         tmp=rkap0*sqrt(min(te00,telim))**5*bx00/bf2
         axx(i,j)=tmp*bx00
         axy(i,j)=tmp*by00
c        axy(i,j)=tmp*by00
c    &         *(1+(tanh((by00-sft)/0.01)+tanh(-(by00+sft)/0.01))/2)
      enddo
      enddo

      do j=1,jx-1
      do i=1,ix
c        te00=0.5*(te(i,j)+te(i,j+1))
         te00=sqrt(te(i,j)*te(i,j+1))
         bx00=0.5*(bx(i,j)+bx(i,j+1))
         by00=0.5*(by(i,j)+by(i,j+1))
         bz00=0.5*(bz(i,j)+bz(i,j+1))
         bf2=(sqrt(bx00**2+by00**2+bz00**2)+sft)**2
         tmp=rkap0*sqrt(min(te00,telim))**5*by00/bf2
         ayy(i,j)=tmp*by00
         ayx(i,j)=tmp*bx00
c        ayx(i,j)=tmp*bx00
c    &         *(1+(tanh((bx00-sft)/0.01)+tanh(-(bx00+sft)/0.01))/2)
      enddo
      enddo

c-----------------------------------------------------------------------
c    coefficient matrices
c    non-uniform section effect included
c    aa0=0 & src=0   for laprace equation
c    aa0=0 & src=!0  for poison equation
c    aa0=ccf(i)*dsh(i) & src = aa0*te    for diffusion equation
c-----------------------------------------------------------------------
       do j=1,jx-1
       do i=1,ix-1
          teh(i,j)=0.25*((te(i,j)+te(i+1,j))+(te(i,j+1)+te(i+1,j+1)))
       enddo
       enddo

       do j=2,jx-1
       do i=2,ix-1
          src(i,j) =a0t(i,j)/dt*te(i,j)
     &      +axy(i  ,j  )/dx(i)/dy(j)*(teh(i  ,j  )-teh(i  ,j-1))
     &      -axy(i-1,j  )/dx(i)/dy(j)*(teh(i-1,j  )-teh(i-1,j-1))
     &      +ayx(i  ,j  )/dy(j)/dx(i)*(teh(i  ,j  )-teh(i-1,j  ))
     &      -ayx(i  ,j-1)/dy(j)/dx(i)*(teh(i  ,j-1)-teh(i-1,j-1))
       enddo
       enddo

      do m=1,5
      do j=1,jx
      do i=1,ix
        cmat(i,j,m)=0.0d0
      enddo
      enddo
      enddo

       do j=2,jx-1
       do i=2,ix-1
          aa0=a0t(i,j)/dt

          aaxp=axx(i,j)/dxm(i)/dx(i)
          aaxm=axx(i-1,j)/dxm(i-1)/dx(i)
          cmat(i,j,2)=-aaxm
          cmat(i,j,3)=-aaxp

          aayp=ayy(i,j)/dym(j)/dy(j)
          aaym=ayy(i,j-1)/dym(j-1)/dy(j)
          cmat(i,j,4)=-aaym
          cmat(i,j,5)=-aayp

          cmat(i,j,1)=aa0+aaxp+aaxm+aayp+aaym
       enddo
       enddo

      return
      end
