c======================================================================|
      subroutine ccfunf_s(cmat,src,rkap0,gm,dt,te,ro
     &                ,x,xm,y,ym,dx,dxm,ix,dy,dym,jx)
c======================================================================|
c
c NAME  ccfunf_s
c
c PURPOSE
c    set coefficient matrix for conduction equation
c        * Uniform
c        * Spherical coordinate, axis-symmetry 
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
c    x(ix), xm(ix): [double] coordinate
c    y(jx), ym(jx): [double] coordinate
c    margin: [integer] margin, i.e. # of grid points outside the boundary
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      dimension cmat(ix,jx,5),src(ix,jx)
      dimension ro(ix,jx)
      dimension te(ix,jx)
      dimension a0t(ix,jx),axx(ix,jx),ayy(ix,jx)
      dimension dx(ix),dxm(ix)
      dimension dy(jx),dym(jx)
      dimension x(ix),xm(ix)
      dimension y(jx),ym(jx)
c----------------------------------------------------------------------|

      do j=1,jx
      do i=1,ix
         a0t(i,j)=ro(i,j)/(gm-1.)*x(i)**2*sin(y(j))
      enddo
      enddo


      do j=1,jx
      do i=1,ix-1
         axx(i,j)=rkap0*xm(i)**2*sin(y(j))
      enddo
      enddo

      do j=1,jx-1
      do i=1,ix
         ayy(i,j)=rkap0*x(i)*sin(ym(j))
      enddo
      enddo

c-----------------------------------------------------------------------
c    coefficient matrices
c    non-uniform section effect included
c    aa0=0 & src=0   for laprace equation
c    aa0=0 & src=!0  for poison equation
c    aa0=ccf(i)*dsh(i) & src = aa0*te    for diffusion equation
c-----------------------------------------------------------------------
       do j=2,jx-1
       do i=2,ix-1
          src(i,j) =a0t(i,j)/dt*te(i,j)
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
