c======================================================================|
      subroutine bbtoaa_s(az,bx,by,dxm,dym,x,y,margin,ix,jx)
c======================================================================|
c
c NAME  bbtoaa_s
c
c PURPOSE
c    calculate magnetic Vector potential
c        * Spherical coordinate, axis-symmetry
c
c OUTPUTS
c    az(ix,jx): [double] magnetic Vector potential
c
c INPUTS
c    bx(ix,jx),by(ix,jx): [double] magnetic field
c    dxm(ix),dy(jx) : [double] grid spacing
c    ix,jx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension dxm(ix),dym(jx)
      dimension x(ix),y(jx)
      dimension bx(ix,jx),by(ix,jx)
      dimension az(ix,jx)
c----------------------------------------------------------------------|

      az(margin+1,margin+1)=0.d0

      i=margin+1
      do j=margin+2,jx
        az(i,j) = az(i,j-1)
     &          + 0.5*x(i)**2*(sin(y(j-1))*bx(i,j-1)
     &                        +sin(y(j  ))*bx(i,j  ))*dym(j-1)
      enddo
      do j=margin,1,-1
        az(i,j) = az(i,j+1)
     &          - 0.5*x(i)**2*(sin(y(j  ))*bx(i,j  )
     &                        +sin(y(j+1))*bx(i,j+1))*dym(j  )
      enddo

      do j=1,jx
      do i=margin+2,ix
         az(i,j) = az(i-1,j)
     &          - 0.5*(x(i-1)*sin(y(j))*by(i-1,j)
     &                +x(i  )*sin(y(j))*by(i  ,j))*dxm(i-1)
      enddo
      do i=margin,1,-1
         az(i,j) = az(i+1,j)
     &          + 0.5*(x(i  )*sin(y(j))*by(i  ,j)
     &                +x(i+1)*sin(y(j))*by(i+1,j))*dxm(i  )
      enddo
      enddo

      do i=1,ix
      do j=1,jx
         az(i,j) = az(i,j)/( x(i)*sin(y(j)) )
      enddo
      enddo


      return
      end
