c======================================================================|
      subroutine bbtocc(cx,cy,cz,bx,by,bz,dx,dy,dz,ix,jx,kx,mfdim)
c======================================================================|
c
c NAME  bbtocx
c
c PURPOSE
c    calculate current density
c        * x-component
c
c OUTPUTS
c    cx(ix,jx,kx): [double] current density
c
c INPUTS
c    bz(ix,jx,kx): [double] magnetic field
c    dy(jx) : [double] grid spacing
c    ix,jx,kx: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      dimension mfdim(3)
      dimension dx(ix),dy(jx),dz(kx)
      dimension bx(ix,jx,kx),by(ix,jx,kx),bz(ix,jx,kx)
      dimension cx(ix,jx,kx),cy(ix,jx,kx),cz(ix,jx,kx)
c----------------------------------------------------------------------|

      do k=1,kx
      do j=1,jx
      do i=1,ix
        cx(i,j,k) = 0.d0
        cy(i,j,k) = 0.d0
        cz(i,j,k) = 0.d0
      enddo
      enddo
      enddo

      if (mfdim(1).eq.1) then
      do k=1,kx
      do j=1,jx
      do i=2,ix-1
        cy(i,j,k) = cy(i,j,k)-(bz(i+1,j,k)-bz(i-1,j,k))/dx(i)/2.d0
        cz(i,j,k) = cz(i,j,k)+(by(i+1,j,k)-by(i-1,j,k))/dx(i)/2.d0
      enddo
      enddo
      enddo
      endif

      if (mfdim(2).eq.1) then
      do k=1,kx
      do j=2,jx-1
      do i=1,ix
        cx(i,j,k) = cx(i,j,k)+(bz(i,j+1,k)-bz(i,j-1,k))/dy(j)/2.d0
        cz(i,j,k) = cz(i,j,k)-(bx(i,j+1,k)-bx(i,j-1,k))/dy(j)/2.d0
      enddo
      enddo
      enddo
      endif

      if (mfdim(3).eq.1) then
      do k=2,kx-1
      do j=1,jx
      do i=1,ix
        cx(i,j,k) = cx(i,j,k)-(by(i,j,k+1)-by(i,j,k-1))/dz(k)/2.d0
        cy(i,j,k) = cy(i,j,k)+(bx(i,j,k+1)-bx(i,j,k-1))/dz(k)/2.d0
      enddo
      enddo
      enddo
      endif

      return
      end
